"""
frame2d 的最小 Web API 層。

設計原則: solver(frame2d 套件本身)完全不改一行,這層只負責
「JSON <-> Frame2D / SolveResult」的轉換,以及提供一支陽春的 HTML 測試頁,
用來驗證這條 API 路徑通不通(不是最終的前處理/後處理 GUI)。
"""
import math
from pathlib import Path

from fastapi import FastAPI, HTTPException, Response
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from frame2d import Frame2D, solve
from frame2d.dofmanager import solve_pdelta
from frame2d.postprocess import member_internal_forces, member_deformed_shape

from .schemas import FrameIn, SolveOut, NodeResultOut, MemberResultOut
from .diagrams import build_diagrams_and_deformed, build_deformed_with_scale
from .storage import LocalFileStorage, InvalidNameError, NotFoundError
from .pdf_export import build_pdf_report, build_fbd_previews, build_fbd_images_archive
from .query_point import query_point

app = FastAPI(title="frame2d API", description="frame2d 2D 矩陣位移法 solver 的 JSON API 外殼")

STATIC_DIR = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# saved_models/ 放在 repo 根目錄(跟 webapi/ webapi_stdlib/ 平行),
# 兩個後端指向同一個路徑,存檔清單不會因為換後端而分裂成兩份。
storage = LocalFileStorage(Path(__file__).parent.parent / "saved_models")


@app.get("/")
def index():
    """陽春測試頁: 貼 JSON -> 按 Solve -> 看結果表格。"""
    return FileResponse(STATIC_DIR / "index.html")


def _build_frame(payload: FrameIn) -> Frame2D:
    """把 JSON payload 轉成 frame2d.Frame2D。逐項對應 model.py 的建模 API,
    不直接操作 dataclass,避免繞過 Frame2D 自己的建模邏輯(例如
    node_index_cache 失效機制)。"""
    f = Frame2D()
    for n in payload.nodes:
        f.add_node(n.id, n.x, n.y)
    for s in payload.sections:
        f.add_section(s.name, E=s.E, I=s.I, A=s.A)
    for m in payload.members:
        f.add_member(m.id, node_i=m.node_i, node_j=m.node_j, section=m.section,
                     member_type=m.member_type, release_i=m.release_i, release_j=m.release_j)
    for sp in payload.supports:
        f.support(sp.node, ux=sp.ux, uy=sp.uy, rot=sp.rot)
    for pl in payload.point_loads:
        f.point_load(pl.node, fx=pl.fx, fy=pl.fy, m=pl.m)
    for dl in payload.distributed_loads:
        f.distributed_load(dl.member, w=dl.w_start, w_end=dl.w_end,
                            x_start=dl.x_start, x_end=dl.x_end,
                            direction=dl.direction, angle_deg=dl.angle_deg)
    for mpl in payload.member_point_loads:
        f.member_point_load(mpl.member, a=mpl.a, fx=mpl.fx, fy=mpl.fy, m=mpl.m,
                             direction=mpl.direction, F=mpl.F, angle_deg=mpl.angle_deg)
    return f


def _build_and_solve(payload: FrameIn, analysis_type: str = None):
    """_build_frame()+solve() 的共用包裝, 統一處理已知會從solve()冒出來
    的錯誤:
    1. KeyError: 「模型內有殘留的參照」(例如分割/刪除桿件後, 還留著
       指向舊桿件編號的均佈載重), frame2d底層會丟出KeyError(9)這種
       只印一個數字的錯誤, 前端顯示出來完全看不懂在講什麼, 這裡統一
       轉成看得懂的訊息。
    2. ValueError: dofmanager.py對「truss/cable桿件不能承受桿件內部
       載重(均佈載重/桿件集中力)」這個力學上的限制, 會丟出說明清楚
       的ValueError(2026-09修正, 見frame2d/dofmanager.py) —— 這個
       訊息本身已經寫得夠清楚, 直接原樣轉成400回傳, 不用像KeyError
       那樣另外翻譯。
    3. RuntimeError: solve_pdelta()疊代不收斂(通常是軸力接近挫屈載重),
       同樣直接把清楚的錯誤訊息轉成400。

    analysis_type: None(預設)時用payload.analysis_type; 呼叫端(例如
    /solve在analysis_type='pdelta'時額外算一次線性基準做比較)可以明確
    覆寫, 不用複製一份payload。"""
    if analysis_type is None:
        analysis_type = payload.analysis_type
    f = _build_frame(payload)
    try:
        if analysis_type == 'pdelta':
            result = solve_pdelta(f)
        else:
            result = solve(f)
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                   f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                   f"還留著指向舊桿件編號), 請檢查並移除",
        )
    except (ValueError, RuntimeError) as e:
        raise HTTPException(status_code=400, detail=str(e))
    return f, result


@app.post("/solve")
def solve_frame(payload: FrameIn):
    f, result = _build_and_solve(payload)

    node_out = []
    for n in payload.nodes:
        ux_i, uy_i, rot_i = f.dofs_of(n.id)
        node_out.append(NodeResultOut(
            node=n.id,
            ux=float(result.displacements[ux_i]),
            uy=float(result.displacements[uy_i]),
            rot=float(result.displacements[rot_i]),
            Rx=float(result.reactions[ux_i]),
            Ry=float(result.reactions[uy_i]),
            M=float(result.reactions[rot_i]),
        ))

    member_out = []
    for mid, mr in result.member_results.items():
        # 重用已驗證過的 postprocess.member_internal_forces 取兩端點的
        # N/V/M(已經是內部驗證過的正負號慣例:拉力為正),不在這層重新推導
        x, N, V, M = member_internal_forces(f, result, mid, n=2)
        member_out.append(MemberResultOut(
            member_id=mid,
            L=float(mr.L),
            angle_deg=float(math.degrees(mr.angle)),
            N1=float(N[0]), V1=float(V[0]), M1=float(M[0]),
            N2=float(N[-1]), V2=float(V[-1]), M2=float(M[-1]),
            slack=bool(mr.slack),
        ))

    diagrams, deformed, deform_scale = build_diagrams_and_deformed(f, result)

    solve_out = SolveOut(nodes=node_out, members=member_out)
    out = {
        **solve_out.model_dump(),
        "diagrams": diagrams,
        "deformed": deformed,
        "deform_scale": deform_scale,
        "analysis_type": payload.analysis_type,
    }

    if payload.analysis_type == 'pdelta':
        # 額外算一次純線性基準(P=0.0)做對照 -- 這裡不catch例外, 因為
        # 前面的_build_and_solve(payload)已經用同一個f、同一個'pdelta'
        # 跑過一次, 線性(P全部為0)理論上比pdelta更不容易發散/機構,
        # 如果連線性都解不出來, 表示模型本身就有問題, 讓例外往上冒出去
        # 用跟上面同一套HTTPException轉換邏輯處理即可。
        _, result_linear = _build_and_solve(payload, analysis_type='linear')
        max_disp_linear = result_linear.max_displacement()
        max_disp_pdelta = result.max_displacement()
        v_linear = max_disp_linear.value if max_disp_linear is not None else 0.0
        v_pdelta = max_disp_pdelta.value if max_disp_pdelta is not None else 0.0
        out["pdelta_comparison"] = {
            "max_disp_linear": float(v_linear),
            "max_disp_pdelta": float(v_pdelta),
            "amplification": float(v_pdelta / v_linear) if v_linear > 1e-12 else None,
            "iterations": result.pdelta_iterations,
        }
        # 疊圖用同一個deform_scale(來自pdelta結果, 通常較大)畫線性變形,
        # 這樣兩條線畫在同一張圖上, 大小差異才是真的P-Delta放大效果,
        # 不是各自套用自動縮放後看起來一樣大。
        out["deformed_linear"] = build_deformed_with_scale(f, result_linear, deform_scale)

    return out


# ---------------- 存檔 / 讀檔(Save / Save As / Load / 刪除) ----------------

@app.get("/models")
def list_models():
    return {"names": storage.list()}


@app.get("/models/{name}")
def load_model(name: str):
    try:
        return storage.load(name)
    except InvalidNameError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except NotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.post("/models/{name}")
def save_model(name: str, payload: FrameIn):
    try:
        storage.save(name, payload.model_dump())
    except InvalidNameError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"saved": name}


@app.delete("/models/{name}")
def delete_model(name: str):
    try:
        storage.delete(name)
    except InvalidNameError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except NotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    return {"deleted": name}


# ---------------- PDF 匯出(重用 plotting.py 的 plot_all) ----------------

@app.post("/export/pdf")
def export_pdf(payload: FrameIn):
    f = _build_frame(payload)
    try:
        pdf_bytes = build_pdf_report(f, units=payload.units, member_ids=payload.member_ids,
                                      include_member_diagrams=not payload.fbd_only)
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                   f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                   f"還留著指向舊桿件編號), 請檢查並移除",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={"Content-Disposition": 'attachment; filename="frame2d_report.pdf"'},
    )


@app.post("/export/fbd_images")
def export_fbd_images(payload: FrameIn):
    """匯出指定桿件的自由體圖圖檔(不是PDF報告): 只有1根桿件時直接
    回傳PNG, 多根桿件打包成zip(每根一個檔案), 給只想要圖片、不需要
    完整報告的情境用。"""
    f = _build_frame(payload)
    if not payload.member_ids:
        raise HTTPException(status_code=400, detail="沒有指定要匯出的桿件(member_ids)")
    try:
        zip_bytes = build_fbd_images_archive(f, payload.member_ids, units=payload.units)
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                   f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                   f"還留著指向舊桿件編號), 請檢查並移除",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    if len(payload.member_ids) == 1:
        # 只有1根桿件時, 直接把zip裡唯一那張PNG解開單獨回傳, 不用
        # 讓使用者多一道「下載zip再解壓縮才能看到一張圖」的手續。
        import io as _io
        import zipfile as _zipfile
        with _zipfile.ZipFile(_io.BytesIO(zip_bytes)) as zf:
            names = zf.namelist()
            if names:
                png_bytes = zf.read(names[0])
                return Response(
                    content=png_bytes, media_type="image/png",
                    headers={"Content-Disposition": f'attachment; filename="member_{payload.member_ids[0]}_fbd.png"'},
                )
    return Response(
        content=zip_bytes, media_type="application/zip",
        headers={"Content-Disposition": 'attachment; filename="frame2d_fbd_images.zip"'},
    )


@app.post("/preview/fbd")
def preview_fbd(payload: FrameIn):
    """匯出自由體圖PDF之前, 先讓前端秀出每根指定桿件的自由體圖預覽
    (base64 PNG), 使用者可以看過之後個別移除不想要的桿件, 確認後
    才真的呼叫/export/pdf。"""
    f = _build_frame(payload)
    if not payload.member_ids:
        raise HTTPException(status_code=400, detail="沒有指定要預覽的桿件(member_ids)")
    try:
        previews = build_fbd_previews(f, payload.member_ids, units=payload.units)
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                   f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                   f"還留著指向舊桿件編號), 請檢查並移除",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"previews": previews}


# ---------------- 查詢桿件內部任意位置的內力/位移 ----------------

class QueryPointIn(FrameIn):
    member: int
    mode: str  # 'relative' | 'absolute'
    value: float


@app.post("/query_point")
def query_point_endpoint(payload: QueryPointIn):
    f, result = _build_and_solve(payload)
    try:
        return query_point(f, result, payload.member, payload.mode, payload.value)
    except (KeyError, ValueError) as e:
        raise HTTPException(status_code=400, detail=str(e))
