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
from frame2d.postprocess import member_internal_forces, member_deformed_shape

from .schemas import FrameIn, SolveOut, NodeResultOut, MemberResultOut
from .diagrams import build_diagrams_and_deformed
from .storage import LocalFileStorage, InvalidNameError, NotFoundError
from .pdf_export import build_pdf_report
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


def _build_and_solve(payload: FrameIn):
    """_build_frame()+solve() 的共用包裝, 統一處理兩類已知會從solve()
    冒出來的錯誤:
    1. KeyError: 「模型內有殘留的參照」(例如分割/刪除桿件後, 還留著
       指向舊桿件編號的均佈載重), frame2d底層會丟出KeyError(9)這種
       只印一個數字的錯誤, 前端顯示出來完全看不懂在講什麼, 這裡統一
       轉成看得懂的訊息。
    2. ValueError: dofmanager.py對「truss/cable桿件不能承受桿件內部
       載重(均佈載重/桿件集中力)」這個力學上的限制, 會丟出說明清楚
       的ValueError(2026-09修正, 見frame2d/dofmanager.py) —— 這個
       訊息本身已經寫得夠清楚, 直接原樣轉成400回傳, 不用像KeyError
       那樣另外翻譯。"""
    f = _build_frame(payload)
    try:
        result = solve(f)
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                   f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                   f"還留著指向舊桿件編號), 請檢查並移除",
        )
    except ValueError as e:
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
    return {
        **solve_out.model_dump(),
        "diagrams": diagrams,
        "deformed": deformed,
        "deform_scale": deform_scale,
    }


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
        pdf_bytes = build_pdf_report(f, units=payload.units, member_ids=payload.member_ids)
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                   f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                   f"還留著指向舊桿件編號), 請檢查並移除",
        )
    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={"Content-Disposition": 'attachment; filename="frame2d_report.pdf"'},
    )


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
