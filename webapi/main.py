"""
frame2d 的最小 Web API 層。

設計原則: solver(frame2d 套件本身)完全不改一行,這層只負責
「JSON <-> Frame2D / SolveResult」的轉換,以及提供一支陽春的 HTML 測試頁,
用來驗證這條 API 路徑通不通(不是最終的前處理/後處理 GUI)。
"""
import math
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Response
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from frame2d import Frame2D, solve
from frame2d.dofmanager import solve_pdelta, initial_hinge_states
from frame2d.pushover import run_pushover, run_pushover_converged, apply_gravity
from frame2d.newton import run_pushover_newton
from frame2d.postprocess import member_internal_forces, member_deformed_shape

from .schemas import FrameIn, SolveOut, NodeResultOut, MemberResultOut
from .diagrams import build_diagrams_and_deformed, build_deformed_with_scale
from .storage import LocalFileStorage, InvalidNameError, NotFoundError
from .pdf_export import build_pdf_report, build_fbd_previews, build_fbd_images_archive, build_pushover_pdf_report, _pushover_final_solve_result
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
                     member_type=m.member_type, release_i=m.release_i, release_j=m.release_j,
                     Mp_i=m.Mp_i, Mp_j=m.Mp_j,
                     R_post_yield_i=m.R_post_yield_i, R_post_yield_j=m.R_post_yield_j)
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


def _prepare_pushover_run(payload: FrameIn):
    """/solve(pushover)跟/export/pdf(pushover)共用的前置邏輯: 檢核欄位、
    建frame、算控制點/底反力DOF、視情況跑重力預載階段。回傳一個dict,
    可以直接展開餵給run_pushover()(**kwargs), 呼叫端自己決定要不要
    額外加include_snapshots/include_final_displacement等旗標——這裡
    只負責"跑pushover需要的固定部分", 不管呼叫端想要哪些額外的回傳值。
    """
    control_nodes = payload.pushover_control_nodes
    if control_nodes is None:
        if payload.pushover_control_node is None:
            raise HTTPException(
                status_code=400,
                detail="pushover需要指定pushover_control_node(單點)或"
                       "pushover_control_nodes(多點)其中一個",
            )
        control_nodes = [payload.pushover_control_node]
    weights = payload.pushover_weights
    if weights is None:
        weights = [1.0] * len(control_nodes)
    elif len(weights) != len(control_nodes):
        raise HTTPException(
            status_code=400,
            detail=f"pushover_weights長度({len(weights)})必須跟"
                   f"pushover_control_nodes長度({len(control_nodes)})一樣",
        )
    if payload.pushover_target is None or payload.pushover_step is None:
        raise HTTPException(
            status_code=400,
            detail="pushover需要指定pushover_target、pushover_step兩個欄位",
        )

    f = _build_frame(payload)
    hinge_states = initial_hinge_states(f)
    if not hinge_states:
        raise HTTPException(
            status_code=400,
            detail="模型裡沒有任何桿件設定Mp_i(塑鉸容量), 無法做pushover -- "
                   "至少要有一根frame桿件設定Mp_i/Mp_j/R_post_yield_i/R_post_yield_j",
        )

    local_idx = {'x': 0, 'y': 1}[payload.pushover_direction]
    try:
        control_dofs = [f.dofs_of(n)[local_idx] for n in control_nodes]
    except KeyError as e:
        raise HTTPException(
            status_code=400,
            detail=f"找不到pushover_control_node(s)指定的節點id {e}",
        )
    base_reaction_dofs = list(dict.fromkeys(
        f.dofs_of(s.node)[local_idx] for s in f.supports))

    initial_cum_forces = None
    has_gravity_loads = bool(f.point_loads or f.distributed_loads or f.member_point_loads)
    if has_gravity_loads:
        try:
            initial_cum_forces, _ = apply_gravity(f, hinge_states)
        except (ValueError, KeyError, RuntimeError) as e:
            raise HTTPException(status_code=400, detail=f"重力預載階段求解失敗: {e}")

    return f, dict(
        frame=f, hinge_states=hinge_states, prescribed_dofs=control_dofs, direction=weights,
        target_total=payload.pushover_target, d_nominal=payload.pushover_step,
        base_reaction_dofs=base_reaction_dofs, initial_cum_forces=initial_cum_forces,
        use_pdelta=payload.pushover_use_pdelta,
        mechanism_ratio_limit=payload.pushover_mechanism_ratio_limit,
        control_mode=payload.pushover_control_mode,
        geometry_update=payload.pushover_geometry_update,
    )


def _run_selected_pushover_solver(payload: FrameIn, run_kwargs: dict, **extra_flags):
    """依payload.pushover_solver呼叫run_pushover()(event-to-event,
    預設)、run_pushover_converged()(幾何平衡疊代版本)、或
    run_pushover_newton()(真正的co-rotational+Newton-Raphson版本)。
    前兩者共用_prepare_pushover_run()準備好的同一組run_kwargs(參數
    名稱完全對得上兩個函式的簽名), 只有converged版本會額外用到
    geom_tol/max_geom_iter。

    newton版本的參數簽名不一樣(沒有use_pdelta/geometry_update/
    mechanism_ratio_limit這些概念——它是完全獨立的另一套實作, 見
    frame2d.newton的說明), 這裡只挑它認得的欄位轉傳; 重力/桿件內部
    載重(distributed_loads/member_point_loads)由run_pushover_newton()
    自己內部處理(見_gravity_fixed_end_forces()), 不需要
    _prepare_pushover_run()算好的initial_cum_forces, 這裡就不傳
    那個欄位過去。

    newton版本回傳的第5個值是converged(True=成功, 語意跟前兩者的
    mechanism_reached(True=失敗提前停止)剛好相反)——這裡統一轉成
    跟前兩者一致的mechanism_reached語意, 讓呼叫端不用另外分支處理。

    extra_flags是呼叫端想額外開的旗標(例如include_snapshots=True),
    三條路徑都會收到(newton版本不認得include_final_reactions/
    include_max_rotation, 呼叫端如果對newton solver要這兩個旗標,
    這裡會直接反映成TypeError, 是刻意讓錯誤明顯冒出來, 不是靜默丟棄)。
    """
    if payload.pushover_solver == 'converged':
        return run_pushover_converged(
            geom_tol=payload.pushover_geom_tol,
            max_geom_iter=payload.pushover_max_geom_iter,
            **run_kwargs, **extra_flags,
        )
    if payload.pushover_solver == 'newton':
        newton_kwargs = {
            k: run_kwargs[k] for k in
            ('frame', 'hinge_states', 'prescribed_dofs', 'direction', 'target_total',
             'd_nominal', 'base_reaction_dofs', 'control_mode')
        }
        raw = run_pushover_newton(
            tol=payload.pushover_newton_tol,
            max_iter=payload.pushover_newton_max_iter,
            **newton_kwargs, **extra_flags,
        )
        # raw[4]是converged(True=成功), 這裡轉成mechanism_reached語意
        # (True=失敗提前停止), 跟run_pushover()/run_pushover_converged()
        # 一致, 其餘欄位原封不動照順序傳回去。
        return raw[:4] + (not raw[4],) + raw[5:]
    return run_pushover(**run_kwargs, **extra_flags)


def _build_hinge_summary(hs_final):
    """把{member_id: HingeState}轉成JSON安全的摘要dict(Mp的inf轉成None,
    見infinity-json-bugfix那個修正)。/solve(pushover)跟/export/pdf
    (pushover)都要用到同一份塑鉸狀態摘要, 抽出來共用, 不要兩邊分別維護
    同一段轉換邏輯。"""
    return {
        str(mid): {
            "Mp": [None if math.isinf(hs.Mp[0]) else float(hs.Mp[0]),
                   None if math.isinf(hs.Mp[1]) else float(hs.Mp[1])],
            "yielded": [bool(hs.yielded[0]), bool(hs.yielded[1])],
            "theta_p": [float(hs.theta_p[0]), float(hs.theta_p[1])],
            "performance_level": [hs.performance_level(0), hs.performance_level(1)],
        }
        for mid, hs in hs_final.items()
    }


def _build_event_log_out(event_log):
    """把run_pushover()回傳的event_log(內部list of dict, 值可能是
    numpy純量)轉成JSON安全的list of dict。"""
    return [
        {"u": float(ev['u']), "F": float(ev['F']),
         "yielded": [[mid, end_idx] for mid, end_idx in ev['yielded']]}
        for ev in event_log
    ]


def _solve_pushover(payload: FrameIn):
    """analysis_type='pushover'的獨立處理路徑, 不跟_build_and_solve()共用
    (那個函式的dispatch邏輯是linear/pdelta二選一, pushover是完全不同的
    位移控制+event-to-event流程, 混在一起會讓兩邊都難懂)。

    回傳跟linear/pdelta的/solve完全不同的形狀(沒有nodes/members/diagrams/
    deformed, 因為那些是"單一狀態的解", pushover的重點是"一整條歷程"):
    容量曲線(history_u/history_F)、每個降伏事件的清單、最終每個塑鉸的
    狀態摘要(Mp、降伏與否、累積塑性轉角、IO/LS/CP分類)、有沒有形成機構。
    """
    f, run_kwargs = _prepare_pushover_run(payload)
    try:
        history_u, history_F, event_log, hs_final, mechanism, snapshots = _run_selected_pushover_solver(
            payload, run_kwargs, include_snapshots=True)
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e))

    hinge_summary = _build_hinge_summary(hs_final)
    event_log_out = _build_event_log_out(event_log)

    # 逐步回放用: 每一步(跟history_u/history_F逐一對應, 已經是純Python
    # list/float, 不含inf)的桿件端點力+塑鉸狀態。member_forces的key
    # 統一轉成字串(JSON物件的key本來就只能是字串, 跟hinge_summary同一個
    #處理方式)。
    snapshots_out = [
        {
            "member_forces": {str(mid): f for mid, f in snap['member_forces'].items()},
            "hinge_states": {str(mid): hs for mid, hs in snap['hinge_states'].items()},
            **({"u_full": snap['u_full']} if 'u_full' in snap else {}),
        }
        for snap in snapshots
    ]

    return {
        "analysis_type": "pushover",
        "history_u": [float(v) for v in history_u],
        "history_F": [float(v) for v in history_F],
        "event_log": event_log_out,
        "mechanism_reached": bool(mechanism),
        "hinge_summary": hinge_summary,
        "history_snapshots": snapshots_out,
        "solver": payload.pushover_solver,
    }


@app.post("/solve")
def solve_frame(payload: FrameIn):
    if payload.analysis_type == 'pushover':
        return _solve_pushover(payload)

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


@app.post("/pushover_step_diagrams")
def pushover_step_diagrams(payload: FrameIn):
    """Pushover逐步回放要看某一步精確N/V/M圖+真正變形後形狀時用:
    不重新跑一次側推, 直接拿前端已經有的那一步快照(pushover_step_
    member_forces/pushover_step_u_full, 來自/solve pushover回應裡
    history_snapshots[i]的原始內容)組一個假的SolveResult, 餵給跟
    線性分析同一套build_diagrams_and_deformed()——這樣算出來的N/V/M
    分佈跟變形曲線, 是用postprocess.py裡已經驗證過的精確積分方法,
    不是逐步回放原本那種"只用桿件端點值線性內插"的簡化畫法, 對桿件
    內部有分布載重的情況才會精確。

    回傳形狀刻意跟/solve(線性)的回應一致(nodes/members/diagrams/
    deformed/deform_scale), 前端可以直接重用既有的組合顯示/N/V/M/
    變形圖繪製邏輯, 不用另外寫一套。

    反力(Rx/Ry/M)這裡固定回傳0——因為快照只存了桿件內力跟位移場,
    沒有存每一步的完整反力向量, 而反力不影響diagrams/deformed的計算
    (只有節點結果表格會顯示反力, 這裡就顯示不出來, 這是刻意的取捨:
    不想為了這個而讓每一步快照的資料量再變大)。
    """
    if payload.pushover_step_member_forces is None or payload.pushover_step_u_full is None:
        raise HTTPException(
            status_code=400,
            detail="pushover_step_diagrams需要pushover_step_member_forces跟"
                   "pushover_step_u_full兩個欄位(從history_snapshots某一步"
                   "的原始內容直接傳過來)",
        )
    f = _build_frame(payload)
    try:
        cum_forces = {int(mid): np.array(vals, dtype=float)
                      for mid, vals in payload.pushover_step_member_forces.items()}
        u_full = np.array(payload.pushover_step_u_full, dtype=float)
    except (ValueError, TypeError) as e:
        raise HTTPException(status_code=400, detail=f"pushover_step資料格式錯誤: {e}")

    missing = set(f.members.keys()) - set(cum_forces.keys())
    if missing:
        raise HTTPException(
            status_code=400,
            detail=f"pushover_step_member_forces缺少桿件{sorted(missing)}的資料, "
                   f"模型跟快照對不起來(常見情況: 傳錯了/solve送出當時的model)",
        )

    result = _pushover_final_solve_result(f, cum_forces, u_full, np.zeros_like(u_full))
    diagrams, deformed, deform_scale = build_diagrams_and_deformed(f, result)

    node_out = []
    for n in payload.nodes:
        ux_i, uy_i, rot_i = f.dofs_of(n.id)
        node_out.append({
            "node": n.id, "ux": float(u_full[ux_i]), "uy": float(u_full[uy_i]), "rot": float(u_full[rot_i]),
            "Rx": 0.0, "Ry": 0.0, "M": 0.0,
        })
    member_out = []
    for mid, mr in result.member_results.items():
        x, N, V, M = member_internal_forces(f, result, mid, n=2)
        member_out.append({
            "member_id": mid, "L": float(mr.L), "angle_deg": float(math.degrees(mr.angle)),
            "N1": float(N[0]), "V1": float(V[0]), "M1": float(M[0]),
            "N2": float(N[-1]), "V2": float(V[-1]), "M2": float(M[-1]),
            "slack": bool(mr.slack),
        })

    return {
        "nodes": node_out, "members": member_out,
        "diagrams": diagrams, "deformed": deformed, "deform_scale": deform_scale,
        "analysis_type": "pushover_step",
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

def _export_pushover_pdf(payload: FrameIn) -> bytes:
    """analysis_type='pushover'時/export/pdf走的路徑: 重新跑一次pushover
    (拿到最終狀態需要的內力/位移/反力, 這幾個/solve的JSON回應裡沒有
    完整保留——history_snapshots雖然有每一步的內力, 但沒有最終的完整
    反力向量跟絕對位移向量, 這裡額外開include_final_displacement/
    include_final_reactions/include_max_rotation三個旗標拿到)。跟
    _solve_pushover()共用_prepare_pushover_run()這段前置邏輯, 不用
    兩邊分別維護。"""
    if payload.pushover_solver == 'newton':
        raise HTTPException(
            status_code=400,
            detail="PDF匯出目前還不支援newton求解器(它還沒有實作"
                   "include_final_reactions/include_max_rotation這兩個"
                   "PDF報告需要的欄位)——這是已知限制, 不是bug, 請改用"
                   "event_to_event或converged求解器匯出PDF, 或直接用"
                   "/solve查看newton求解器的結果。",
        )
    f, run_kwargs = _prepare_pushover_run(payload)
    try:
        (history_u, history_F, event_log, hs_final, mechanism,
         u_full_cum, snapshots, cum_reaction, max_rotation) = _run_selected_pushover_solver(
            payload, run_kwargs, include_final_displacement=True, include_snapshots=True,
            include_final_reactions=True, include_max_rotation=True)
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e))

    pushover_result = {
        "history_u": history_u, "history_F": history_F,
        "event_log": _build_event_log_out(event_log),
        "hinge_summary": _build_hinge_summary(hs_final),
        "mechanism_reached": bool(mechanism),
        "cum_forces": snapshots[-1]["member_forces"],
        "u_full_cum": u_full_cum,
        "cum_reaction": cum_reaction,
        "solver": payload.pushover_solver,
        "max_rotation": max_rotation,
    }
    return build_pushover_pdf_report(f, pushover_result, units=payload.units)


@app.post("/export/pdf")
def export_pdf(payload: FrameIn):
    if payload.analysis_type == 'pushover':
        try:
            pdf_bytes = _export_pushover_pdf(payload)
        except KeyError as e:
            raise HTTPException(
                status_code=400,
                detail=f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照, 請檢查並移除",
            )
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        return Response(
            content=pdf_bytes,
            media_type="application/pdf",
            headers={"Content-Disposition": 'attachment; filename="frame2d_pushover_report.pdf"'},
        )

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
