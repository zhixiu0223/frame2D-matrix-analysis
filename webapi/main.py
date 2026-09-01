"""
frame2d 的最小 Web API 層。

設計原則: solver(frame2d 套件本身)完全不改一行,這層只負責
「JSON <-> Frame2D / SolveResult」的轉換,以及提供一支陽春的 HTML 測試頁,
用來驗證這條 API 路徑通不通(不是最終的前處理/後處理 GUI)。
"""
import math
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

from .schemas import FrameIn, SolveOut, NodeResultOut, MemberResultOut

app = FastAPI(title="frame2d API", description="frame2d 2D 矩陣位移法 solver 的 JSON API 外殼")

STATIC_DIR = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


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


@app.post("/solve", response_model=SolveOut)
def solve_frame(payload: FrameIn) -> SolveOut:
    f = _build_frame(payload)
    result = solve(f)

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

    return SolveOut(nodes=node_out, members=member_out)
