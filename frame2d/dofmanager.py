"""
DOFManager 版本的求解器 -- 內部鉸接的另一種實作方式, 用來跟solve.py的
靜力凝縮(static condensation)版本交叉驗證。

跟solve.py的關鍵差異: 這裡完全不用組裝前的靜力凝縮, 而是真的讓每個
release端擁有自己專屬、不跟別人共用的轉角自由度。因為DOF系統本身
「認得」哪裡有鉸接, 桿件的局部勁度矩陣跟固定端反力公式可以直接用最
標準的版本(跟完全剛接時一模一樣的公式), 不需要另外推導release專屬
的公式——鉸接處M=0這件事, 是解聯立方程組時自然浮現的結果(那個自由度
沒有其他任何東西連著, 平衡方程式自然會讓它转到不需要傳遞彎矩的角度)。

這在數學上跟靜力凝縮是等價的(高斯消去法, 只是消去的時機不同: 靜力
凝縮在組裝前手動消去, 這裡讓solver在解Ku=F的時候自然消去), 所以兩者
答案理論上要精確一致(浮點誤差等級)。這份檔案刻意獨立於solve.py寫,
不是包一層去呼叫它, 這樣交叉驗證抓到的才是「各自實作是否有bug」,
不是同一套邏輯自己驗自己。

**額外的好處**: 因為不需要「release專屬固定端反力公式」, DOFManager版本
可以直接支援solve.py目前還不支援的組合, 例如「release端的桿件上加
局部段均佈載重/桿件內部集中力」, 不需要為每種載重類型都推導一次
release專屬公式。見 tests/test_dofmanager_vs_condensation.py 的示範。
"""
from dataclasses import dataclass
import numpy as np

from .model import Frame2D
from .elements import (
    member_geometry,
    member_stiffness_local,
    member_stiffness_local_truss,
    transformation_matrix,
    fixed_end_forces_udl,
    fixed_end_forces_partial_udl,
    fixed_end_forces_point_load,
    fixed_end_forces_point_moment,
    fixed_end_forces_axial_point_load,
)


@dataclass
class DofManagerResult:
    displacements: dict          # global_dof_index -> 位移值 (只存有值的, 用dict因為DOF數量動態)
    reactions: np.ndarray        # 長度3*n_nodes, 只有節點自己的3個DOF(不含額外release DOF)
    member_end_forces: dict      # member_id -> [Fx1,Fy1,M1,Fx2,Fy2,M2] (局部座標)
    n_node_dof: int
    n_extra_dof: int


def build_dof_map(frame: Frame2D):
    """回傳 (member_dofs, n_node_dof, n_extra_dof)。
    member_dofs[mid] = (ux_i,uy_i,rot_i, ux_j,uy_j,rot_j) 六個全域DOF編號。
    節點自己的ux,uy,rot固定是3*node_id..+2(共用); 每個release端額外分配
    一個從3*n_nodes開始往後編號、不跟任何人共用的專屬轉角DOF。"""
    n_nodes = len(frame.nodes)
    n_node_dof = 3 * n_nodes
    next_extra = n_node_dof
    member_dofs = {}
    for mid, m in frame.members.items():
        ux_i, uy_i, rot_i = 3 * m.node_i, 3 * m.node_i + 1, 3 * m.node_i + 2
        ux_j, uy_j, rot_j = 3 * m.node_j, 3 * m.node_j + 1, 3 * m.node_j + 2
        if m.member_type == 'frame' and m.release_i:
            rot_i = next_extra
            next_extra += 1
        if m.member_type == 'frame' and m.release_j:
            rot_j = next_extra
            next_extra += 1
        member_dofs[mid] = (ux_i, uy_i, rot_i, ux_j, uy_j, rot_j)
    n_extra_dof = next_extra - n_node_dof
    return member_dofs, n_node_dof, n_extra_dof


def solve_dofmanager(frame: Frame2D) -> DofManagerResult:
    """跟solve.py的_solve_once()是完全獨立寫的第二套求解路徑, 用來交叉驗證
    release(內部鉸接)的靜力凝縮實作對不對。不支援cable鬆弛迭代(release
    交叉驗證用不到), 遇到cable直接當truss處理(如果有cable元素混在同一個
    模型裡, 請改用solve.py的主要求解器)。"""
    member_dofs, n_node_dof, n_extra_dof = build_dof_map(frame)
    n = n_node_dof + n_extra_dof
    K = np.zeros((n, n))
    F = np.zeros(n)

    member_T = {}
    member_L = {}
    fixed_end_local = {}

    # ---- 1. 組裝勁度矩陣: 全部都用「標準」局部勁度矩陣(不做靜力凝縮!),
    #      release的物理效果完全靠DOF不共用來體現 ----
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        node_i = frame.nodes[m.node_i]
        node_j = frame.nodes[m.node_j]
        L, angle = member_geometry(node_i, node_j)
        T = transformation_matrix(angle)
        member_T[mid] = T
        member_L[mid] = L

        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            k_local = member_stiffness_local(section.E, section.I, section.A, L)  # 標準版, release=False,False

        k_global = T.T @ k_local @ T
        idx = np.array(member_dofs[mid])
        K[np.ix_(idx, idx)] += k_global

    # ---- 2. 分佈載重: 全部都用「標準」固定端反力公式(不用release專屬公式!) ----
    for dl in frame.distributed_loads:
        m = frame.members[dl.member]
        L = member_L[dl.member]
        x_start = 0.0 if dl.x_start is None else dl.x_start
        x_end = L if dl.x_end is None else dl.x_end
        if x_start <= 1e-9 and x_end >= L - 1e-9:
            f_FE_local = fixed_end_forces_udl(dl.w_start, dl.w_end, L)
        else:
            f_FE_local = fixed_end_forces_partial_udl(dl.w_start, dl.w_end, x_start, x_end, L)
        fixed_end_local[dl.member] = fixed_end_local.get(dl.member, np.zeros(6)) + f_FE_local
        F[np.array(member_dofs[dl.member])] += member_T[dl.member].T @ f_FE_local

    # ---- 2b. 桿件內部集中力/力矩: 同樣用標準公式 ----
    for pl_m in frame.member_point_loads:
        L = member_L[pl_m.member]
        a = min(max(pl_m.a, 0.0), L)
        f_FE_local = np.zeros(6)
        if abs(pl_m.fx) > 0:
            f_FE_local += fixed_end_forces_axial_point_load(pl_m.fx, a, L)
        if abs(pl_m.fy) > 0:
            f_FE_local += fixed_end_forces_point_load(pl_m.fy, a, L)
        if abs(pl_m.m) > 0:
            f_FE_local += fixed_end_forces_point_moment(pl_m.m, a, L)
        fixed_end_local[pl_m.member] = fixed_end_local.get(pl_m.member, np.zeros(6)) + f_FE_local
        F[np.array(member_dofs[pl_m.member])] += member_T[pl_m.member].T @ f_FE_local

    # ---- 3. 節點集中力 ----
    for pl in frame.point_loads:
        F[3 * pl.node] += pl.fx
        F[3 * pl.node + 1] += pl.fy
        F[3 * pl.node + 2] += pl.m

    # ---- 4. 邊界條件(只作用在節點自己的DOF, 額外的release DOF一定是自由的) ----
    fixed_mask = np.zeros(n, dtype=bool)
    u_prescribed = np.zeros(n)
    for s in frame.supports:
        base = 3 * s.node
        for k, val in enumerate([s.ux, s.uy, s.rot]):
            if val is not None:
                fixed_mask[base + k] = True
                u_prescribed[base + k] = val

    free = np.where(~fixed_mask)[0]
    sup = np.where(fixed_mask)[0]

    # 跟solve.py主要求解器一樣的邊界情況: 某個自由度完全沒有任何桿件貢獻
    # 勁度(這裡常見於: 某節點的"共用"轉角DOF, 因為連接的桿件全部都用了
    # 自己專屬的release DOF, 沒有人用那個共用的了), 直接排除不解。
    diag = np.diag(K)
    inactive_mask = (~fixed_mask) & np.isclose(diag, 0.0)
    if np.any(inactive_mask) and np.any(np.abs(F[inactive_mask]) > 1e-9):
        raise ValueError("有自由度完全沒有勁度貢獻卻被施加外力, 系統無法平衡")
    free = np.where((~fixed_mask) & (~inactive_mask))[0]

    u = np.zeros(n)
    u[sup] = u_prescribed[sup]
    K_ff = K[np.ix_(free, free)]
    F_f = F[free]
    if len(sup) > 0:
        F_f = F_f - K[np.ix_(free, sup)] @ u[sup]
    u[free] = np.linalg.solve(K_ff, F_f)

    R = K @ u - F

    # ---- 5. 桿端內力回代 ----
    member_end_forces = {}
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        L = member_L[mid]
        T = member_T[mid]
        u_local = T @ u[np.array(member_dofs[mid])]
        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            k_local = member_stiffness_local(section.E, section.I, section.A, L)
        f_FE = fixed_end_local.get(mid, np.zeros(6))
        member_end_forces[mid] = k_local @ u_local - f_FE

    return DofManagerResult(
        displacements={i: u[i] for i in range(n)},
        reactions=R[:n_node_dof],
        member_end_forces=member_end_forces,
        n_node_dof=n_node_dof,
        n_extra_dof=n_extra_dof,
    )
