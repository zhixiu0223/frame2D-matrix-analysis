"""
DOFManager 版本的求解器 -- 內部鉸接的另一種實作方式, 現在是主要求解器
(見 frame2d/__init__.py: `solve` 這個公開名字指向這裡的 solve_dofmanager)。

跟 solve.py(靜力凝縮版本, 現在改名叫 solve_condensation, 保留當參考/回歸
測試實作)的關鍵差異: 這裡完全不用組裝前的靜力凝縮, 而是真的讓每個
release端擁有自己專屬、不跟別人共用的轉角自由度。因為DOF系統本身
「認得」哪裡有鉸接, 桿件的局部勁度矩陣跟固定端反力公式可以直接用最
標準的版本(跟完全剛接時一模一樣的公式), 不需要另外推導release專屬
的公式——鉸接處M=0這件事, 是解聯立方程組時自然浮現的結果(那個自由度
沒有其他任何東西連著, 平衡方程式自然會讓它转到不需要傳遞彎矩的角度)。

這在數學上跟靜力凝縮是等價的(高斯消去法, 只是消去的時機不同: 靜力
凝縮在組裝前手動消去, 這裡讓solver在解Ku=F的時候自然消去), 所以兩者
答案理論上要精確一致(浮點誤差等級)。這份檔案刻意獨立於solve.py寫,
不是包一層去呼叫它, 這樣交叉驗證抓到的才是「各自實作是否有bug」,
不是同一套邏輯自己驗自己(見 tests/test_dofmanager_vs_condensation.py)。

**額外的好處**: 因為不需要「release專屬固定端反力公式」, 這個版本可以
直接支援solve_condensation()還不支援的組合, 例如「release端的桿件上加
局部段均佈載重/桿件內部集中力」, 不需要為每種載重類型都推導一次
release專屬公式。

cable鬆弛迭代邏輯(結構上比照solve.py的模式, 底層勁度組裝維持獨立寫):
移除受壓cable的勁度貢獻 -> 重新求解 -> 檢查還有沒有cable受壓 -> 重複,
直到收斂或超過max_iterations。
"""
import numpy as np

from .model import Frame2D
from .result import MemberResult, SolveResult
from .elements import (
    member_geometry,
    member_stiffness_local,
    member_stiffness_local_truss,
    transformation_matrix,
    fixed_end_forces_udl,
    fixed_end_forces_axial_udl,
    fixed_end_forces_axial_udl_varying,
    fixed_end_forces_partial_udl,
    fixed_end_forces_axial_partial_udl,
    fixed_end_forces_point_load,
    fixed_end_forces_point_moment,
    fixed_end_forces_axial_point_load,
)


def build_dof_map(frame: Frame2D):
    """回傳 (member_dofs, n_node_dof, n_extra_dof)。
    member_dofs[mid] = (ux_i,uy_i,rot_i, ux_j,uy_j,rot_j) 六個全域DOF編號。
    節點自己的ux,uy,rot透過frame.dofs_of()取得(共用, 支援不連續/不從0開始
    的node_id); 每個release端額外分配一個從3*n_nodes開始往後編號、不跟
    任何人共用的專屬轉角DOF。"""
    n_nodes = len(frame.nodes)
    n_node_dof = 3 * n_nodes
    next_extra = n_node_dof
    member_dofs = {}
    for mid, m in frame.members.items():
        ux_i, uy_i, rot_i = frame.dofs_of(m.node_i)
        ux_j, uy_j, rot_j = frame.dofs_of(m.node_j)
        if m.member_type == 'frame' and m.release_i:
            rot_i = next_extra
            next_extra += 1
        if m.member_type == 'frame' and m.release_j:
            rot_j = next_extra
            next_extra += 1
        member_dofs[mid] = (ux_i, uy_i, rot_i, ux_j, uy_j, rot_j)
    n_extra_dof = next_extra - n_node_dof
    return member_dofs, n_node_dof, n_extra_dof


def _solve_once_dofmanager(frame: Frame2D, slack_cables: set) -> SolveResult:
    """跑一次線性求解, slack_cables裡的cable member直接跳過(不貢獻勁度,
    end_forces_local設為0向量), 跟solve.py的_solve_once()是完全獨立寫的
    第二套求解路徑。"""
    member_dofs, n_node_dof, n_extra_dof = build_dof_map(frame)
    n = n_node_dof + n_extra_dof
    K = np.zeros((n, n))
    F = np.zeros(n)

    member_T = {}
    member_L = {}
    fixed_end_local = {}

    # ---- 1. 組裝勁度矩陣: 全部都用「標準」局部勁度矩陣(不做靜力凝縮!),
    #      release的物理效果完全靠DOF不共用來體現; 鬆弛的cable跳過不組裝 ----
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        node_i = frame.nodes[m.node_i]
        node_j = frame.nodes[m.node_j]
        L, angle = member_geometry(node_i, node_j)
        T = transformation_matrix(angle)
        member_T[mid] = T
        member_L[mid] = L

        if mid in slack_cables:
            continue   # 鬆弛的cable: 不貢獻勁度, 視同不存在

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
        if m.member_type in ('truss', 'cable'):
            raise ValueError(
                f"member {dl.member} 是{m.member_type}元素, 兩端鉸接、沒有彎曲勁度,"
                " 不能承受分佈載重(fixed_end_forces_udl假設的是有彎曲能力的frame"
                " 元素; truss/cable自己的桿端彎矩M1/M2永遠必須是0, 加了均佈載重會讓"
                " 固定端彎矩項被靜默塞進求解, 得到錯誤但看起來合理的答案——這是"
                " 2026-09開發過程中發現的真實bug, 已修正為明確報錯)。如果要模擬"
                " 自重, 改成在兩端節點各加一半重量的point_load。")
        if dl.direction == 'global_y':
            # 全域垂直方向均佈/線性變化載重(大小以沿桿件長度量測, 例如
            # 屋頂重力/雪載重的標準表示方式): 依桿件角度分解成局部x(軸向)
            # +局部y(橫向)兩個分量的固定端反力, 直接相加。推導見elements.py
            # 的fixed_end_forces_axial_udl()說明, 已用slop-roof案例對照
            # SW FEA驗證過, 見tests/test_sloped_roof_global_udl.py。
            # 全域載重向量(0,-w)旋轉到局部座標: 直接用該桿件的T矩陣左上角
            # 2x2(平移自由度的旋轉部分)做向量旋轉, 不用另外手動算cos/sin。
            # w_start/w_end可以不同(非均勻/線性變化, 例如不均勻雪載重):
            # 桿件是直的, 角度沿桿長不變, 所以兩端各自投影到局部座標後,
            # 局部x、局部y分量沿桿長仍然各自是線性函數, 分開分解、分開
            # 代入線性變化公式即可(見model.py的DistributedLoad說明)。
            T_mat = member_T[dl.member]
            R = T_mat[0:2, 0:2]
            local_start = R @ np.array([0.0, -dl.w_start])
            local_end = R @ np.array([0.0, -dl.w_end])
            wx_start, wy_start = local_start[0], local_start[1]
            wx_end, wy_end = local_end[0], local_end[1]
            f_FE_local = (fixed_end_forces_udl(wy_start, wy_end, L)
                          + fixed_end_forces_axial_udl_varying(wx_start, wx_end, L))
        elif dl.direction == 'global':
            # 全域"任意角度"均佈載重(global_y的推廣版, 見model.py的
            # DistributedLoad說明): 用angle_deg指定角度(標準數學慣例),
            # 支援局部段(x_start/x_end)+線性變化的任意組合, 用跟
            # fixed_end_forces_partial_udl()同一套高斯積分手法分開處理
            # 局部x(fixed_end_forces_axial_partial_udl)、局部y(既有
            # fixed_end_forces_partial_udl)。
            x_start = 0.0 if dl.x_start is None else dl.x_start
            x_end = L if dl.x_end is None else dl.x_end
            ang = np.radians(dl.angle_deg)
            u_global = np.array([np.cos(ang), np.sin(ang)])
            T_mat = member_T[dl.member]
            R = T_mat[0:2, 0:2]
            local_start = R @ (u_global * dl.w_start)
            local_end = R @ (u_global * dl.w_end)
            wx_start, wy_start = local_start[0], local_start[1]
            wx_end, wy_end = local_end[0], local_end[1]
            f_FE_local = (fixed_end_forces_partial_udl(wy_start, wy_end, x_start, x_end, L)
                          + fixed_end_forces_axial_partial_udl(wx_start, wx_end, x_start, x_end, L))
        else:
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
        m = frame.members[pl_m.member]
        if m.member_type in ('truss', 'cable'):
            raise ValueError(
                f"member {pl_m.member} 是{m.member_type}元素, 兩端鉸接、沒有彎曲勁度,"
                " 不能承受桿件內部集中力/力矩(理由同distributed_load的檢查, 見"
                " 該處註解)。如果要模擬桁架/纜線中間的集中力, 改成拆成兩段桿件、"
                " 在新節點上用point_load。")
        L = member_L[pl_m.member]
        a = min(max(pl_m.a, 0.0), L)
        f_FE_local = np.zeros(6)
        if pl_m.direction == 'global':
            # 全域任意角度集中力(跟distributed_load的direction='global'
            # 同一套分解邏輯): 用該桿件的T矩陣把全域方向向量旋轉成局部分量
            ang = np.radians(pl_m.angle_deg)
            u_global = np.array([np.cos(ang), np.sin(ang)])
            T_mat = member_T[pl_m.member]
            R = T_mat[0:2, 0:2]
            local_vec = R @ (u_global * pl_m.F)
            fx_local, fy_local = local_vec[0], local_vec[1]
        else:
            fx_local, fy_local = pl_m.fx, pl_m.fy
        if abs(fx_local) > 0:
            f_FE_local += fixed_end_forces_axial_point_load(fx_local, a, L)
        if abs(fy_local) > 0:
            f_FE_local += fixed_end_forces_point_load(fy_local, a, L)
        if abs(pl_m.m) > 0:
            f_FE_local += fixed_end_forces_point_moment(pl_m.m, a, L)
        fixed_end_local[pl_m.member] = fixed_end_local.get(pl_m.member, np.zeros(6)) + f_FE_local
        F[np.array(member_dofs[pl_m.member])] += member_T[pl_m.member].T @ f_FE_local

    # ---- 3. 節點集中力 ----
    for pl in frame.point_loads:
        ux, uy, rot = frame.dofs_of(pl.node)
        F[ux] += pl.fx
        F[uy] += pl.fy
        F[rot] += pl.m

    # ---- 4. 邊界條件(只作用在節點自己的DOF, 額外的release DOF一定是自由的) ----
    fixed_mask = np.zeros(n, dtype=bool)
    u_prescribed = np.zeros(n)
    for s in frame.supports:
        ux, uy, rot = frame.dofs_of(s.node)
        for dof_idx, val in zip((ux, uy, rot), (s.ux, s.uy, s.rot)):
            if val is not None:
                fixed_mask[dof_idx] = True
                u_prescribed[dof_idx] = val

    # 某個自由度完全沒有任何桿件貢獻勁度(常見於: 純桁架/纜線節點的轉角,
    # 或某節點的"共用"轉角DOF因為連接的桿件全部都用了自己專屬的release
    # DOF, 或該DOF對應的桿件剛好全部鬆弛), 直接排除不解。
    diag = np.diag(K)
    inactive_mask = (~fixed_mask) & np.isclose(diag, 0.0)
    if np.any(inactive_mask) and np.any(np.abs(F[inactive_mask]) > 1e-9):
        raise ValueError("有自由度完全沒有勁度貢獻卻被施加外力, 系統無法平衡")
    free = np.where((~fixed_mask) & (~inactive_mask))[0]
    sup = np.where(fixed_mask)[0]

    u = np.zeros(n)
    u[sup] = u_prescribed[sup]
    K_ff = K[np.ix_(free, free)]
    F_f = F[free]
    if len(sup) > 0:
        F_f = F_f - K[np.ix_(free, sup)] @ u[sup]
    if len(free) > 0:
        try:
            u[free] = np.linalg.solve(K_ff, F_f)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                "勁度矩陣奇異, 無法求解 -- 常見原因: 結構是機構(支承/桿件不足以"
                "抵抗所有可能的位移模式), 或所有能穩住某個節點的cable都鬆弛了。"
                "請檢查模型的支承跟桿件連接。")

    R = K @ u - F

    # ---- 5. 桿端內力回代 ----
    member_results = {}
    for mid, m in frame.members.items():
        L, angle = member_L[mid], member_geometry(frame.nodes[m.node_i], frame.nodes[m.node_j])[1]
        if mid in slack_cables:
            member_results[mid] = MemberResult(
                member_id=mid, L=L, angle=angle, end_forces_local=np.zeros(6), slack=True)
            continue
        section = frame.sections[m.section]
        T = member_T[mid]
        u_local = T @ u[np.array(member_dofs[mid])]
        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            k_local = member_stiffness_local(section.E, section.I, section.A, L)
        f_FE = fixed_end_local.get(mid, np.zeros(6))
        end_forces_local = k_local @ u_local - f_FE
        member_results[mid] = MemberResult(
            member_id=mid, L=L, angle=angle, end_forces_local=end_forces_local, slack=False)

    return SolveResult(
        displacements=u[:n_node_dof],
        reactions=R[:n_node_dof],
        member_results=member_results,
        slack_cables=set(slack_cables),
        frame=frame,
    )


def solve_dofmanager(frame: Frame2D, max_iterations: int = 20) -> SolveResult:
    """求解frame(DOFManager版本, 主要求解器)。如果模型裡有cable元素,
    會自動反覆偵測+移除受壓的cable(鬆弛退出作用)並重新求解, 直到收斂
    (沒有cable受壓)或達到max_iterations。跟solve_condensation()回傳同一種
    SolveResult格式, plotting.py/postprocess.py不用關心是哪一種算法算的。"""
    has_cable = any(m.member_type == 'cable' for m in frame.members.values())
    if not has_cable:
        return _solve_once_dofmanager(frame, slack_cables=set())

    slack_cables = set()
    for _ in range(max_iterations):
        try:
            result = _solve_once_dofmanager(frame, slack_cables)
        except (ValueError, RuntimeError) as e:
            if slack_cables:
                raise RuntimeError(
                    f"移除鬆弛的cable {sorted(slack_cables)} 之後, 結構無法承受"
                    f"這個載重(可能變成機構, 或有節點完全失去支撐)。原始錯誤: {e}")
            raise
        newly_slack = set()
        for mid, m in frame.members.items():
            if m.member_type != 'cable' or mid in slack_cables:
                continue
            N = -result.member_results[mid].end_forces_local[0]   # 拉力為正
            if N < -1e-9:
                newly_slack.add(mid)
        if not newly_slack:
            return result
        slack_cables |= newly_slack

    raise RuntimeError(
        f"cable鬆弛迭代超過 {max_iterations} 次仍未收斂 -- 可能是模型設計本身有問題"
        " (例如載重下所有cable都會鬆弛、結構變成機構), 請檢查模型。")
