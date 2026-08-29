"""
組裝 K, F -> 處理邊界條件 -> 解 Ku=F -> 回代桿端內力

均佈載重的等效節點載重規則(採 work-equivalent consistent load 定義,
用 f_FE = fixed_end_forces_udl() 算出的向量與局部y方向載重w同號):
  1. f_FE 是「work-equivalent 等效節點載重」本身(不是反力), 直接疊加進全域 F (不取負號)
  2. 解出位移後, 桿端真實內力 = k_local @ u_local - f_FE_local
  3. 已用「全固定端(u=0)退化情況」驗證: 此時桿端內力/反力 = -f_FE, 大小=wL/2, 方向與w相反,
     符合物理直覺(均佈載重方向下, 支承反力方向上)

纜線(cable)只能受拉、不能受壓的迭代處理:
  cable桿件只能受拉力, 如果一次線性求解算出某條cable是壓力, 物理上代表這條纜線
  在這個載重下會鬆弛(slack)、退出作用, 不是真的能傳遞壓力。solve()因此採用
  「移除受壓cable的勁度貢獻 -> 重新求解 -> 檢查還有沒有cable受壓 -> 重複」的
  迭代流程(經典的tension-only member處理方式), 直到沒有cable受壓為止,
  或超過max_iterations(此時代表模型有問題, 例如全部cable都可能鬆弛導致機構,
  拋出RuntimeError)。純truss桿件不受此限制, 可以同時傳拉力或壓力。
"""
from dataclasses import dataclass
import numpy as np

from .model import Frame2D
from .elements import (
    member_geometry,
    member_stiffness_local,
    member_local_stiffness_dispatch,
    member_stiffness_global,
    transformation_matrix,
    fixed_end_forces_udl,
    fixed_end_forces_partial_udl,
    fixed_end_forces_udl_release_i,
    fixed_end_forces_udl_release_j,
    fixed_end_forces_point_load,
    fixed_end_forces_point_moment,
    fixed_end_forces_axial_point_load,
)


@dataclass
class MemberResult:
    member_id: int
    L: float
    angle: float
    # 局部座標系桿端內力: [Fx1, Fy1, M1, Fx2, Fy2, M2]
    # Fx = 軸力(正=拉力), Fy = 剪力, M = 彎矩
    end_forces_local: np.ndarray
    slack: bool = False   # True = 這是一條cable, 且這次求解判定為鬆弛(不受力)


@dataclass
class SolveResult:
    displacements: np.ndarray          # 全域自由度位移向量, 長度 3*n_nodes
    reactions: np.ndarray              # 全域自由度反力向量 (僅支承自由度有意義)
    member_results: dict               # member_id -> MemberResult
    slack_cables: set = None           # 最終解裡判定為鬆弛的cable member_id集合


def _support_dof_values(frame: Frame2D):
    """回傳(mask, values): mask是長度n_dof的布林陣列(True=被支承拘束,
    不管是fix的0.0還是強制位移的非0), values是同長度陣列, 在mask=True的
    位置填入指定的位移值。"""
    n = frame.n_dof()
    mask = np.zeros(n, dtype=bool)
    values = np.zeros(n)
    for s in frame.supports:
        ux, uy, rot = frame.dofs_of(s.node)
        if s.ux is not None:
            mask[ux] = True
            values[ux] = s.ux
        if s.uy is not None:
            mask[uy] = True
            values[uy] = s.uy
        if s.rot is not None:
            mask[rot] = True
            values[rot] = s.rot
    return mask, values


def _solve_once(frame: Frame2D, slack_cables: set) -> SolveResult:
    """跑一次線性求解, slack_cables裡的cable member直接跳過(不貢獻勁度,
    end_forces_local設為0向量)。"""
    n = frame.n_dof()
    K = np.zeros((n, n))
    F = np.zeros(n)

    fixed_end_local = {}
    member_dofs = {}
    member_T = {}
    member_L = {}
    member_angle = {}

    # ---- 1. 組裝勁度矩陣 (跳過slack_cables) ----
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        node_i = frame.nodes[m.node_i]
        node_j = frame.nodes[m.node_j]
        L, angle = member_geometry(node_i, node_j)
        T = transformation_matrix(angle)
        dofs = frame.dofs_of(m.node_i) + frame.dofs_of(m.node_j)
        member_dofs[mid] = dofs
        member_T[mid] = T
        member_L[mid] = L
        member_angle[mid] = angle

        if mid in slack_cables:
            continue   # 鬆弛的cable: 不貢獻勁度, 視同不存在

        k_global, _, _, _ = member_stiffness_global(
            section, node_i, node_j, m.member_type, m.release_i, m.release_j)
        idx = np.array(dofs)
        K[np.ix_(idx, idx)] += k_global

    # ---- 2. 分佈載重 -> 固定端反力 -> 等效節點載重疊加進F ----
    for dl in frame.distributed_loads:
        m = frame.members[dl.member]
        if m.member_type in ('truss', 'cable'):
            raise ValueError(
                f"member {dl.member} 是{m.member_type}元素, 兩端鉸接、沒有彎曲勁度,"
                " 不能承受垂直分佈載重(fixed_end_forces_udl假設的是有彎曲能力的"
                " frame元素)。如果要模擬自重, 改成在兩端節點各加一半重量的"
                " point_load。")
        if m.release_i or m.release_j:
            if dl.x_start is not None or dl.x_end is not None:
                raise ValueError(
                    f"member {dl.member} 有端點鉸接, 目前只支援整根桿件的均佈"
                    " 載重(不支援局部段), 且w_start必須等於w_end。")
            if m.release_i and m.release_j:
                raise ValueError(f"member {dl.member} 兩端都鉸接, 不支援分佈載重。")
            if m.release_j:
                f_FE_local = fixed_end_forces_udl_release_j(dl.w_start, dl.w_end, member_L[dl.member])
            else:
                f_FE_local = fixed_end_forces_udl_release_i(dl.w_start, dl.w_end, member_L[dl.member])
            fixed_end_local[dl.member] = fixed_end_local.get(
                dl.member, np.zeros(6)) + f_FE_local
            f_FE_global = member_T[dl.member].T @ f_FE_local
            F[np.array(member_dofs[dl.member])] += f_FE_global
            continue
        L = member_L[dl.member]
        T = member_T[dl.member]
        x_start = 0.0 if dl.x_start is None else dl.x_start
        x_end = L if dl.x_end is None else dl.x_end
        if x_start < -1e-9 or x_end > L + 1e-9 or x_start > x_end:
            raise ValueError(
                f"member {dl.member} 的分佈載重範圍 [{x_start},{x_end}] 超出"
                f"桿件長度[0,{L}]或起訖顛倒")
        if x_start <= 1e-9 and x_end >= L - 1e-9:
            f_FE_local = fixed_end_forces_udl(dl.w_start, dl.w_end, L)
        else:
            f_FE_local = fixed_end_forces_partial_udl(dl.w_start, dl.w_end, x_start, x_end, L)
        fixed_end_local[dl.member] = fixed_end_local.get(
            dl.member, np.zeros(6)) + f_FE_local

        f_FE_global = T.T @ f_FE_local
        idx = np.array(member_dofs[dl.member])
        F[idx] += f_FE_global

    # ---- 2b. 桿件內部集中力/力矩 -> 固定端反力 -> 等效節點載重疊加進F ----
    for pl_m in frame.member_point_loads:
        m = frame.members[pl_m.member]
        if m.member_type in ('truss', 'cable'):
            raise ValueError(
                f"member {pl_m.member} 是{m.member_type}元素, 兩端鉸接、沒有彎曲勁度,"
                " 不能承受橫向集中力/力矩。如果要模擬桁架/纜線中間的集中力,"
                " 改成拆成兩段桿件、在新節點上用point_load。")
        if m.release_i or m.release_j:
            raise ValueError(
                f"member {pl_m.member} 有端點鉸接(release_i/release_j), 目前沒有"
                " 處理鉸接端的固定端反力修正公式, 暫不支援在有鉸接的桿件上加"
                " 桿件內部集中力/力矩。")
        L = member_L[pl_m.member]
        T = member_T[pl_m.member]
        a = pl_m.a
        if not (0.0 <= a <= L + 1e-9):
            raise ValueError(f"member {pl_m.member} 的member_point_load位置a={a}超出桿件範圍[0,{L}]")
        a = min(max(a, 0.0), L)

        f_FE_local = np.zeros(6)
        if abs(pl_m.fx) > 0:
            f_FE_local += fixed_end_forces_axial_point_load(pl_m.fx, a, L)
        if abs(pl_m.fy) > 0:
            f_FE_local += fixed_end_forces_point_load(pl_m.fy, a, L)
        if abs(pl_m.m) > 0:
            f_FE_local += fixed_end_forces_point_moment(pl_m.m, a, L)

        fixed_end_local[pl_m.member] = fixed_end_local.get(
            pl_m.member, np.zeros(6)) + f_FE_local

        f_FE_global = T.T @ f_FE_local
        idx = np.array(member_dofs[pl_m.member])
        F[idx] += f_FE_global

    # ---- 3. 集中載重疊加進F ----
    for pl in frame.point_loads:
        ux, uy, rot = frame.dofs_of(pl.node)
        F[ux] += pl.fx
        F[uy] += pl.fy
        F[rot] += pl.m

    # ---- 4. 邊界條件: partition method, 支援非零指定位移(K_ff u_f = F_f - K_fc u_c) ----
    fixed_mask, u_prescribed = _support_dof_values(frame)

    # 純桁架/纜線節點(沒有任何frame桿件連接, 或連接的cable全部鬆弛)的轉角自由度,
    # 完全沒有任何桿件貢獻勁度, 也沒有被支承拘束的話, K的那一列/行會是全零,
    # 直接丟進去解會是奇異矩陣。這種自由度物理上不作功, 直接跳過不解即可。
    diag = np.diag(K)
    inactive_mask = (~fixed_mask) & (np.isclose(diag, 0.0))
    if np.any(inactive_mask) and np.any(np.abs(F[inactive_mask]) > 1e-9):
        raise ValueError(
            "有自由度(通常是純桁架/纜線節點的轉角, 或所有連接的cable都鬆弛的節點)"
            "完全沒有勁度貢獻, 卻被施加了外力。這種自由度沒有任何元素可以抵抗該外力,"
            "系統無法平衡, 請檢查模型。")

    free = np.where((~fixed_mask) & (~inactive_mask))[0]
    sup = np.where(fixed_mask)[0]

    u = np.zeros(n)
    u[sup] = u_prescribed[sup]   # 指定位移值(fix/pin/roller是0.0, 沉陷分析是非0)

    K_ff = K[np.ix_(free, free)]
    F_f = F[free]
    if len(sup) > 0:
        K_fc = K[np.ix_(free, sup)]
        F_f = F_f - K_fc @ u[sup]   # 已知位移對自由自由度的等效載重貢獻

    if len(free) > 0:
        try:
            u_f = np.linalg.solve(K_ff, F_f)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                "勁度矩陣奇異, 無法求解 -- 常見原因: 結構是機構(支承/桿件不足以"
                "抵抗所有可能的位移模式), 或所有能穩住某個節點的cable都鬆弛了。"
                "請檢查模型的支承跟桿件連接。")
        u[free] = u_f

    # ---- 5. 反力 = K @ u - F (僅支承自由度有意義) ----
    R = K @ u - F

    # ---- 6. 回代桿端內力 ----
    member_results = {}
    for mid, m in frame.members.items():
        if mid in slack_cables:
            member_results[mid] = MemberResult(
                member_id=mid, L=member_L[mid], angle=member_angle[mid],
                end_forces_local=np.zeros(6), slack=True)
            continue

        dofs = member_dofs[mid]
        T = member_T[mid]
        u_member_global = u[np.array(dofs)]
        u_local = T @ u_member_global

        section = frame.sections[m.section]
        k_local = member_local_stiffness_dispatch(
            m.member_type, section.E, section.I, section.A, member_L[mid], m.release_i, m.release_j)

        f_FE = fixed_end_local.get(mid, np.zeros(6))
        end_forces_local = k_local @ u_local - f_FE

        member_results[mid] = MemberResult(
            member_id=mid, L=member_L[mid], angle=member_angle[mid],
            end_forces_local=end_forces_local, slack=False)

    return SolveResult(displacements=u, reactions=R, member_results=member_results,
                        slack_cables=set(slack_cables))


def solve(frame: Frame2D, max_iterations: int = 20) -> SolveResult:
    """求解frame。如果模型裡有cable元素, 會自動反覆偵測+移除受壓的cable
    (鬆弛退出作用)並重新求解, 直到收斂(沒有cable受壓)或達到max_iterations。"""
    has_cable = any(m.member_type == 'cable' for m in frame.members.values())
    if not has_cable:
        return _solve_once(frame, slack_cables=set())

    slack_cables = set()
    for _ in range(max_iterations):
        try:
            result = _solve_once(frame, slack_cables)
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
