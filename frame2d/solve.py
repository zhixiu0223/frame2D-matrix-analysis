"""
組裝 K, F -> 處理邊界條件 -> 解 Ku=F -> 回代桿端內力

均佈載重的等效節點載重規則(採 work-equivalent consistent load 定義,
用 f_FE = fixed_end_forces_udl() 算出的向量與局部y方向載重w同號):
  1. f_FE 是「work-equivalent 等效節點載重」本身(不是反力), 直接疊加進全域 F (不取負號)
  2. 解出位移後, 桿端真實內力 = k_local @ u_local - f_FE_local
  3. 已用「全固定端(u=0)退化情況」驗證: 此時桿端內力/反力 = -f_FE, 大小=wL/2, 方向與w相反,
     符合物理直覺(均佈載重方向下, 支承反力方向上)
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
)


@dataclass
class MemberResult:
    member_id: int
    L: float
    angle: float
    # 局部座標系桿端內力: [Fx1, Fy1, M1, Fx2, Fy2, M2]
    # Fx = 軸力(正=拉力), Fy = 剪力, M = 彎矩
    end_forces_local: np.ndarray


@dataclass
class SolveResult:
    displacements: np.ndarray          # 全域自由度位移向量, 長度 3*n_nodes
    reactions: np.ndarray              # 全域自由度反力向量 (僅支承自由度有意義)
    member_results: dict               # member_id -> MemberResult


def _support_dof_mask(frame: Frame2D) -> np.ndarray:
    """回傳長度n_dof的布林陣列, True=被支承拘束(位移=0)"""
    n = frame.n_dof()
    fixed = np.zeros(n, dtype=bool)
    for s in frame.supports:
        ux, uy, rot = frame.dofs_of(s.node)
        if s.ux:
            fixed[ux] = True
        if s.uy:
            fixed[uy] = True
        if s.rot:
            fixed[rot] = True
    return fixed


def solve(frame: Frame2D) -> SolveResult:
    n = frame.n_dof()
    K = np.zeros((n, n))
    F = np.zeros(n)

    # 每根桿件的固定端反力(局部座標), 之後回代內力要用
    fixed_end_local = {}

    # ---- 1. 組裝勁度矩陣 ----
    member_dofs = {}
    member_T = {}
    member_L = {}
    member_angle = {}
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        node_i = frame.nodes[m.node_i]
        node_j = frame.nodes[m.node_j]
        k_global, L, angle, T = member_stiffness_global(section, node_i, node_j, m.member_type)

        dofs = frame.dofs_of(m.node_i) + frame.dofs_of(m.node_j)  # 6個編號
        member_dofs[mid] = dofs
        member_T[mid] = T
        member_L[mid] = L
        member_angle[mid] = angle

        idx = np.array(dofs)
        K[np.ix_(idx, idx)] += k_global

    # ---- 2. 分佈載重 -> 固定端反力 -> 等效節點載重疊加進F ----
    for dl in frame.distributed_loads:
        m = frame.members[dl.member]
        if m.member_type == 'truss':
            raise ValueError(
                f"member {dl.member} 是桁架元素(truss), 兩端鉸接、沒有彎曲勁度,"
                " 不能承受垂直分佈載重(fixed_end_forces_udl假設的是有彎曲能力的"
                " frame元素)。如果要模擬桁架自重, 改成在兩端節點各加一半重量的"
                " point_load。")
        L = member_L[dl.member]
        T = member_T[dl.member]
        f_FE_local = fixed_end_forces_udl(dl.w_start, dl.w_end, L)
        fixed_end_local[dl.member] = fixed_end_local.get(
            dl.member, np.zeros(6)) + f_FE_local

        f_FE_global = T.T @ f_FE_local
        dofs = member_dofs[dl.member]
        idx = np.array(dofs)
        F[idx] += f_FE_global   # 等效節點載重(work-equivalent, 直接疊加, 不取負號)

    # ---- 3. 集中載重疊加進F ----
    for pl in frame.point_loads:
        ux, uy, rot = frame.dofs_of(pl.node)
        F[ux] += pl.fx
        F[uy] += pl.fy
        F[rot] += pl.m

    # ---- 4. 邊界條件: partition method (劃掉被拘束的自由度) ----
    fixed_mask = _support_dof_mask(frame)

    # 純桁架節點(只有truss桿件連接, 沒有任何frame桿件)的轉角自由度(rot),
    # 完全沒有任何桿件貢獻勁度(truss只傳軸力), 也沒有被支承拘束的話,
    # K的那一列/行會是全零, 直接丟進去解會是奇異矩陣。
    # 這種"沒有任何勁度、也沒有外加彎矩"的自由度物理上不作功, 直接跳過不解,
    # 位移設0即可(既不影響其他自由度的解, 也不會有反力)。
    diag = np.diag(K)
    inactive_mask = (~fixed_mask) & (np.isclose(diag, 0.0))
    if np.any(inactive_mask) and np.any(np.abs(F[inactive_mask]) > 1e-9):
        raise ValueError(
            "有自由度(通常是純桁架節點的轉角)完全沒有勁度貢獻, 卻被施加了外力"
            "(例如對一個只連接truss桿件的節點加彎矩)。這種自由度沒有任何元素"
            "可以抵抗該外力, 系統無法平衡, 請檢查模型。")

    free = np.where((~fixed_mask) & (~inactive_mask))[0]
    sup = np.where(fixed_mask)[0]

    K_ff = K[np.ix_(free, free)]
    F_f = F[free]

    u = np.zeros(n)
    if len(free) > 0:
        u_f = np.linalg.solve(K_ff, F_f)
        u[free] = u_f

    # ---- 5. 反力 = K @ u - F (僅支承自由度有意義) ----
    R = K @ u - F

    # ---- 6. 回代桿端內力 ----
    member_results = {}
    for mid, m in frame.members.items():
        dofs = member_dofs[mid]
        T = member_T[mid]
        u_member_global = u[np.array(dofs)]
        u_local = T @ u_member_global

        section = frame.sections[m.section]
        k_local = member_local_stiffness_dispatch(
            m.member_type, section.E, section.I, section.A, member_L[mid])

        f_FE = fixed_end_local.get(mid, np.zeros(6))
        end_forces_local = k_local @ u_local - f_FE

        member_results[mid] = MemberResult(
            member_id=mid,
            L=member_L[mid],
            angle=member_angle[mid],
            end_forces_local=end_forces_local,
        )

    return SolveResult(displacements=u, reactions=R, member_results=member_results)
