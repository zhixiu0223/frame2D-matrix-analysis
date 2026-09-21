"""
質量矩陣 (動力分析 D1)。

`assemble_M(frame, kind)` 產生跟 `assembly.assemble_K()` **同一套DOF編號**
(節點DOF在前、release端專屬轉角DOF在後)的全域質量矩陣, 所以 K 跟 M 可以
直接放進同一個廣義特徵值問題 (K - ω²M)φ = 0。

單位(核心不做單位換算, 只要求一致)
----------------------------------
質量單位 = 力單位·s²/長度單位:

  SI(網頁後端固定用這套): 力N、長度m -> 質量kg, 密度kg/m³, 轉動慣量kg·m²
  kN、m手算/測試:        力kN、長度m -> 質量ton, 密度ton/m³, 轉動慣量ton·m²

單位寫錯不會報錯, 只會讓所有頻率差一個常數倍(例如把ton當kg, 頻率差√1000
倍)——所以驗證時一定要拿一個已知解析解的案例(例如頂端質量懸臂柱
ω²=3EI/(mL³))確認量級。

兩種質量矩陣
------------
kind='lumped'(集中質量, 預設):
  桿件質量 ρAL 各分一半到兩端節點的 ux, uy(全域方向); **不計轉動慣量**
  (轉角DOF質量為0)。節點質量 add_mass(mx, my, Iz) 直接加到對應DOF。
  矩陣是對角的。轉角DOF與release端專屬DOF沒有質量 -> M奇異, D2特徵值
  分析要先對無質量DOF做靜力凝縮(對無質量DOF是精確的)。

kind='consistent'(一致質量):
  Euler-Bernoulli 樑: 軸向用線性形函數、橫向用Hermite三次形函數, 局部6x6
  矩陣 (m̄ = ρA):
      軸向:  m̄L/6 · [[2,1],[1,2]]
      橫向:  m̄L/420 · [[156, 22L, 54, -13L],
                        [22L, 4L², 13L, -3L²],
                        [54, 13L, 156, -22L],
                        [-13L, -3L², -22L, 4L²]]
  再用 M_g = Tᵀ m T 轉到全域。不含轉動慣量項 ρI (Euler-Bernoulli 忽略
  剪切與轉動慣量, 跟勁度矩陣的假設一致)。桁架元素用 4x4 平動一致質量
  m̄L/6·[[2,0,1,0],[0,2,0,1],[1,0,2,0],[0,1,0,2]] (旋轉不變)。
  release端: 標準6x6直接放在該端的專屬DOF上, 這是一致有限元素對鉸接端的
  正確處理(那個DOF有自己的慣性項, 不是硬把端點彎矩當0)。

動力分析初期(D1~D5)的明確拒絕(見ROADMAP.md「既有限制」表格)
------------------------------------------------------------
- 含cable: 鬆弛是狀態相依的, 模態分析沒有定義好的基準狀態
- 含equal_dofs: 懲罰法的高勁度彈簧會製造虛假的高頻模態
- 非零指定支承位移: 動力分析只接受0.0/None; 地震輸入用等效力 -M·r·a_g
"""
import numpy as np

from .dofmanager import build_dof_map
from .elements import member_geometry, transformation_matrix
from .model import Frame2D

VALID_KINDS = ('lumped', 'consistent')


def check_dynamic_supported(frame: Frame2D):
    """動力分析的模型限制檢查, 遇到尚未支援的組合明確 raise ValueError。"""
    cables = sorted(mid for mid, m in frame.members.items() if m.member_type == 'cable')
    if cables:
        raise ValueError(
            f"動力分析目前不支援cable元素(member {cables}): cable的鬆弛是狀態相依的, "
            "模態分析需要先定義基準狀態(例如重力下拉緊)。請先把它們換成truss, "
            "或等後續階段支援。")
    if frame.equal_dofs:
        raise ValueError(
            "動力分析目前不支援equal_dof: 目前用懲罰法(高勁度彈簧)實現, 會在特徵值"
            "分析中製造虛假的高頻模態並讓矩陣條件數惡化。等改成精確的master-slave"
            "消去再支援。")
    for s in frame.supports:
        for name, val in (('ux', s.ux), ('uy', s.uy), ('rot', s.rot)):
            if val is not None and val != 0.0:
                raise ValueError(
                    f"節點{s.node}的支承{name}={val}是非零指定位移(沉陷/強制位移), 動力分析"
                    "只接受0.0或None。地震輸入請用等效力 -M·r·a_g, 不要用支承位移。")
    for nm in frame.node_masses:
        if nm.node not in frame.nodes:
            raise ValueError(f"add_mass()指定的節點{nm.node}不存在")


def consistent_mass_local(rho_A: float, L: float) -> np.ndarray:
    """Euler-Bernoulli 樑元素的局部6x6一致質量矩陣, DOF順序
    (u_i, v_i, θ_i, u_j, v_j, θ_j), rho_A = 單位長度質量。"""
    m = np.zeros((6, 6))
    a = rho_A * L / 6.0
    m[0, 0] = m[3, 3] = 2 * a
    m[0, 3] = m[3, 0] = a
    b = rho_A * L / 420.0
    bend = b * np.array([
        [156.0, 22 * L, 54.0, -13 * L],
        [22 * L, 4 * L * L, 13 * L, -3 * L * L],
        [54.0, 13 * L, 156.0, -22 * L],
        [-13 * L, -3 * L * L, -22 * L, 4 * L * L],
    ])
    idx = [1, 2, 4, 5]
    m[np.ix_(idx, idx)] = bend
    return m


def consistent_mass_local_truss(rho_A: float, L: float) -> np.ndarray:
    """桁架元素的局部6x6一致質量矩陣(只有平動DOF有值, 軸向與橫向慣性相同)。"""
    m = np.zeros((6, 6))
    a = rho_A * L / 6.0
    idx = [0, 1, 3, 4]
    m[np.ix_(idx, idx)] = a * np.array([
        [2.0, 0.0, 1.0, 0.0],
        [0.0, 2.0, 0.0, 1.0],
        [1.0, 0.0, 2.0, 0.0],
        [0.0, 1.0, 0.0, 2.0],
    ])
    return m


def assemble_M(frame: Frame2D, kind: str = 'lumped') -> np.ndarray:
    """組裝全域質量矩陣M, DOF編號跟 assembly.assemble_K() 完全相同。

    kind: 'lumped'(集中, 預設) 或 'consistent'(一致), 見模組說明。
    rho為None的斷面, 桿件沒有分佈質量(質量矩陣不含該桿件); 節點質量
    (add_mass)不受影響。
    """
    if kind not in VALID_KINDS:
        raise ValueError(f"kind必須是{VALID_KINDS}其中之一, 收到'{kind}'")
    check_dynamic_supported(frame)

    member_dofs, n_node_dof, n_extra_dof = build_dof_map(frame)
    n = n_node_dof + n_extra_dof
    M = np.zeros((n, n))

    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        if section.rho is None or section.rho == 0.0:
            continue
        L, angle = member_geometry(frame.nodes[m.node_i], frame.nodes[m.node_j])
        rho_A = section.rho * section.A
        idx = np.array(member_dofs[mid])
        if kind == 'lumped':
            half = rho_A * L / 2.0
            for dof in (idx[0], idx[1], idx[3], idx[4]):    # 兩端的ux, uy
                M[dof, dof] += half
        else:
            if m.member_type == 'truss':
                m_local = consistent_mass_local_truss(rho_A, L)
            else:
                m_local = consistent_mass_local(rho_A, L)
            T = transformation_matrix(angle)
            M[np.ix_(idx, idx)] += T.T @ m_local @ T

    for nm in frame.node_masses:
        ux, uy, rot = frame.dofs_of(nm.node)
        M[ux, ux] += nm.mx
        M[uy, uy] += nm.my
        M[rot, rot] += nm.Iz
    return M


def influence_vector(frame: Frame2D, direction: str = 'x') -> np.ndarray:
    """地面運動的影響向量 r: 該方向的節點平動DOF為1, 其餘(轉角、release
    專屬DOF、另一個方向)為0。長度跟 assemble_M() 的維度一致。"""
    if direction not in ('x', 'y'):
        raise ValueError(f"direction必須是'x'或'y', 收到'{direction}'")
    _, n_node_dof, n_extra_dof = build_dof_map(frame)
    r = np.zeros(n_node_dof + n_extra_dof)
    k = 0 if direction == 'x' else 1
    for nid in frame.nodes:
        r[frame.dofs_of(nid)[k]] = 1.0
    return r


def total_mass(frame: Frame2D, direction: str = 'x') -> float:
    """該方向的總質量 = 所有桿件分佈質量 ρAL + 節點質量(mx 或 my)。
    (質量矩陣的 rᵀ M r 必須等於這個數字, 這是D1最基本的檢核。)"""
    if direction not in ('x', 'y'):
        raise ValueError(f"direction必須是'x'或'y', 收到'{direction}'")
    total = 0.0
    for m in frame.members.values():
        section = frame.sections[m.section]
        if section.rho is None:
            continue
        L, _ = member_geometry(frame.nodes[m.node_i], frame.nodes[m.node_j])
        total += section.rho * section.A * L
    for nm in frame.node_masses:
        total += nm.mx if direction == 'x' else nm.my
    return float(total)
