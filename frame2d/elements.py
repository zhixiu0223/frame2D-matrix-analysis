"""
2D Euler-Bernoulli 樑柱元素公式

局部自由度順序: [u1, v1, theta1, u2, v2, theta2]
  u = 沿桿軸方向 (軸向)
  v = 垂直桿軸方向 (撓曲)
  theta = 轉角

這是整個專案唯一有「公式對不對」風險的地方,獨立寫成小函式方便單元測試。
"""
import numpy as np


def member_geometry(node_i, node_j):
    """回傳 (L, angle) : 長度與桿件相對於全域x軸的夾角(弧度)"""
    dx = node_j.x - node_i.x
    dy = node_j.y - node_i.y
    L = float(np.hypot(dx, dy))
    angle = float(np.arctan2(dy, dx))
    return L, angle


def member_stiffness_local(E, I, A, L, release_i=False, release_j=False):
    """局部座標系下的 6x6 勁度矩陣 (軸向 + 彎曲耦合已分離,標準組合式)。

    release_i/release_j: 該端是否有內部鉸接(彎矩釋放, M=0)。用靜力凝縮
    (static condensation)處理, 不改變全域DOF系統(dofs_of()完全不用動,
    這是這個功能不需要先做DOFManager升級就能實作的原因)。公式已用sympy
    驗證跟標準教科書「一端鉸接梁」的3EI/L係數一致(標準4EI/L,6EI/L²改成
    3EI/L,3EI/L², 12EI/L³不變)。
    只支援單端釋放或不釋放; 兩端同時釋放(等同truss, 但這裡沒有處理桿件
    內部有載重時的固定端反力修正公式)不支援, 兩端都釋放時只用軸向+已經
    是0的彎曲勁度, 不保證載重情況正確, 請改用add_truss()。
    """
    EA_L = E * A / L
    EI = E * I
    L2 = L * L
    L3 = L2 * L

    k = np.zeros((6, 6))

    # 軸向 (u1, u2) -> 索引 0, 3 (跟release無關)
    k[0, 0] = EA_L
    k[0, 3] = -EA_L
    k[3, 0] = -EA_L
    k[3, 3] = EA_L

    # 彎曲 (v1, theta1, v2, theta2) -> 索引 1,2,4,5
    if not release_i and not release_j:
        kb = EI * np.array([
            [12 / L3,   6 / L2,  -12 / L3,   6 / L2],
            [6 / L2,    4 / L,   -6 / L2,    2 / L],
            [-12 / L3, -6 / L2,   12 / L3,  -6 / L2],
            [6 / L2,    2 / L,   -6 / L2,    4 / L],
        ])
    elif release_j and not release_i:
        # J端(theta2)鉸接: 標準4EI/L,6EI/L²改成3EI/L,3EI/L², theta2那一列/行全為0
        kb = EI * np.array([
            [3 / L3,   3 / L2,  -3 / L3,   0.0],
            [3 / L2,   3 / L,   -3 / L2,   0.0],
            [-3 / L3, -3 / L2,   3 / L3,   0.0],
            [0.0,      0.0,      0.0,      0.0],
        ])
    elif release_i and not release_j:
        # I端(theta1)鉸接: 鏡像版本, theta1那一列/行全為0
        kb = EI * np.array([
            [3 / L3,   0.0,  -3 / L3,   3 / L2],
            [0.0,      0.0,   0.0,      0.0],
            [-3 / L3,  0.0,   3 / L3,  -3 / L2],
            [3 / L2,   0.0,  -3 / L2,   3 / L],
        ])
    else:
        # 兩端都釋放: 彎曲勁度全為0(等同truss, 但不處理內部載重的固定端反力)
        kb = np.zeros((4, 4))

    bend_idx = [1, 2, 4, 5]
    for a, ia in enumerate(bend_idx):
        for b, ib in enumerate(bend_idx):
            k[ia, ib] = kb[a, b]

    return k


def member_stiffness_local_truss(E, A, L):
    """局部座標系下的 6x6 勁度矩陣, 純桁架(軸力)元素:
    只有軸向(u1,u2, 索引0,3)有勁度, 彎曲/剪力相關項全為0
    (跟frame共用同一個6x6格式+同一個transformation_matrix, 組裝程式碼不用改,
    只是彎曲block是零矩陣 -- 物理上代表這根桿件兩端都是鉸接, 只能傳軸力)
    """
    EA_L = E * A / L
    k = np.zeros((6, 6))
    k[0, 0] = EA_L
    k[0, 3] = -EA_L
    k[3, 0] = -EA_L
    k[3, 3] = EA_L
    return k


def member_local_stiffness_dispatch(member_type, E, I, A, L, release_i=False, release_j=False):
    """依member_type選擇對應的局部勁度矩陣公式。
    'truss'跟'cable'共用同一條軸力公式(兩端鉸接、只傳軸力), 差別在
    solve.py會不會把受壓的cable桿件當成鬆弛移除, truss則不管拉壓都保留。
    'frame'則依release_i/release_j決定要不要做端點鉸接的靜力凝縮。"""
    if member_type in ('truss', 'cable'):
        return member_stiffness_local_truss(E, A, L)
    return member_stiffness_local(E, I, A, L, release_i, release_j)


def transformation_matrix(angle):
    """全域 -> 局部 的 6x6 旋轉矩陣 T,使 k_global = T^T @ k_local @ T"""
    c, s = np.cos(angle), np.sin(angle)
    r = np.array([
        [c,  s, 0],
        [-s, c, 0],
        [0,  0, 1],
    ])
    T = np.zeros((6, 6))
    T[0:3, 0:3] = r
    T[3:6, 3:6] = r
    return T


def member_stiffness_global(section, node_i, node_j, member_type='frame', release_i=False, release_j=False):
    """組出全域座標系下的 6x6 勁度矩陣,回傳 (k_global, L, angle, T)"""
    L, angle = member_geometry(node_i, node_j)
    k_local = member_local_stiffness_dispatch(member_type, section.E, section.I, section.A, L, release_i, release_j)
    T = transformation_matrix(angle)
    k_global = T.T @ k_local @ T
    return k_global, L, angle, T


def fixed_end_forces_udl(w_start, w_end, L):
    """局部座標系下,垂直均佈/線性變化載重的固定端反力(彎矩+剪力)。
    符號慣例: w 沿局部 +y 方向為正。
    回傳長度6向量 [Fx1, Fy1, M1, Fx2, Fy2, M2] (這是"固定端反力",
    要當作等效節點載重時,組裝進F時要取負號)。

    使用一般線性變化分布載重公式 (梯形載重), w_start=w_end 時退化為均佈載重標準式:
      均佈 w: V = wL/2 , M1 = wL^2/12 , M2 = -wL^2/12
    """
    if abs(w_start - w_end) < 1e-12:
        w = w_start
        V = w * L / 2.0
        M1 = w * L * L / 12.0
        M2 = -w * L * L / 12.0
        return np.array([0.0, V, M1, 0.0, V, M2])

    # 線性變化載重: 拆成均佈分量 w_start + 三角形分量 (w_end-w_start)
    w0 = w_start
    dw = w_end - w_start

    # 均佈部分 w0
    V1_u = w0 * L / 2.0
    M1_u = w0 * L * L / 12.0
    V2_u = w0 * L / 2.0
    M2_u = -w0 * L * L / 12.0

    # 三角形部分 (0 -> dw), 標準固定端反力公式
    V1_t = 3.0 * dw * L / 20.0
    M1_t = dw * L * L / 30.0
    V2_t = 7.0 * dw * L / 20.0
    M2_t = -dw * L * L / 20.0

    V1 = V1_u + V1_t
    M1 = M1_u + M1_t
    V2 = V2_u + V2_t
    M2 = M2_u + M2_t

    return np.array([0.0, V1, M1, 0.0, V2, M2])


def fixed_end_forces_point_load(P, a, L):
    """局部座標系下, 桿件內部任意位置(距node_i為a, 0<=a<=L)的橫向點載重P
    (沿局部+y方向為正)之固定端反力。回傳格式跟 fixed_end_forces_udl 一致
    (直接當work-equivalent等效節點載重使用, 不取負號)。

    公式來源: 用EI*v''=M(x)分段雙重積分, 配合固定-固定邊界條件, 從梁的
    微分方程直接推導(非憑記憶抄書, 避免重蹈M(x)公式曾經正負號抄反的錯誤),
    並用簡支梁決定性反力(R1=Pb/L, R2=Pa/L, 純靜力學、不受此公式正確性影響)
    交叉驗證過, 見 tests/test_member_point_load.py。
    """
    b = L - a
    V1 = P * b**2 * (3 * a + b) / L**3
    M1 = P * a * b**2 / L**2
    V2 = P * a**2 * (3 * b + a) / L**3
    M2 = -P * a**2 * b / L**2
    return np.array([0.0, V1, M1, 0.0, V2, M2])


def fixed_end_forces_point_moment(M0, a, L):
    """局部座標系下, 桿件內部任意位置(距node_i為a)的集中力矩M0
    (逆時針為正, 跟theta同號約定)之固定端反力。回傳格式同上。

    **這裡曾經有一個真的正負號bug**: 原本的實作(V1,M1,V2,M2全部反號)是把
    fixed_end_forces_point_load()「sympy原始推導結果要整體取負號才對」這個
    規律, 直接套用到力矩案例上, 但兩種案例的sympy推導約定不是同一回事,
    不能照搬。用「懸臂梁自由端直接加外力矩」(固定端反力矩應該=-M0, 梁彎矩
    應為常數+M0, 這兩個都是最基本、不依賴這條公式本身的物理事實)、以及
    「節點分割法」(在a處插入真實節點, 改用經過獨立驗證的point_load的m參數,
    比較兩者是否給出相同答案)這兩種方法互相佐證, 抓到並修正: 這裡不需要
    對sympy原始推導結果取負號, 直接用sympy的原始輸出即可(跟point_load的
    情況相反)。見 tests/test_member_moment_sign_consistency.py。
    """
    b = L - a
    V1 = -6.0 * M0 * a * b / L**3
    M1 = -M0 * b * (2 * a - b) / L**2
    V2 = 6.0 * M0 * a * b / L**3
    M2 = -M0 * a * (2 * b - a) / L**2
    return np.array([0.0, V1, M1, 0.0, V2, M2])


def fixed_end_forces_axial_point_load(P, a, L):
    """局部座標系下, 桿件內部任意位置(距node_i為a)的軸向集中力P
    (沿局部+x方向為正)之固定端反力。

    推導: 固定-固定桿件受軸向點載重, 等同兩根軸向彈簧(勁度EA/a, EA/b)
    在中間節點串聯承受外力P, 節點位移 u=P*a*b/(EA*L), 兩端固定端反力
    F1=EA/a*u=P*b/L, F2=EA/b*u=P*a/L (簡單物理量, 不需要sympy推導)。
    """
    b = L - a
    F1 = P * b / L
    F2 = P * a / L
    return np.array([F1, 0.0, 0.0, F2, 0.0, 0.0])


_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(6)


def fixed_end_forces_partial_udl(w_start, w_end, c, d, L):
    """局部座標系下, 桿件內部局部段 [c,d] (0<=c<=d<=L, 可以不是整根桿件)
    的均佈/線性變化載重之固定端反力。

    推導方式: 不手動謄寫龐大的封閉式展開式(降低抄寫出錯風險), 改用高斯-勒讓德
    數值積分, 對已經驗證過的 fixed_end_forces_point_load() 在 [c,d] 區間
    積分(等於把分佈載重拆成無限多個點載重疊加)。被積函數是 s 的4次多項式
    (point load公式本身是a的3次式, 乘上w(s)這個s的1次式), 用6點高斯積分
    (精確積分到11次多項式)對這個被積函數是"數值精確解", 不是近似。
    c=0, d=L 時應該退化成 fixed_end_forces_udl() 的結果(見
    tests/test_partial_udl.py的交叉驗證)。
    """
    if d <= c:
        return np.zeros(6)
    jac = 0.5 * (d - c)
    s_vals = jac * _GAUSS_NODES + 0.5 * (d + c)
    f_FE = np.zeros(6)
    for si, wi in zip(s_vals, _GAUSS_WEIGHTS):
        w_si = w_start + (w_end - w_start) * (si - c) / (d - c)
        f_FE += wi * jac * fixed_end_forces_point_load(w_si, si, L)
    return f_FE


def fixed_end_forces_udl_release_j(w_start, w_end, L):
    """局部座標系下, J端有鉸接釋放(M2=0)的桿件, 受整根桿件均佈載重的
    固定端反力(標準"一端固接一端鉸接梁", propped cantilever, 公式)。
    用EI*v''=M(x)配合天然邊界條件M(L)=0(取代一般固接情況的v'(L)=0)
    推導, 均佈載重的標準結果是wL²/8(比固接情況的wL²/12大, 符合物理直覺:
    少了一端的轉角拘束, 固接端要扛的彎矩更大)。只支援均佈(w_start=w_end),
    線性變化的鉸接固定端反力公式還沒推導。
    """
    if abs(w_start - w_end) > 1e-9:
        raise NotImplementedError("目前只支援均佈(w_start=w_end)的鉸接固定端反力公式")
    w = w_start
    V1 = 5 * w * L / 8
    M1 = w * L**2 / 8
    V2 = w * L - V1
    return np.array([0.0, V1, M1, 0.0, V2, 0.0])


def fixed_end_forces_udl_release_i(w_start, w_end, L):
    """局部座標系下, I端有鉸接釋放(M1=0)的桿件, 受整根桿件均佈載重的
    固定端反力(鏡像版本)。"""
    if abs(w_start - w_end) > 1e-9:
        raise NotImplementedError("目前只支援均佈(w_start=w_end)的鉸接固定端反力公式")
    w = w_start
    V2 = 5 * w * L / 8
    M2 = -w * L**2 / 8
    V1 = w * L - V2
    return np.array([0.0, V1, 0.0, 0.0, V2, M2])
