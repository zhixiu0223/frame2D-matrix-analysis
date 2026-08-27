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


def member_stiffness_local(E, I, A, L):
    """局部座標系下的 6x6 勁度矩陣 (軸向 + 彎曲耦合已分離,標準組合式)"""
    EA_L = E * A / L
    EI = E * I
    L2 = L * L
    L3 = L2 * L

    k = np.zeros((6, 6))

    # 軸向 (u1, u2) -> 索引 0, 3
    k[0, 0] = EA_L
    k[0, 3] = -EA_L
    k[3, 0] = -EA_L
    k[3, 3] = EA_L

    # 彎曲 (v1, theta1, v2, theta2) -> 索引 1,2,4,5
    kb = EI * np.array([
        [12 / L3,   6 / L2,  -12 / L3,   6 / L2],
        [6 / L2,    4 / L,   -6 / L2,    2 / L],
        [-12 / L3, -6 / L2,   12 / L3,  -6 / L2],
        [6 / L2,    2 / L,   -6 / L2,    4 / L],
    ])
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


def member_local_stiffness_dispatch(member_type, E, I, A, L):
    """依member_type選擇對應的局部勁度矩陣公式"""
    if member_type == 'truss':
        return member_stiffness_local_truss(E, A, L)
    return member_stiffness_local(E, I, A, L)


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


def member_stiffness_global(section, node_i, node_j, member_type='frame'):
    """組出全域座標系下的 6x6 勁度矩陣,回傳 (k_global, L, angle, T)"""
    L, angle = member_geometry(node_i, node_j)
    k_local = member_local_stiffness_dispatch(member_type, section.E, section.I, section.A, L)
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
