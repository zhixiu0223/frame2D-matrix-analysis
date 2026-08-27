"""
後處理模組: 結構圖 / 受力圖 / 軸力圖 / 剪力圖 / 彎矩圖 / 變形圖

設計原則(呼應規劃階段的討論): 不要讓繪圖程式自己重新計算力學,
所有數值都是從 SolveResult 這份「唯一的正確答案」推導出來的不同視角。

桿端內力 V(x)/M(x) 公式推導與驗證:
  V(x) = Fy1_local + W(x)              W(x) = 0到x的累積分佈載重
  M(x) = -M1_local + Fy1_local*x + ∫[0,x] W(s) ds

  注意 M(x) 前面是 -M1_local (不是+M1_local): 這是 FEM 節點力慣例
  (M1,M2都是CCW正, 兩端同號) 跟古典彎矩圖慣例(全長連續, 上下拉力側正負一致)
  之間的標準轉換關係 M(0)=-M1, M(L)=+M2。之前只用簡支梁(M1=M2≈0)測過
  這條公式, 兩端剛好都是0, 正負號錯了也測不出來; 後來拿懸臂梁(M1=40非零)
  反推才抓到 M(L)應該等於M2, 而不是+M1+Fy1*L那樣算出來的80。
  已重新驗證: 懸臂梁 M(L)=0=M2 (原本錯誤版本算出80), 簡支梁跨中M=22.5=wL²/8 依然吻合。

變形內插用 Hermite cubic shape function (homogeneous解), 節點值精確,
但若桿件跨中有分佈載重且只用單一元素代表整根桿件, 內插的跨中撓度會略微
低估實際下垂量(因為真實解在均佈載重下是四次多項式, 不是三次)——
若要更精確的變形視覺化, 建議把該桿件切成多個 Member 元素。
"""
import numpy as np


def _w_at(w_start, w_end, x, L):
    return w_start + (w_end - w_start) * x / L


def _cumulative_W(w_start, w_end, x, L):
    """W(x) = integral of w(s) ds from 0 to x"""
    dw = w_end - w_start
    return w_start * x + dw * x**2 / (2 * L)


def _cumulative_M_from_W(w_start, w_end, x, L):
    """integral of W(s) ds from 0 to x (雙重積分, 用於M(x))"""
    dw = w_end - w_start
    return w_start * x**2 / 2 + dw * x**3 / (6 * L)


def member_internal_forces(frame, result, member_id, n=21):
    """回傳局部座標系下沿桿長分佈的 (x, N, V, M) 陣列"""
    mr = result.member_results[member_id]
    L = mr.L
    Fx1, Fy1, M1 = mr.end_forces_local[0], mr.end_forces_local[1], mr.end_forces_local[2]

    w_start = w_end = 0.0
    for dl in frame.distributed_loads:
        if dl.member == member_id:
            w_start += dl.w_start
            w_end += dl.w_end

    x = np.linspace(0, L, n)
    N = np.full(n, Fx1)
    W = _cumulative_W(w_start, w_end, x, L)
    V = Fy1 + W
    M = -M1 + Fy1 * x + _cumulative_M_from_W(w_start, w_end, x, L)
    return x, N, V, M


def _hermite_shape(x, L):
    """回傳 Hermite cubic shape function 在 x 處的值 [N1,N2,N3,N4] (對應 v1,theta1,v2,theta2)"""
    xi = x / L
    N1 = 1 - 3 * xi**2 + 2 * xi**3
    N2 = L * (xi - 2 * xi**2 + xi**3)
    N3 = 3 * xi**2 - 2 * xi**3
    N4 = L * (-xi**2 + xi**3)
    return N1, N2, N3, N4


def member_deformed_shape(frame, result, member_id, scale=1.0, n=21):
    """回傳桿件變形後的全域座標 (X, Y) 陣列 (含放大係數scale)"""
    from .elements import member_geometry, transformation_matrix
    m = frame.members[member_id]
    node_i = frame.nodes[m.node_i]
    node_j = frame.nodes[m.node_j]
    L, angle = member_geometry(node_i, node_j)
    T = transformation_matrix(angle)

    dofs = frame.dofs_of(m.node_i) + frame.dofs_of(m.node_j)
    u_global_member = result.displacements[np.array(dofs)]
    u_local = T @ u_global_member
    u1, v1, th1, u2, v2, th2 = u_local

    x = np.linspace(0, L, n)
    X_local = np.zeros(n)
    Y_local = np.zeros(n)
    for k, xx in enumerate(x):
        N1, N2, N3, N4 = _hermite_shape(xx, L)
        v = N1 * v1 + N2 * th1 + N3 * v2 + N4 * th2
        u_ax = u1 + (u2 - u1) * (xx / L)   # 軸向位移線性內插
        X_local[k] = xx + scale * u_ax
        Y_local[k] = scale * v

    # 轉回全域座標 (局部->全域是 T的轉置, 再加上桿件原點座標)
    c, s = np.cos(angle), np.sin(angle)
    X_global = node_i.x + X_local * c - Y_local * s
    Y_global = node_i.y + X_local * s + Y_local * c
    return X_global, Y_global
