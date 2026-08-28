"""
驗證案例14: Phase 3(支承強制位移)對照 SW FEA 第三方工具

使用者上傳兩個案例(phase-03-rot, phase-03-sinking), 都是同一個固接-固接
門型鋼架: 節點0=(0,0)、節點1=(6,0)為基礎, 節點2=(0,4)、節點3=(6,4)為柱頂,
F0=柱(0-2)、F1=梁(2-3)、F2=柱(3-1)。斷面 E=200GPa, I=8.0E7mm^4(=8e-5m^4),
A=0.01m^2, 跟本專案的kN,m慣例直接對應, 不需要單位換算。

案例A(phase-03-rot): 節點0被強制轉0.3 rad(沒有外加載重, 純粹是支承
  強制轉角造成的內力)。
案例B(phase-03-sinking): 節點0下沉10mm(轉角仍固定=0)。

驗證用SW FEA報告的完整逐點資料(每根桿件11個位置點的SF/BM及局部座標
dX/dY), 這是目前為止最完整的一次交叉驗證(2案例 x 3桿件 x 11點 x 多個
物理量)。

**重要note**: SW FEA每根桿件自己報告表格裡的dX/dY欄位, 是桿件"局部座標"
(dX=沿桿軸方向, dY=垂直桿軸方向), 不是全域座標——這點在比對垂直桿件
(例如F0是柱子)時特別容易搞混(垂直桿件的局部座標系跟全域XY差了90度旋轉),
已經用這點自己驗證過一次(見下方比對邏輯)。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces
from frame2d.elements import member_geometry, transformation_matrix

E, I, A = 200e6, 8e-5, 0.01


def _build(rot0, uy0):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')   # F0
    f.add_member(1, 2, 3, 's')   # F1
    f.add_member(2, 3, 1, 's')   # F2
    f.support(0, ux=0.0, uy=uy0, rot=rot0)
    f.support(1, ux=0.0, uy=0.0, rot=0.0)
    return f


def _local_deformed(f, r, mid, n=11):
    """回傳桿件局部座標系的(x, dX軸向, dY橫向)位移分布, 對照SW FEA每根桿件
    報告表格用的同一套局部座標慣例。"""
    m = f.members[mid]
    ni, nj = f.nodes[m.node_i], f.nodes[m.node_j]
    L, angle = member_geometry(ni, nj)
    T = transformation_matrix(angle)
    dofs = f.dofs_of(m.node_i) + f.dofs_of(m.node_j)
    u1, v1, th1, u2, v2, th2 = T @ r.displacements[np.array(dofs)]
    EI = E * I
    x_full, N_f, V_f, M_f = member_internal_forces(f, r, mid, n=101)
    kappa = M_f / EI
    slope = np.concatenate([[0.0], np.cumsum((kappa[1:] + kappa[:-1]) / 2 * np.diff(x_full))])
    defl = np.concatenate([[0.0], np.cumsum((slope[1:] + slope[:-1]) / 2 * np.diff(x_full))])
    C2, C1 = v1, (v2 - v1 - defl[-1]) / L
    x = np.linspace(0, L, n)
    dY = np.interp(x, x_full, defl + C1 * x_full + C2)
    dX = u1 + (u2 - u1) * (x / L)
    return x, dX, dY


def _check_case(f, r, sw_reactions, sw_members, label):
    print(f"=== {label} ===")
    for nid, (Rx, Ry, M) in sw_reactions.items():
        ux, uy, rot = f.dofs_of(nid)
        assert abs(r.reactions[ux] - Rx) < 0.02, f"node{nid} Rx不吻合"
        assert abs(r.reactions[uy] - Ry) < 0.02, f"node{nid} Ry不吻合"
        assert abs(r.reactions[rot] - M) < 0.02, f"node{nid} M不吻合"
    print("  反力: PASS")

    max_err_sf, max_err_bm, max_err_dx, max_err_dy = 0, 0, 0, 0
    for mid, points in sw_members.items():
        x, N, V, M = member_internal_forces(f, r, mid, n=11)
        xl, dX, dY = _local_deformed(f, r, mid, n=11)
        for i, (sf, bm, dxsw, dysw) in enumerate(points):
            max_err_sf = max(max_err_sf, abs(V[i] - sf))
            max_err_bm = max(max_err_bm, abs(M[i] - bm))
            max_err_dx = max(max_err_dx, abs(dX[i] * 1000 - dxsw))
            max_err_dy = max(max_err_dy, abs(dY[i] * 1000 - dysw))
    print(f"  逐點最大誤差: SF={max_err_sf:.4f}  BM={max_err_bm:.4f}  dX={max_err_dx:.4f}mm  dY={max_err_dy:.4f}mm")
    assert max_err_sf < 0.01 and max_err_bm < 0.01 and max_err_dx < 0.01 and max_err_dy < 0.01
    print("  逐點SF/BM/dX/dY: PASS\n")


# ---- 案例A: 節點0強制轉角0.3 rad (無外加載重) ----
fA = _build(rot0=0.3, uy0=0.0)
rA = solve(fA)
sw_reactions_A = {0: (-561.395, 159.886, 2426.501), 1: (561.395, -159.886, -1467.183)}
sw_members_A = {
    0: [(561.395, -2426.501, 0.000, 0.000), (561.395, -2201.943, -0.032, 108.242),
        (561.395, -1977.385, -0.064, 194.464), (561.395, -1752.827, -0.096, 260.913),
        (561.395, -1528.269, -0.128, 309.833), (561.395, -1303.712, -0.160, 343.470),
        (561.395, -1079.154, -0.192, 364.071), (561.395, -854.596, -0.224, 373.880),
        (561.395, -630.038, -0.256, 375.142), (561.395, -405.480, -0.288, 370.105),
        (561.395, -180.922, -0.320, 361.013)],
    1: [(159.886, -180.922, -361.013, -0.320), (159.886, -84.990, -360.844, -17.552),
        (159.886, 10.942, -360.676, -36.697), (159.886, 106.873, -360.507, -55.595),
        (159.886, 202.805, -360.339, -72.089), (159.886, 298.737, -360.171, -84.020),
        (159.886, 394.669, -360.002, -89.229), (159.886, 490.600, -359.834, -85.558),
        (159.886, 586.532, -359.665, -70.848), (159.886, 682.464, -359.497, -42.942),
        (159.886, 778.396, -359.328, 0.320)],
    2: [(-561.395, 778.396, -0.320, -359.328), (-561.395, 553.838, -0.288, -321.371),
        (-561.395, 329.280, -0.256, -277.876), (-561.395, 104.722, -0.224, -231.088),
        (-561.395, -119.836, -0.192, -183.252), (-561.395, -344.394, -0.160, -136.615),
        (-561.395, -568.952, -0.128, -93.422), (-561.395, -793.510, -0.096, -55.918),
        (-561.395, -1018.067, -0.064, -26.350), (-561.395, -1242.625, -0.032, -6.962),
        (-561.395, -1467.183, 0.000, 0.000)],
}
_check_case(fA, rA, sw_reactions_A, sw_members_A, "案例A: 支承強制轉角0.3rad")


# ---- 案例B: 節點0下沉10mm ----
fB = _build(rot0=0.0, uy0=-0.01)
rB = solve(fB)
sw_reactions_B = {0: (-0.000, -1.777, -5.330), 1: (-0.000, 1.777, -5.330)}
sw_members_B = {
    0: [(0.000, 5.330, -10.000, 0.000), (0.000, 5.330, -10.000, 0.027), (0.000, 5.330, -9.999, 0.107),
        (0.000, 5.330, -9.999, 0.240), (0.000, 5.330, -9.999, 0.426), (0.000, 5.330, -9.998, 0.666),
        (0.000, 5.330, -9.998, 0.959), (0.000, 5.330, -9.998, 1.306), (0.000, 5.330, -9.997, 1.705),
        (0.000, 5.330, -9.997, 2.158), (0.000, 5.330, -9.996, 2.665)],
    1: [(-1.777, 5.330, -2.665, -9.996), (-1.777, 4.264, -2.665, -9.141), (-1.777, 3.198, -2.665, -8.190),
        (-1.777, 2.132, -2.665, -7.166), (-1.777, 1.066, -2.665, -6.095), (-1.777, 0.000, -2.665, -5.000),
        (-1.777, -1.066, -2.665, -3.905), (-1.777, -2.132, -2.665, -2.834), (-1.777, -3.198, -2.665, -1.810),
        (-1.777, -4.264, -2.665, -0.859), (-1.777, -5.330, -2.665, -0.004)],
    2: [(0.000, -5.330, 0.004, -2.665), (0.000, -5.330, 0.003, -2.158), (0.000, -5.330, 0.003, -1.705),
        (0.000, -5.330, 0.002, -1.306), (0.000, -5.330, 0.002, -0.959), (0.000, -5.330, 0.002, -0.666),
        (0.000, -5.330, 0.001, -0.426), (0.000, -5.330, 0.001, -0.240), (0.000, -5.330, 0.001, -0.107),
        (0.000, -5.330, 0.000, -0.027), (0.000, -5.330, 0.000, 0.000)],
}
_check_case(fB, rB, sw_reactions_B, sw_members_B, "案例B: 支承沉陷10mm")

print("PASS: Phase 3(支承強制位移) 對照SW FEA, 2案例x3桿件x11點, 全數驗證完成")
