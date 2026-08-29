"""
驗證案例9: 桿件內部(不在節點上)的集中力/集中力矩 (Load System v2 Phase 1)

驗證策略: 優先用「簡支梁的支承反力」當基準——簡支梁是靜定結構, 反力可以
純用靜力學(力平衡+力矩平衡)算出來, 完全不依賴 fixed_end_forces_point_load
公式本身對不對, 是最不會被自己的bug混淆的驗證方式。
再用懸臂梁的解析解撓度公式做第二層交叉確認。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 1e-2


def check(label, fem, exact, tol=1e-9):
    rel = abs(fem - exact) / max(abs(exact), 1e-12)
    print(f"  {label}: FEM={fem:.8f}  解析解={exact:.8f}  rel_err={rel:.2e}")
    assert rel < tol, f"{label} 不吻合"


# ---- 案例A: 簡支梁, 內部橫向點載重(不對稱位置a≠b), 對照純靜力學決定性反力 ----
print("=== 案例A: 簡支梁內部點載重, 對照純靜力學反力(不依賴固定端公式本身) ===")
L, P, a = 10.0, 30.0, 3.5   # a=3.5 (刻意不對稱, b=6.5, 避免對稱案例矇混過去)
b = L - a
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.member_point_load(0, a=a, fy=-P)   # 向下
r = solve(f)

R1_exact = P * b / L    # 純靜力學: 對節點1取矩 -> R0*L = P*b
R2_exact = P * a / L
R1_fem = r.reactions[f.dofs_of(0)[1]]
R2_fem = r.reactions[f.dofs_of(1)[1]]
check("R0(pin)", R1_fem, R1_exact)
check("R1(roller)", R2_fem, R2_exact)
print("PASS: 簡支梁反力精確吻合純靜力學決定性解\n")


# ---- 案例B: 懸臂梁, 內部橫向點載重, 對照解析解撓度公式 ----
print("=== 案例B: 懸臂梁內部點載重, 對照解析解撓度公式 ===")
L, P, a = 8.0, 25.0, 3.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.fix(0)
f.member_point_load(0, a=a, fy=-P)
r = solve(f)

# 標準公式: 懸臂梁(固定於0), 點載重P(向下)在x=a, 自由端(x=L)撓度:
#   v(L) = -P*a^2*(3L-a) / (6EI)
v_tip_exact = -P * a**2 * (3 * L - a) / (6 * E * I)
v_tip_fem = r.displacements[f.dofs_of(1)[1]]
check("tip撓度", v_tip_fem, v_tip_exact)

# 固定端反力: Fy=P(平衡), M=P*a(力臂)
Fy_exact = P
M_exact = P * a
check("固定端Fy反力", r.reactions[f.dofs_of(0)[1]], Fy_exact)
check("固定端M反力", r.reactions[f.dofs_of(0)[2]], M_exact)
print("PASS: 懸臂梁內部點載重撓度+反力完全吻合解析解\n")


# ---- 案例C: 簡支梁, 內部集中力矩, 對照純靜力學決定性反力 ----
print("=== 案例C: 簡支梁內部集中力矩, 對照純靜力學反力 ===")
L, M0, a = 12.0, 40.0, 5.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.member_point_load(0, a=a, m=M0)
r = solve(f)

# 純靜力學: 對節點0取矩 -> R1*L + M0 = 0 (CCW力矩M0在任一點的力矩貢獻都是+M0,
# 不因位置改變, 力矩平衡: R1*L + M0 = 0 -> R1 = -M0/L; 垂直平衡: R0 = -R1 = +M0/L)
R2_exact = -M0 / L
R1_exact = M0 / L
check("R0(pin)", r.reactions[f.dofs_of(0)[1]], R1_exact)
check("R1(roller)", r.reactions[f.dofs_of(1)[1]], R2_exact)
print("PASS: 簡支梁內部集中力矩反力精確吻合純靜力學決定性解\n")


# ---- 案例D: 桿件內部軸向點載重, 對照彈簧串聯解析解 ----
print("=== 案例D: 桿件內部軸向點載重, 對照彈簧串聯解析解 ===")
L, P, a = 6.0, 50.0, 2.0
b = L - a
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.fix(0)
f.fix(1)
f.member_point_load(0, a=a, fx=P)
r = solve(f)

F1_exact = P * b / L
F2_exact = P * a / L
check("F0(固定端軸力反力)", r.reactions[f.dofs_of(0)[0]], -F1_exact)
check("F1(固定端軸力反力)", r.reactions[f.dofs_of(1)[0]], -F2_exact)
print("PASS: 軸向點載重反力吻合彈簧串聯解析解\n")

print("PASS: Phase 1(桿件中間集中力/力矩) 全數驗證通過")


# ---- 案例E: 內部力/力矩診斷圖 V(x)/M(x) 是否正確反映跳躍不連續 ----
print()
print("=== 案例E: postprocess的V(x)/M(x)診斷圖是否正確反映點載重/點力矩的不連續 ===")
from frame2d.postprocess import member_internal_forces

L, P, a = 10.0, 30.0, 3.5
b = L - a
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.member_point_load(0, a=a, fy=-P)
r = solve(f)
x, N, V, M = member_internal_forces(f, r, 0, n=41)
R0 = P * b / L
for xi, Vi, Mi in zip(x, V, M):
    if xi < a - 1e-9:
        V_exact, M_exact = R0, R0 * xi
    else:
        V_exact, M_exact = R0 - P, R0 * xi - P * (xi - a)
    assert abs(Vi - V_exact) < 1e-8, f"V(x={xi}) 不吻合"
    assert abs(Mi - M_exact) < 1e-8, f"M(x={xi}) 不吻合"
print("PASS: 點載重造成的V(x)剪力跳躍、M(x)斜率變化, 41個取樣點全數吻合純靜力學公式")

L, M0, a = 12.0, 40.0, 5.0
b = L - a
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's')
f2.pin(0)
f2.roller_y(1)
f2.member_point_load(0, a=a, m=M0)
r2 = solve(f2)
x2, N2, V2, M2 = member_internal_forces(f2, r2, 0, n=41)
R0_2 = M0 / L
R1_2 = -M0 / L
for xi, Mi in zip(x2, M2):
    if xi < a - 1e-9:
        M_exact = R0_2 * xi
    else:
        M_exact = R1_2 * (L - xi)
    assert abs(Mi - M_exact) < 1e-8, f"M(x={xi}) 不吻合(點力矩案例)"
print("PASS: 點力矩造成的M(x)彎矩跳躍, 41個取樣點全數吻合純靜力學公式")

