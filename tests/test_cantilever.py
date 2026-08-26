"""
驗證案例1: 懸臂梁自由端點載重
解析解 (教科書標準公式):
  tip deflection v = -P*L^3/(3EI)
  tip rotation    theta = -P*L^2/(2EI)
"""
import numpy as np
from frame2d import Frame2D, solve

E = 200e6   # kPa (=200 GPa, 單位用 kN/m^2)
I = 8e-5    # m^4
A = 1e-2    # m^2 (不影響此案例, 只有彎曲)
L = 4.0     # m
P = 10.0    # kN, 向下 (local +y 為正, 向下取負)

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('sec', E=E, I=I, A=A)
f.add_member(0, node_i=0, node_j=1, section='sec')
f.fix(0)
f.point_load(1, fy=-P)

result = solve(f)

v_tip = result.displacements[f.dofs_of(1)[1]]
theta_tip = result.displacements[f.dofs_of(1)[2]]

v_exact = -P * L**3 / (3 * E * I)
theta_exact = -P * L**2 / (2 * E * I)

print(f"tip v      : FEM={v_tip:.8e}  exact={v_exact:.8e}  rel_err={(v_tip-v_exact)/v_exact:.3e}")
print(f"tip theta  : FEM={theta_tip:.8e}  exact={theta_exact:.8e}  rel_err={(theta_tip-theta_exact)/theta_exact:.3e}")

assert abs((v_tip - v_exact) / v_exact) < 1e-10, "撓度不符解析解"
assert abs((theta_tip - theta_exact) / theta_exact) < 1e-10, "轉角不符解析解"

# 整體平衡檢查: 支承反力(此工具定義為"結構受到的外力", 即K@u-F) + 外加載重 應為0
Fy_reaction = result.reactions[f.dofs_of(0)[1]]
M_reaction = result.reactions[f.dofs_of(0)[2]]
assert abs(Fy_reaction - P) < 1e-8, f"垂直方向不平衡: reaction={Fy_reaction}, 應為{P}"
assert abs(M_reaction - P * L) < 1e-8, f"彎矩不平衡: reaction={M_reaction}, 應為{P*L}"

# 迴歸測試: postprocess.member_internal_forces() 沿桿長的M(x)公式,
# 用這個案例(M1=40非零)專門鎖住之前抓到的正負號bug——
# 簡支梁(M1=M2≈0)測不出這個bug, 因為正負號怎麼取結果都是0
from frame2d.postprocess import member_internal_forces
x, N, V, M = member_internal_forces(f, result, 0, n=5)
ef = result.member_results[0].end_forces_local
M1_local = ef[2]
assert abs(M[0] - (-M1_local)) < 1e-8, f"M(0)應為-M1={-M1_local}, 得到{M[0]}"
assert abs(M[-1] - ef[5]) < 1e-6, f"M(L)應為M2={ef[5]}, 得到{M[-1]}"
print(f"PASS: member_internal_forces() M(x)公式正確 (M(0)=-M1, M(L)=M2, 非零M1情況下驗證)")

print("PASS: 懸臂梁點載重 完全吻合解析解")
