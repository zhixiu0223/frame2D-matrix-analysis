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

print("PASS: 懸臂梁點載重 完全吻合解析解")
