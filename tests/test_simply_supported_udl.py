"""
驗證案例2: 簡支梁均佈載重 (分成2個元素,中點剛好落在節點上)
解析解:
  端點反力 R = wL/2
  中點撓度 v_mid = 5*w*L^4/(384*EI)   (向下, w向下取負)
  端彎矩應為0 (簡支,無彎矩支承)
"""
import numpy as np
from frame2d import Frame2D, solve

E = 200e6
I = 8e-5
A = 1e-2
L = 6.0     # 總跨度
w = -5.0    # kN/m, 向下 (local +y為正, 向下取負)

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L / 2, 0)
f.add_node(2, L, 0)
f.add_section('sec', E=E, I=I, A=A)
f.add_member(0, node_i=0, node_j=1, section='sec')
f.add_member(1, node_i=1, node_j=2, section='sec')

f.pin(0)
f.roller_y(2)
f.distributed_load(0, w=w)
f.distributed_load(1, w=w)

result = solve(f)

v_mid = result.displacements[f.dofs_of(1)[1]]
v_mid_exact = 5 * w * L**4 / (384 * E * I)   # w是負值,結果自然是負的(向下)

print(f"mid deflection: FEM={v_mid:.8e}  exact={v_mid_exact:.8e}  rel_err={(v_mid-v_mid_exact)/v_mid_exact:.3e}")
assert abs((v_mid - v_mid_exact) / v_mid_exact) < 1e-9, "中點撓度不符解析解"

R0 = result.reactions[f.dofs_of(0)[1]]
R2 = result.reactions[f.dofs_of(2)[1]]
R_exact = -w * L / 2   # 反力方向與載重相反
print(f"reactions: R0={R0:.6f}  R2={R2:.6f}  exact={R_exact:.6f}")
assert abs(R0 - R_exact) < 1e-8
assert abs(R2 - R_exact) < 1e-8

# 兩根桿件在跨中節點交接處, 端彎矩應該連續且等於簡支梁該處理論彎矩 wL^2/8 (取絕對值比較)
M_mid_from_member0 = result.member_results[0].end_forces_local[5]  # member0的j端(=node1)彎矩
M_mid_exact = abs(w) * L**2 / 8
print(f"mid moment: FEM={M_mid_from_member0:.6f}  exact(abs)={M_mid_exact:.6f}")
assert abs(abs(M_mid_from_member0) - M_mid_exact) / M_mid_exact < 1e-9

# 簡支端(pin/roller)彎矩應為0
M0 = result.member_results[0].end_forces_local[2]
M2 = result.member_results[1].end_forces_local[5]
assert abs(M0) < 1e-8, f"pin端彎矩應為0, 得到{M0}"
assert abs(M2) < 1e-8, f"roller端彎矩應為0, 得到{M2}"

print("PASS: 簡支梁均佈載重 完全吻合解析解 (含分佈載重固定端反力公式)")
