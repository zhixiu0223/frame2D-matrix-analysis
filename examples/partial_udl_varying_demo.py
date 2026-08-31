"""
Part 1: 簡單模型驗證 -- 簡支梁中間段受"不均勻(線性變化/梯形)"分布力。

這個功能(distributed_load的x_start/x_end局部段 + w_start!=w_end線性變化)
其實在Load System v2 Phase 2就已經支援並驗證過了(見tests/test_partial_udl.py
案例C: 簡支梁局部段梯形載重, 反力對照純靜力學決定性解, 精確吻合)。
這支腳本只是把同一個驗證案例畫出來, 當作Part 2(不均勻雪載重屋頂模型)的
暖身/對照組。

跑法: PYTHONPATH=. python examples/partial_udl_varying_demo.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

E, I, A = 200e6, 8e-5, 1e-2
L = 10.0

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
# 中間段[2,6]受不均勻(梯形)分布力: 從-10kN/m線性變化到-30kN/m
f.distributed_load(0, w=-10.0, w_end=-30.0, x_start=2.0, x_end=6.0)

r = solve(f)

# 對照純靜力學決定性反力(見tests/test_partial_udl.py案例C的推導)
c, d = 2.0, 6.0
w1v, w2v = -10.0, -30.0
span = d - c
W_uniform = w1v * span
centroid_uniform = (c + d) / 2
dw = w2v - w1v
W_tri = 0.5 * dw * span
centroid_tri = c + 2.0 / 3.0 * span
W_total = W_uniform + W_tri
moment_about_0 = W_uniform * centroid_uniform + W_tri * centroid_tri
centroid_total = moment_about_0 / W_total
R0_exact = -W_total * (L - centroid_total) / L
R1_exact = -W_total * centroid_total / L

print("節點反力 (簡支梁, 中間段不均勻分布力 -10→-30 kN/m, x=[2,6]):")
ux0, uy0, rot0 = f.dofs_of(0)
ux1, uy1, rot1 = f.dofs_of(1)
print(f"  node0(pin):    Ry={r.reactions[uy0]:.4f}  (純靜力學解析解={R0_exact:.4f})")
print(f"  node1(roller): Ry={r.reactions[uy1]:.4f}  (純靜力學解析解={R1_exact:.4f})")

fig = plot_all(f, r, figsize=(15, 9))
fig.suptitle("Part 1: Simply-Supported Beam, Non-Uniform (Trapezoidal) "
             "Midspan Load -10 to -30 kN/m", fontsize=11)
fig.tight_layout()
fig.savefig(f'{OUT}/partial_udl_varying.png', dpi=130)
plt.close(fig)
print(f"\n六合一圖已存到: {OUT}/partial_udl_varying.png")
