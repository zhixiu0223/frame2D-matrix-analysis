"""
Part 2: 不均勻(不平衡)雪載重屋頂模型 -- 跟你手機SW FEA截圖同一個幾何
(node0=(0,0)/node1=(6,0)固接柱腳, node2=(0,4)/node3=(6,4)柱頂,
node4=(3,7)屋頂尖端), 但改成不均勻雪載重(而不是SW FEA截圖裡均佈的情況)。

雪載重用distributed_load(direction='global_y', w_start!=w_end)表示:
大小以沿桿件長度量測, 方向永遠垂直向下(不管屋頂本身斜不斜), 這是
工程上表示屋頂重力/雪載重的標準方式。這裡示範一個典型的"不平衡雪載重"
形狀(風吹雪在屋脊附近堆積、簷口較薄的簡化示意, 不是特定規範公式):
左屋頂簷口4kN/m -> 屋脊20kN/m, 右屋頂屋脊20kN/m -> 簷口8kN/m
(刻意左右不對稱, 呈現真正的"不均勻"載重, 不是單純對稱的兩側鏡射)。

這個功能是這次新增的(frame2d/elements.py新增
fixed_end_forces_axial_udl_varying(), frame2d/dofmanager.py更新global_y
分解邏輯支援兩端不同值), 獨立驗證見
tests/test_global_y_varying_snow_load.py(細網格分段模型收斂交叉驗證)。

跑法: PYTHONPATH=. python examples/snow_load_demo.py
"""
import os
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.sans-serif'] = ['Noto Sans CJK JP', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
import numpy as np
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

E, I, A = 200e6, 8e-5, 0.01

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, 6, 0)
f.add_node(2, 0, 4)
f.add_node(3, 6, 4)
f.add_node(4, 3, 7)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 2, 's')   # F0: 左柱
f.add_member(1, 2, 4, 's')   # F1: 左屋頂斜梁 (簷口node2 -> 屋脊node4)
f.add_member(2, 4, 3, 's')   # F2: 右屋頂斜梁 (屋脊node4 -> 簷口node3)
f.add_member(3, 3, 1, 's')   # F3: 右柱
f.fix(0)
f.fix(1)

W_LEFT_EAVE, W_LEFT_RIDGE = 4.0, 20.0
W_RIGHT_RIDGE, W_RIGHT_EAVE = 20.0, 8.0
f.distributed_load(1, w=W_LEFT_EAVE, w_end=W_LEFT_RIDGE, direction='global_y')
f.distributed_load(2, w=W_RIGHT_RIDGE, w_end=W_RIGHT_EAVE, direction='global_y')

r = solve(f)

print("節點反力 (不均勻雪載重: 左簷口4→屋脊20 kN/m, 右屋脊20→簷口8 kN/m):")
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    print(f"  node{n}: Rx={r.reactions[ux]:.3f}kN  Ry={r.reactions[uy]:.3f}kN  "
          f"M={r.reactions[rot]:.3f}kNm")

# 整體靜力平衡檢查(必要條件, 跟載重多不均勻無關恆成立)
L_slope = np.sqrt(3.0 ** 2 + 3.0 ** 2)
W_total_exact = 0.5 * (W_LEFT_EAVE + W_LEFT_RIDGE) * L_slope \
    + 0.5 * (W_RIGHT_RIDGE + W_RIGHT_EAVE) * L_slope
ux0, uy0, _ = f.dofs_of(0)
ux1, uy1, _ = f.dofs_of(1)
Ry_total = r.reactions[uy0] + r.reactions[uy1]
Rx_total = r.reactions[ux0] + r.reactions[ux1]
print(f"\n總垂直反力: {Ry_total:.3f}kN  (總雪載重解析值: {W_total_exact:.3f}kN, "
      f"差={abs(Ry_total - W_total_exact):.2e})")
print(f"總水平反力: {Rx_total:.2e}kN (應接近0)")

print(f"\n屋頂尖端(node4)全域垂直位移: {r.displacements[f.dofs_of(4)[1]] * 1000:.3f}mm")

fig = plot_all(f, r, figsize=(15, 9))
fig.suptitle("Part2: 不均勻雪載重屋頂 (簷口輕→屋脊重的不平衡分布)", fontsize=11)
fig.tight_layout()
fig.savefig(f'{OUT}/snow_load_unbalanced.png', dpi=130)
plt.close(fig)
print(f"\n六合一圖已存到: {OUT}/snow_load_unbalanced.png")
