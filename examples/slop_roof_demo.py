"""
斜屋頂示範腳本 -- distributed_load(direction='global_y'), 對照SW FEA

跑法: PYTHONPATH=. python examples/slop_roof_demo.py

模型: 單層三角形屋架, 節點0=(0,0)、節點1=(6,0)固接, 節點2=(0,4)、
節點3=(6,4)為柱頂, 節點4=(3,7)為屋頂尖端。F0=左柱、F7=左屋頂斜梁、
F8=右屋頂斜梁、F2=右柱。

屋頂載重10kN/m用direction='global_y': 大小以沿桿件長度量測, 方向永遠
垂直向下(不管桿件本身斜不斜)——這是工程上表示屋頂重力/雪載重的標準
方式, 不是"垂直於斜屋頂表面"那種局部座標方向。

跟SW FEA的完整交叉驗證(反力+4桿件x11點N/V/M逐點比對)見
tests/test_sloped_roof_global_udl.py。
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
f.add_member(1, 2, 4, 's')   # F7: 左屋頂斜梁
f.add_member(2, 4, 3, 's')   # F8: 右屋頂斜梁
f.add_member(3, 3, 1, 's')   # F2: 右柱
f.fix(0)
f.fix(1)
f.distributed_load(1, w=10.0, direction='global_y')
f.distributed_load(2, w=10.0, direction='global_y')

r = solve(f)

print("節點反力:")
sw_reactions = {0: (11.047, 42.426, -20.157), 1: (-11.047, 42.426, 20.157)}
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    print(f"  node{n}: Rx={r.reactions[ux]:.3f}  Ry={r.reactions[uy]:.3f}  "
          f"M={r.reactions[rot]:.3f}  (SW FEA={sw_reactions[n]})")

print(f"\n屋頂尖端(node4)全域垂直位移: {r.displacements[f.dofs_of(4)[1]]*1000:.3f}mm")
print("(SW FEA報告: -2.867mm)")

fig = plot_all(f, r, figsize=(15, 9))
fig.suptitle("slop-roof: distributed_load(direction='global_y'), matches SW FEA",
             fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUT}/slop_roof.png', dpi=130)
plt.close(fig)

print(f"\n六合一圖已存到: {OUT}/slop_roof.png")
