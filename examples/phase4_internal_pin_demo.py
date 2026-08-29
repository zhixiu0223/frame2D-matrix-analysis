"""
Phase 4(內部鉸接 / element release)示範腳本 -- 對照SW FEA

跑法: PYTHONPATH=. python examples/phase4_internal_pin_demo.py

模型: 固接-固接門型鋼架, 節點0=(0,0)、節點1=(6,0)為基礎,
節點2=(0,4)、節點3=(6,4)為柱頂。F0=左柱(剛接)、F1=梁、F2=右柱(剛接),
F1跟F2在節點3(右上角)是內部鉸接。節點2受水平力10kN向右。

跟SW FEA的交叉驗證(反力+三桿件x11個位置點BM)見
tests/test_element_release_vs_swfea.py。
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
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 2, 's')                     # F0: 左柱, 剛接
f.add_member(1, 2, 3, 's', release_j=True)      # F1: 梁, J端(node3)內部鉸接
f.add_member(2, 3, 1, 's')                      # F2: 右柱, 剛接
f.fix(0)
f.fix(1)
f.point_load(2, fx=10.0)   # 10kN向右
r = solve(f)

print("節點反力:")
for nid in [0, 1]:
    ux, uy, rot = f.dofs_of(nid)
    print(f"  node{nid}: Rx={r.reactions[ux]:.3f}  Ry={r.reactions[uy]:.3f}  M={r.reactions[rot]:.3f}")

print("\n鉸接端(F1的J端, node3)彎矩:", f"{r.member_results[1].end_forces_local[5]:.2e}", "(應精確為0)")

fig = plot_all(f, r, figsize=(15, 9))
fig.suptitle('Phase4 vs SW FEA: portal frame with internal pin at node3 (top-right)', fontsize=11)
fig.tight_layout()
fig.savefig(f'{OUT}/phase4_internal_pin_vs_swfea.png', dpi=130)
plt.close(fig)

print(f"\n六合一圖已存到: {OUT}/phase4_internal_pin_vs_swfea.png")
