"""
Phase 3(支承強制位移)示範腳本 -- 對照SW FEA的兩個案例

跑法: PYTHONPATH=. python examples/phase3_support_displacement_demo.py

模型: 固接-固接門型鋼架, 節點0=(0,0)、節點1=(6,0)為基礎,
節點2=(0,4)、節點3=(6,4)為柱頂。E=200GPa, I=8.0E7mm^4(=8e-5m^4), A=0.01m^2。

案例A: 節點0被強制轉0.3 rad(沒有外加載重)
案例B: 節點0下沉10mm(轉角仍固定=0)

跟SW FEA的交叉驗證(反力+兩案例x三桿件x11個位置點的SF/BM/dX/dY)見
tests/test_support_displacement_vs_swfea.py。
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


def build(rot0, uy0):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')   # 左柱
    f.add_member(1, 2, 3, 's')   # 梁
    f.add_member(2, 3, 1, 's')   # 右柱
    f.support(0, ux=0.0, uy=uy0, rot=rot0)
    f.support(1, ux=0.0, uy=0.0, rot=0.0)
    return f


def print_result(f, r, title):
    print(f"\n{'='*60}\n{title}\n{'='*60}")
    for nid in [0, 1]:
        ux, uy, rot = f.dofs_of(nid)
        print(f"  node{nid}: Rx={r.reactions[ux]:.3f}  Ry={r.reactions[uy]:.3f}  M={r.reactions[rot]:.3f}")


# ---- 案例A: 強制轉角0.3 rad, 無外加載重 ----
fA = build(rot0=0.3, uy0=0.0)
rA = solve(fA)
print_result(fA, rA, "案例A: 支承強制轉角0.3rad")

figA = plot_all(fA, rA, figsize=(15, 9))
figA.suptitle('Phase3 vs SW FEA: forced rotation 0.3rad at node0 (no external load)', fontsize=11)
figA.tight_layout()
figA.savefig(f'{OUT}/phase3_rot_vs_swfea.png', dpi=130)
plt.close(figA)

# ---- 案例B: 支承沉陷10mm ----
fB = build(rot0=0.0, uy0=-0.01)
rB = solve(fB)
print_result(fB, rB, "案例B: 支承沉陷10mm")

figB = plot_all(fB, rB, figsize=(15, 9))
figB.suptitle('Phase3 vs SW FEA: support settlement 10mm at node0', fontsize=11)
figB.tight_layout()
figB.savefig(f'{OUT}/phase3_sinking_vs_swfea.png', dpi=130)
plt.close(figB)

print(f"\n六合一圖已存到: {OUT}/phase3_rot_vs_swfea.png, {OUT}/phase3_sinking_vs_swfea.png")
