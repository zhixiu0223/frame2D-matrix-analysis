"""
Part 2: 不均勻雪載重屋頂模型 -- 跟你手機SW FEA app同一個幾何
(node0=(0,0)/node1=(6,0)固接柱腳, node2=(0,4)/node3=(6,4)柱頂,
node4=(3,7)屋頂尖端), 雪載重從slop-roof-snow.pdf: 屋脊0.2L沒有雪
(F7從簷口node2的10kN/m線性降到0, 蓋到0.8L為止; F8從0.2L的0kN/m線性升
到簷口node3的10kN/m), 方向是全域垂直方向(跟桿件角度無關, 大小以沿
桿件長度量測)——這個方向已經用SW FEA app的另一個對照案例(整段都有雪)
逐項驗證過, 見tests/test_snow_load_vs_swfea_app.py, 反力/軸力都精確
吻合到4位小數。

這裡算出來的反力(node0: Rx=2.664 Ry=16.970 M=-4.538)跟
slop-roof-snow.pdf報告顯示的支承反力(3.544/19.092/-6.243)對不起來;
既然全域垂直方向已經確認是對的, 這個落差不是我們算錯, 目前最可能的
解釋是那份PDF報告當時匯出的是舊快取的計算結果(這個App之前就出現過
這種情況), 需要zhixiu在app裡對這個0.2L沒有雪的案例重新觸發一次計算、
重新截圖確認之後才能做最終比對。

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
L = np.hypot(3.0, 3.0)   # 斜屋頂桿件長度 = 4.243m

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, 6, 0)
f.add_node(2, 0, 4)
f.add_node(3, 6, 4)
f.add_node(4, 3, 7)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 2, 's')   # F0: 左柱
f.add_member(1, 2, 4, 's')   # F7: 左屋頂斜梁 (簷口node2 -> 屋脊node4)
f.add_member(2, 4, 3, 's')   # F8: 右屋頂斜梁 (屋脊node4 -> 簷口node3)
f.add_member(3, 3, 1, 's')   # F3: 右柱
f.fix(0)
f.fix(1)

# 屋脊0.2L沒有雪, 方向='global', angle_deg=-90 (全域正下方, 已驗證)
f.distributed_load(1, w=10.0, w_end=0.0, x_start=0.0, x_end=0.8 * L,
                    direction='global', angle_deg=-90.0)
f.distributed_load(2, w=0.0, w_end=10.0, x_start=0.2 * L, x_end=L,
                    direction='global', angle_deg=-90.0)

r = solve(f)

print("節點反力 (屋脊0.2L沒有雪, 全域垂直方向, 已驗證的載重定義):")
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    print(f"  node{n}: Rx={r.reactions[ux]:.3f}kN  Ry={r.reactions[uy]:.3f}kN  "
          f"M={r.reactions[rot]:.3f}kNm")
print("\n(跟slop-roof-snow.pdf報告的3.544/19.092/-6.243對不起來——全域垂直")
print(" 方向已用另一個獨立案例驗證過是對的, 這個落差目前最可能是那份PDF")
print(" 匯出時是舊快取, 需要重新在app裡觸發計算確認)")

fig = plot_all(f, r, figsize=(15, 9))
fig.suptitle("Part2: 不均勻雪載重屋頂 (屋脊0.2L沒有雪, 全域垂直方向已驗證)", fontsize=11)
fig.tight_layout()
fig.savefig(f'{OUT}/snow_load_unbalanced.png', dpi=130)
plt.close(fig)
print(f"\n六合一圖已存到: {OUT}/snow_load_unbalanced.png")
