"""
Internal-Hinge Benchmark Suite -- 幾何結構總覽圖

畫出兩層樓門型鋼架, 六種不同的內部鉸接配置(pin2_0~pin2_5), 用有意義的
分割長度(0.8~1.2m)示意鉸接實際位置, 方便目視確認幾何正確。

跑法: PYTHONPATH=. python examples/pin_configs_overview_demo.py

對應的完整數值驗證(反力+逐點BM比對)見:
  tests/test_meaningful_length_split_vs_swfea.py  (有意義長度版本)
  tests/test_tiny_offset_reproduction.py          (SW FEA原始極小offset版本)
  tests/test_triple_release_swfea_bug.py          (pin2_4/pin2_5數學等價性證明)
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from frame2d import Frame2D
from frame2d.plotting import plot_structure

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

E, I, A = 200e6, 8e-5, 0.01
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}


def base(f):
    """六個案例共用的節點/斷面/支承(2樓門型鋼架, 節點0,1固接)。"""
    for nid, (x, y) in NODES.items():
        f.add_node(nid, x, y)
    f.add_section('s', E=E, I=I, A=A)
    f.fix(0)
    f.fix(1)


# 每個案例用什麼模式建立(None=無鉸接; 其餘見下方if/elif分支)
CONFIGS = [
    ('pin2_0: no hinge (control)', None),
    ('pin2_1: F6 (1F beam) hinge 1.2m from N3', 'F6_1.2'),
    ('pin2_2: F5 (2F right col) hinge 0.8m from N3', 'F5_0.8'),
    ('pin2_3: F2 (1F right col) hinge 0.8m from N3', 'F2_0.8'),
    ('pin2_4: F2+F6+F5 all released at N3', 'ALL'),
    ('pin2_5: only F2+F5 (beam rigid)', 'F2F5'),
]


def build(mode):
    """依mode建立對應的frame2d模型。有鉸接的案例都用"雙桿件分割"表示
    (在鉸接位置插入真實節點, 拆成兩段, 其中一段的分割點端release), 這樣
    畫圖時鉸接位置(紅色空心圓)才會落在桿件實際位置上, 不是疊在節點正上方。"""
    f = Frame2D()
    base(f)
    if mode is None:
        f.add_member(0, 0, 2, 's')
        f.add_member(1, 2, 3, 's')
        f.add_member(2, 3, 1, 's')
        f.add_member(3, 4, 2, 's')
        f.add_member(4, 5, 4, 's')
        f.add_member(5, 5, 3, 's')
    elif mode == 'F6_1.2':
        f.add_node(6, 4.8, 4)   # F6上, 距節點3為1.2m處(F6從節點2(0,4)到節點3(6,4))
        f.add_member(0, 0, 2, 's')
        f.add_member(1, 2, 6, 's', release_j=True)
        f.add_member(9, 6, 3, 's')
        f.add_member(2, 3, 1, 's')
        f.add_member(3, 4, 2, 's')
        f.add_member(4, 5, 4, 's')
        f.add_member(5, 5, 3, 's')
    elif mode == 'F5_0.8':
        f.add_node(7, 6, 4.8)   # F5上, 距節點3為0.8m處(F5從節點5(6,8)到節點3(6,4))
        f.add_member(0, 0, 2, 's')
        f.add_member(1, 2, 3, 's')
        f.add_member(2, 3, 1, 's')
        f.add_member(3, 4, 2, 's')
        f.add_member(4, 5, 4, 's')
        f.add_member(5, 5, 7, 's', release_j=True)
        f.add_member(10, 7, 3, 's')
    elif mode == 'F2_0.8':
        f.add_node(8, 6, 3.2)   # F2上, 距節點3為0.8m處(F2從節點3(6,4)到節點1(6,0))
        f.add_member(0, 0, 2, 's')
        f.add_member(1, 2, 3, 's')
        f.add_member(2, 3, 8, 's')
        f.add_member(9, 8, 1, 's', release_i=True)
        f.add_member(3, 4, 2, 's')
        f.add_member(4, 5, 4, 's')
        f.add_member(5, 5, 3, 's')
    elif mode == 'ALL':
        f.add_node(6, 4.8, 4)
        f.add_node(7, 6, 4.8)
        f.add_node(8, 6, 3.2)
        f.add_member(0, 0, 2, 's')
        f.add_member(1, 2, 6, 's', release_j=True)
        f.add_member(9, 6, 3, 's')
        f.add_member(2, 3, 8, 's')
        f.add_member(11, 8, 1, 's', release_i=True)
        f.add_member(3, 4, 2, 's')
        f.add_member(4, 5, 4, 's')
        f.add_member(5, 5, 7, 's', release_j=True)
        f.add_member(10, 7, 3, 's')
    elif mode == 'F2F5':
        f.add_node(7, 6, 4.8)
        f.add_node(8, 6, 3.2)
        f.add_member(0, 0, 2, 's')
        f.add_member(1, 2, 3, 's')
        f.add_member(2, 3, 8, 's')
        f.add_member(9, 8, 1, 's', release_i=True)
        f.add_member(3, 4, 2, 's')
        f.add_member(4, 5, 4, 's')
        f.add_member(5, 5, 7, 's', release_j=True)
        f.add_member(10, 7, 3, 's')
    return f


fig, axes = plt.subplots(2, 3, figsize=(16, 11))
for ax, (title, mode) in zip(axes.flat, CONFIGS):
    f = build(mode)
    plot_structure(f, ax=ax, show_member_ids=False, show_node_ids=True)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-1.5, 7.5)
    ax.set_ylim(-1, 9)

fig.suptitle(
    'Internal-Hinge Benchmark Suite: all tested two-story frame hinge configurations\n'
    '(red circle = hinge location; shown here at meaningful offsets 0.8-1.2m; '
    'exact-node cases also validated separately)',
    fontsize=11,
)
fig.tight_layout()
fig.savefig(f'{OUT}/all_pin_configs_overview.png', dpi=130)
plt.close(fig)

print(f"已存圖: {OUT}/all_pin_configs_overview.png")
