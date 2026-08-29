"""
兩層樓門型鋼架, 三種內部鉸接情境示範 -- 對照SW FEA(pin2_0/pin2_1/pin2_2)

跑法: PYTHONPATH=. python examples/two_story_internal_pin_demo.py

驗證細節見 tests/test_two_story_release_vs_swfea.py。這支腳本只負責畫圖跟
印反力, 方便直接看六合一圖比較三種情境。
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from frame2d import Frame2D, solve

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

E, I, A = 200e6, 8e-5, 0.01
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}


def build(release_beam_at5=False, release_upper_right_col_at3=False):
    f = Frame2D()
    for nid, (x, y) in NODES.items():
        f.add_node(nid, x, y)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')                                          # F0: 1樓左柱
    f.add_member(1, 2, 3, 's')                                          # F1: 1樓梁
    f.add_member(2, 3, 1, 's')                                          # F2: 1樓右柱
    f.add_member(3, 4, 2, 's')                                          # F3: 2樓左柱
    f.add_member(4, 5, 4, 's', release_i=release_beam_at5)              # F4: 2樓梁
    f.add_member(5, 5, 3, 's', release_j=release_upper_right_col_at3)   # F5: 2樓右柱
    f.fix(0)
    f.fix(1)
    f.point_load(2, fx=10.0)
    f.point_load(4, fx=15.0)
    return f


from frame2d.plotting import plot_all

cases = [
    ('pin2_0_no_release', build(), 'pin2_0: no release (control, matches SW FEA)'),
    ('pin2_1_beam_at_node5', build(release_beam_at5=True),
     'pin2_1: beam released at node5 (frame2d self-consistent; SW FEA disagrees & self-inconsistent)'),
    ('pin2_2_col_at_node3', build(release_upper_right_col_at3=True),
     'pin2_2: upper-right column released at node3 (matches SW FEA)'),
]
for fname, f, title in cases:
    r = solve(f)
    fig = plot_all(f, r, figsize=(15, 9))
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{OUT}/{fname}.png', dpi=120)
    plt.close(fig)
    print('saved', fname)
    ux0, uy0, rot0 = f.dofs_of(0)
    ux1, uy1, rot1 = f.dofs_of(1)
    print(f'  node0: Rx={r.reactions[ux0]:.3f} Ry={r.reactions[uy0]:.3f} M={r.reactions[rot0]:.3f}')
    print(f'  node1: Rx={r.reactions[ux1]:.3f} Ry={r.reactions[uy1]:.3f} M={r.reactions[rot1]:.3f}')
