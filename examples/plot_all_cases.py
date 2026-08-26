import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all

EI = dict(E=1.0, I=15000.0, A=1e8)
OUT = '/home/claude/frame2d_plots'


def save(f, title, fname):
    r = solve(f)
    fig = plot_all(f, r, figsize=(15, 9))
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(f'{OUT}/{fname}.png', dpi=120)
    plt.close(fig)
    print('saved', fname)


# ---- Case-01: propped cantilever ----
L, w = 6.0, 20.0
f = Frame2D()
f.add_node(0, 0, 0); f.add_node(1, L, 0)
f.add_section('s', **EI)
f.add_member(0, 0, 1, 's')
f.fix(0); f.roller_y(1)
f.distributed_load(0, w=-w)
save(f, 'Case-01: Propped Cantilever (L=6, w=20)', 'case01')

# ---- Case-02: two-span beam ----
L1, L2, w1, w2 = 5.0, 6.0, 15.0, 20.0
f = Frame2D()
f.add_node(0, 0, 0); f.add_node(1, L1, 0); f.add_node(2, L1 + L2, 0)
f.add_section('s', **EI)
f.add_member(0, 0, 1, 's'); f.add_member(1, 1, 2, 's')
f.fix(0); f.roller_y(1); f.roller_y(2)
f.distributed_load(0, w=-w1); f.distributed_load(1, w=-w2)
save(f, 'Case-02: Two-Span Beam (L1=5,L2=6, w1=15,w2=20)', 'case02')

# ---- Case-03: no-sway frame ----
H, L, w = 4.0, 6.0, 24.0
f = Frame2D()
f.add_node(0, 0, 0); f.add_node(1, 0, H); f.add_node(2, L, H); f.add_node(3, L, 0)
f.add_section('s', **EI)
f.add_member(0, 0, 1, 's'); f.add_member(1, 1, 2, 's'); f.add_member(2, 3, 2, 's')
f.fix(0); f.fix(3)
f.distributed_load(1, w=-w)
save(f, 'Case-03: No-Sway Frame (H=4,L=6, w=24)', 'case03')


def build_sway(H, L, P=0.0, w=0.0):
    f = Frame2D()
    f.add_node(0, 0, 0); f.add_node(1, 0, H); f.add_node(2, L, H); f.add_node(3, L, 0)
    f.add_section('s', **EI)
    f.add_member(0, 0, 1, 's'); f.add_member(1, 1, 2, 's'); f.add_member(2, 3, 2, 's')
    f.fix(0); f.fix(3)
    if P: f.point_load(1, fx=P)
    if w: f.distributed_load(1, w=-w)
    return f


# ---- Case-04 ----
f = build_sway(4.0, 6.0, P=12.0)
save(f, 'Case-04: Sway Frame, no UDL (H=4,L=6,P=12)', 'case04')

# ---- Case-04.5 ----
f = build_sway(4.0, 6.0, P=12.0, w=24.0)
save(f, 'Case-04.5: Sway Frame + UDL (H=4,L=6,P=12,w=24)', 'case04_5')


def build_two_story(H1, H2, L, w1=0.0, w2=0.0, P1=0.0, P2=0.0):
    Htot = H1 + H2
    f = Frame2D()
    f.add_node(0, 0, 0); f.add_node(1, 0, H1); f.add_node(2, 0, Htot)
    f.add_node(3, L, Htot); f.add_node(4, L, H1); f.add_node(5, L, 0)
    f.add_section('s', **EI)
    f.add_member(0, 0, 1, 's'); f.add_member(1, 1, 2, 's'); f.add_member(2, 2, 3, 's')
    f.add_member(3, 3, 4, 's'); f.add_member(4, 4, 5, 's'); f.add_member(5, 1, 4, 's')
    f.fix(0); f.fix(5)
    if w1: f.distributed_load(5, w=-w1)
    if w2: f.distributed_load(2, w=-w2)
    if P1: f.point_load(1, fx=P1)
    if P2: f.point_load(2, fx=P2)
    return f


# ---- Case-05 ----
f = build_two_story(4.0, 3.5, 6.0, w1=24.0, w2=18.0)
save(f, 'Case-05: Two-Story Frame, no lateral load (w1=24,w2=18)', 'case05')

# ---- Case-06 ----
f = build_two_story(4.0, 3.5, 6.0, P1=15.0, P2=10.0)
save(f, 'Case-06: Two-Story Sway, no UDL (P1=15,P2=10)', 'case06')

# ---- Case-06.5 ----
f = build_two_story(4.0, 3.5, 6.0, w1=24.0, w2=18.0, P1=15.0, P2=10.0)
save(f, 'Case-06.5: Two-Story Sway + UDL (P1=15,P2=10,w1=24,w2=18)', 'case06_5')


def build_two_bay(H, L1, L2, w1=0.0, w2=0.0, P=0.0):
    f = Frame2D()
    f.add_node(0, 0, 0); f.add_node(1, 0, H); f.add_node(2, L1, H)
    f.add_node(3, L1, 0); f.add_node(4, L1 + L2, H); f.add_node(5, L1 + L2, 0)
    f.add_section('s', **EI)
    f.add_member(0, 0, 1, 's'); f.add_member(1, 3, 2, 's'); f.add_member(2, 5, 4, 's')
    f.add_member(3, 1, 2, 's'); f.add_member(4, 2, 4, 's')
    f.fix(0); f.fix(3); f.fix(5)
    if w1: f.distributed_load(3, w=-w1)
    if w2: f.distributed_load(4, w=-w2)
    if P: f.point_load(1, fx=P)
    return f


# ---- Case-07 ----
f = build_two_bay(4.0, 5.0, 7.0, w1=20.0, w2=15.0)
save(f, 'Case-07: Two-Bay Frame, no lateral load (w1=20,w2=15)', 'case07')

# ---- Case-07.5 ----
f = build_two_bay(4.0, 5.0, 7.0, w1=20.0, w2=15.0, P=10.0)
save(f, 'Case-07.5: Two-Bay Frame + P (w1=20,w2=15,P=10)', 'case07_5')

# ---- Case-08 (無UDL, 純點載重版本, 對照sd_framework預設值) ----
f = Frame2D()
nodes = {0: (0, 0), 1: (5, 0), 2: (12, 0),
         3: (0, 4), 4: (5, 4), 5: (12, 4),
         6: (0, 7.5), 7: (5, 7.5), 8: (12, 7.5)}
for nid, (x, y) in nodes.items():
    f.add_node(nid, x, y)
f.add_section('s', **EI)
members = [(0, 0, 3), (1, 1, 4), (2, 5, 2), (3, 3, 4), (4, 4, 5),
           (5, 6, 3), (6, 7, 4), (7, 8, 5), (8, 6, 7), (9, 7, 8)]
for mid, ni, nj in members:
    f.add_member(mid, ni, nj, 's')
for n in [0, 1, 2]:
    f.fix(n)
f.point_load(6, fx=10.0)
f.point_load(3, fx=15.0)
save(f, 'Case-08: Two-Story Two-Bay (P1=15,P2=10, no UDL)', 'case08')

# ---- Case-08 UDL版 (對照SW FEA報告那組, w1=10,w2=12) ----
f2 = Frame2D()
for nid, (x, y) in nodes.items():
    f2.add_node(nid, x, y)
f2.add_section('s', E=150e9, I=1e-6, A=0.05)
for mid, ni, nj in members:
    f2.add_member(mid, ni, nj, 's')
for n in [0, 1, 2]:
    f2.fix(n)
f2.point_load(6, fx=-10000)
f2.point_load(3, fx=-15000)
f2.distributed_load(8, w=12000)
f2.distributed_load(9, w=12000)
f2.distributed_load(3, w=10000)
f2.distributed_load(4, w=10000)
save(f2, 'Case-08 (UDL): Two-Story Two-Bay vs SW FEA (P1=15kN,P2=10kN,w1=10kN/m,w2=12kN/m)', 'case08_udl')

print('全部完成')
