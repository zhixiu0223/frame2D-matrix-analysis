"""
Load System v2 Phase 1 + Phase 2 示範腳本

跑法: PYTHONPATH=. python examples/load_system_v2_demo.py
(或裝成套件後直接 python examples/load_system_v2_demo.py)

會印出完整反力/桿端內力數值(方便直接抓去跟SW FEA或其他工具比對),
並且輸出六合一圖到 examples/output_plots/。
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all
from frame2d.postprocess import member_internal_forces

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

E, I, A = 200e6, 8e-5, 1e-2   # 任意單位, 只是demo用的斷面數值


def print_result(f, r, title):
    print(f"\n{'='*60}\n{title}\n{'='*60}")
    print("節點反力 (Rx, Ry, M):")
    for nid in f.nodes:
        ux, uy, rot = f.dofs_of(nid)
        Rx, Ry, M = r.reactions[ux], r.reactions[uy], r.reactions[rot]
        if abs(Rx) > 1e-9 or abs(Ry) > 1e-9 or abs(M) > 1e-9:
            print(f"  node{nid}: Rx={Rx:10.4f}  Ry={Ry:10.4f}  M={M:10.4f}")
    print("\n桿件端點內力 (N, V, M, 局部座標, index0=node_i端, index-1=node_j端):")
    for mid in f.members:
        x, N, V, M = member_internal_forces(f, r, mid, n=2)
        print(f"  member{mid}: node_i端 N={N[0]:9.4f} V={V[0]:9.4f} M={M[0]:9.4f}"
              f"   |   node_j端 N={N[-1]:9.4f} V={V[-1]:9.4f} M={M[-1]:9.4f}")


# ==================== Phase 1: 桿件中間的集中力 ====================
# 簡支梁, 跨度10m, 節點0=(0,0)鉸接, 節點1=(10,0)滾支承,
# 在距節點0 3.5m處(桿件內部, 不在節點上)加垂直向下30kN的集中力。
L, P, a = 10.0, 30.0, 3.5
f1 = Frame2D()
f1.add_node(0, 0, 0)
f1.add_node(1, L, 0)
f1.add_section('s', E=E, I=I, A=A)
f1.add_member(0, 0, 1, 's')
f1.pin(0)
f1.roller_y(1)
f1.member_point_load(0, a=a, fy=-P)
r1 = solve(f1)
print_result(f1, r1, "Phase 1 demo: 簡支梁, 桿件中間點載重 (L=10, P=30 向下, a=3.5)")

fig1 = plot_all(f1, r1, figsize=(15, 9))
fig1.suptitle('Phase1 demo: simply supported beam, interior point load at a=3.5', fontsize=11)
fig1.tight_layout()
fig1.savefig(f'{OUT}/phase1_point_load_demo.png', dpi=120)
plt.close(fig1)

# 注意: 上面六合一圖裡「變形圖」的 Δmax 標籤, 是用單一桿件的Hermite內插算的,
# 桿件內部有點載重時內插撓度會低估實際值(見postprocess.py的說明)。
# 如果要精確的撓度數值(例如要跟SW FEA這類逐點計算的工具比對), 建議直接在
# 載重點插入一個真實節點, 把桿件拆成兩段 -- 該點的位移就是FEM節點解, 沒有
# 內插誤差。這裡示範精確算法:
f1_exact = Frame2D()
f1_exact.add_node(0, 0, 0)
f1_exact.add_node(1, a, 0)     # 直接在載重位置a放節點
f1_exact.add_node(2, L, 0)
f1_exact.add_section('s', E=E, I=I, A=A)
f1_exact.add_member(0, 0, 1, 's')
f1_exact.add_member(1, 1, 2, 's')
f1_exact.pin(0)
f1_exact.roller_y(2)
f1_exact.point_load(1, fy=-P)   # 載重直接當節點載重, 不用member_point_load
r1_exact = solve(f1_exact)
d_exact = r1_exact.displacements[f1_exact.dofs_of(1)[1]]
print(f"\n[精確撓度] 在載重點x={a}處插入節點求解(沒有Hermite內插誤差): {d_exact*1000:.4f} mm")
print("  (這個數字才適合拿去跟SW FEA逐點計算的撓度值比對, 不是六合一圖裡的Δmax)")


# ==================== Phase 2: 局部段均佈載重 ====================
# 簡支梁, 跨度12m, 均佈載重20kN/m只加在[3,7]這一段(不是整根梁都有)。
L2, w, c, d = 12.0, -20.0, 3.0, 7.0
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L2, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's')
f2.pin(0)
f2.roller_y(1)
f2.distributed_load(0, w=w, x_start=c, x_end=d)
r2 = solve(f2)
print_result(f2, r2, "Phase 2 demo: 簡支梁, 局部段均佈載重 (L=12, w=-20 在[3,7])")

fig2 = plot_all(f2, r2, figsize=(15, 9))
fig2.suptitle('Phase2 demo: simply supported beam, partial UDL on [3,7] of a 12m span', fontsize=11)
fig2.tight_layout()
fig2.savefig(f'{OUT}/phase2_partial_udl_demo.png', dpi=120)
plt.close(fig2)

print(f"\n六合一圖已存到: {OUT}/phase1_point_load_demo.png, {OUT}/phase2_partial_udl_demo.png")
