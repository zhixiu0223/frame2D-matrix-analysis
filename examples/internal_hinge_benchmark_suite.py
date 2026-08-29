"""
Internal-Hinge Benchmark Suite (PIN2-0 ~ PIN2-6)

兩層樓門型鋼架, 節點3(1樓右側樓板)是三根桿件(F1=1樓梁、F2=1樓右柱、
F5=2樓右柱)的交界處。這組案例涵蓋內部鉸接的各種組合, 不是單純模仿
SW FEA的3個檔案, 是自己建立的完整測試家族:

  PIN2-0: 無鉸接(控制組, 對照SW FEA)
  PIN2-1: 2樓梁(F4)在節點5端釋放(對照SW FEA -- 但SW FEA這個案例本身
          的報告有內部矛盾, 見tests/test_two_story_release_vs_swfea.py
          的說明, 這裡取frame2d自己驗證過的結果)
  PIN2-2: 2樓右柱(F5)在節點3端釋放(對照SW FEA)
  PIN2-3: 1樓右柱(F2)在節點3端釋放(跟PIN2-2是節點3的另一根桿件, 沒有
          對應的SW FEA案例, 純frame2d自己的驗證)
  PIN2-4: 節點3的三根桿件(F1,F2,F5)全部釋放(chatGPT建議的benchmark,
          驗證"三向鉸接"不會變成機構, 整體力平衡成立)
  PIN2-5: 內部鉸接節點上同時施加節點外力(fy, m), 驗證M=0(鉸接的定義)
          跟Fx,Fy可以非零互不衝突
  PIN2-6: 桿件內部集中力剛好加在鉸接位置(這個組合靜力凝縮版本目前不
          支援, 但主要求解器solve()=DOFManager可以直接處理)

跑法: PYTHONPATH=. python examples/internal_hinge_benchmark_suite.py
每個案例都會印反力+畫六合一圖存到output_plots/。
逐點跟SW FEA的詳細比對見 tests/test_two_story_release_vs_swfea.py。
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
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}


def base_frame():
    """PIN2-0~PIN2-4共用的基礎兩層樓門型鋼架(節點/斷面/支承/載重都一樣,
    差別只在哪根桿件的哪一端有release)"""
    f = Frame2D()
    for nid, (x, y) in NODES.items():
        f.add_node(nid, x, y)
    f.add_section('s', E=E, I=I, A=A)
    f.fix(0)
    f.fix(1)
    f.point_load(2, fx=10.0)
    f.point_load(4, fx=15.0)
    return f


def report(name, f, r, title, fname):
    print(f"\n=== {name} ===")
    for n in [0, 1]:
        ux, uy, rot = f.dofs_of(n)
        print(f"  node{n}: Rx={r.reactions[ux]:.3f}  Ry={r.reactions[uy]:.3f}  M={r.reactions[rot]:.3f}")
    fig = plot_all(f, r, figsize=(15, 9))
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{OUT}/{fname}.png', dpi=120)
    plt.close(fig)
    print(f"  已存圖: {fname}.png")


# ---- PIN2-0: 無鉸接 ----
f0 = base_frame()
f0.add_member(0, 0, 2, 's')
f0.add_member(1, 2, 3, 's')
f0.add_member(2, 3, 1, 's')
f0.add_member(3, 4, 2, 's')
f0.add_member(4, 5, 4, 's')
f0.add_member(5, 5, 3, 's')
report("PIN2-0: 無鉸接(控制組)", f0, solve(f0),
       'PIN2-0: no internal hinge (control, matches SW FEA)', 'PIN2-0_no_release')

# ---- PIN2-1: 2樓梁(F4)在節點5端釋放 ----
f1 = base_frame()
f1.add_member(0, 0, 2, 's')
f1.add_member(1, 2, 3, 's')
f1.add_member(2, 3, 1, 's')
f1.add_member(3, 4, 2, 's')
f1.add_member(4, 5, 4, 's', release_i=True)
f1.add_member(5, 5, 3, 's')
report("PIN2-1: 2樓梁在節點5端釋放", f1, solve(f1),
       'PIN2-1: beam released at node5', 'PIN2-1_beam_at_node5')

# ---- PIN2-2: 2樓右柱(F5)在節點3端釋放 ----
f2 = base_frame()
f2.add_member(0, 0, 2, 's')
f2.add_member(1, 2, 3, 's')
f2.add_member(2, 3, 1, 's')
f2.add_member(3, 4, 2, 's')
f2.add_member(4, 5, 4, 's')
f2.add_member(5, 5, 3, 's', release_j=True)
report("PIN2-2: 2樓右柱在節點3端釋放", f2, solve(f2),
       'PIN2-2: upper column released at node3', 'PIN2-2_upper_col_at_node3')

# ---- PIN2-3: 1樓右柱(F2)在節點3端釋放(節點3的"另一根"桿件) ----
f3 = base_frame()
f3.add_member(0, 0, 2, 's')
f3.add_member(1, 2, 3, 's')
f3.add_member(2, 3, 1, 's', release_i=True)
f3.add_member(3, 4, 2, 's')
f3.add_member(4, 5, 4, 's')
f3.add_member(5, 5, 3, 's')
report("PIN2-3: 1樓右柱在節點3端釋放", f3, solve(f3),
       'PIN2-3: lower column released at node3', 'PIN2-3_lower_col_at_node3')

# ---- PIN2-4: 節點3三根桿件全部釋放 ----
f4 = base_frame()
f4.add_member(0, 0, 2, 's')
f4.add_member(1, 2, 3, 's', release_j=True)
f4.add_member(2, 3, 1, 's', release_i=True)
f4.add_member(3, 4, 2, 's')
f4.add_member(4, 5, 4, 's')
f4.add_member(5, 5, 3, 's', release_j=True)
r4 = solve(f4)
report("PIN2-4: 節點3三根桿件全部釋放", f4, r4,
       'PIN2-4: all 3 members released at node3 (still stable, not a mechanism)', 'PIN2-4_all_released')
for mid, end_idx, label in [(1, 5, "梁J端"), (2, 2, "1樓右柱I端"), (5, 5, "2樓右柱J端")]:
    print(f"    {label}彎矩: {r4.member_results[mid].end_forces_local[end_idx]:.2e} (應精確為0)")

# ---- PIN2-5: 鉸接節點上同時加節點外力 ----
f5 = base_frame()
f5.add_member(0, 0, 2, 's')
f5.add_member(1, 2, 3, 's')
f5.add_member(2, 3, 1, 's')
f5.add_member(3, 4, 2, 's')
f5.add_member(4, 5, 4, 's', release_i=True)
f5.add_member(5, 5, 3, 's')
f5.point_load(5, fy=-8.0, m=3.0)   # 鉸接端節點5額外加垂直力+節點力矩
r5 = solve(f5)
report("PIN2-5: 鉸接節點(node5)上額外施加fy=-8,m=3", f5, r5,
       'PIN2-5: nodal load (fy,m) applied directly at the hinge node', 'PIN2-5_hinge_with_nodal_load')
print(f"    鉸接端(F4的I端)彎矩: {r5.member_results[4].end_forces_local[2]:.2e} (應精確為0, 不受節點外力影響)")

# ---- PIN2-6: 桿件內部集中力加在鉸接附近(但不是剛好在鉸接點上, 這樣
#      才看得出"載重在桿件中間"這個情境, 不是退化成直接加在支承上) ----
f6 = Frame2D()
f6.add_node(0, 0, 0)
f6.add_node(1, 8, 0)
f6.add_section('s', E=E, I=I, A=A)
f6.add_member(0, 0, 1, 's', release_j=True)
f6.fix(0)
f6.pin(1)
f6.member_point_load(0, a=6.0, fy=-10.0)   # 桿件中段(不是端點), J端(node1)仍是鉸接
r6 = solve(f6)
report("PIN2-6: 桿件中段集中力, J端有鉸接", f6, r6,
       'PIN2-6: member point load at midspan, J-end released', 'PIN2-6_point_load_near_hinge')
print(f"    鉸接端彎矩: {r6.member_results[0].end_forces_local[5]:.2e} (應精確為0)")
print("    (這個組合solve_condensation()會報錯, 只有主要求解器solve()=DOFManager能直接處理)")

print(f"\n全部七張圖已存到: {OUT}/")
