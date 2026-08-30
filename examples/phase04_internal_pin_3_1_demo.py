"""
phase-04-internal-pin_3_1 示範腳本 -- 桿件內部點載重+局部段均佈載重+
兩個桿件力矩+內部鉸接, 全部組合在一起, 對照SW FEA

跑法: PYTHONPATH=. python examples/phase04_internal_pin_3_1_demo.py

模型: 單層門型鋼架, 節點0=(0,0)、節點1=(6,0)固接, 節點2=(0,4)、
節點3=(6,4)。F0=左柱、F1=梁(J端在節點3有內部鉸接)、F2=右柱。

載重:
  F1(梁): a=2.0m處10kN向下的桿件內部點載重
  F1(梁): [2.0,4.0]這段10kN/m向下的局部段均佈載重
  F1(梁): a=0(=節點2位置)處10kNm逆時針的桿件內部力矩
  F2(右柱): a=2.0m處10kNm逆時針的桿件內部力矩

**開發過程中的重要記錄**: fixed_end_forces_point_moment()這條公式原本有
正負號bug(推導時誤套用了點載重公式的取負號規律), 導致SW FEA UI標示的
「Counter-Clockwise」一度被誤判成跟我們的CCW正慣例相反——這個誤判已經
撤回並修正了公式本身, 現在m=+10直接對應逆時針, 不需要反號。詳細記錄見
tests/test_combined_release_loads_vs_swfea.py、
tests/test_member_point_load.py。

跟SW FEA的完整交叉驗證(反力+3桿件x11點BM逐點比對+桿件a=0力矩與真正
節點力矩的等價性驗證)見 tests/test_combined_release_loads_vs_swfea.py。
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
f.add_member(0, 0, 2, 's')                     # F0: 左柱
f.add_member(1, 2, 3, 's', release_j=True)     # F1: 梁, J端(節點3)內部鉸接
f.add_member(2, 3, 1, 's')                     # F2: 右柱
f.fix(0)
f.fix(1)
f.member_point_load(1, a=2.0, fy=-10.0)                  # F1: a=2.0m, 10kN向下
f.distributed_load(1, w=-10.0, x_start=2.0, x_end=4.0)   # F1: [2,4]這段, 10kN/m向下
f.member_point_load(2, a=2.0, m=10.0)                    # F2: a=2.0m, 10kNm逆時針
f.member_point_load(1, a=0.0, m=10.0)                    # F1: a=0(=節點2), 10kNm逆時針

r = solve(f)

print("節點反力:")
sw_reactions = {0: (3.765, 20.434, -2.457), 1: (-3.765, 9.566, 5.060)}
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    print(f"  node{n}: Rx={r.reactions[ux]:.3f}  Ry={r.reactions[uy]:.3f}  "
          f"M={r.reactions[rot]:.3f}  (SW FEA={sw_reactions[n]})")

fig = plot_all(f, r, figsize=(15, 9))
fig.suptitle('phase-04-internal-pin_3_1: point load + partial UDL + 2 member moments '
             '+ release, matches SW FEA', fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUT}/phase04_internal_pin_3_1.png', dpi=130)
plt.close(fig)

print(f"\n六合一圖已存到: {OUT}/phase04_internal_pin_3_1.png")
