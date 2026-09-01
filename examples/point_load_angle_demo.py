"""
point_load / member_point_load 角度便利介面 -- 視覺化示範

跑法: PYTHONPATH=. python examples/point_load_angle_demo.py

驗證邏輯(退化案例對照, 對應tests/test_point_load_angle.py):
在斜屋頂的F7桿件(45度斜梁)跨中(a=2.0)加一個全域垂直向下10kN的集中力,
分別用「手動分解成局部fx/fy」跟「新介面direction='global'自動分解」
兩種寫法建模, 印出反力確認完全一致, 再畫出六合一圖。
"""
import math
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

E, I, A = 200e6, 8e-5, 0.01


def base_roof():
    """單層三角形屋架, 節點0/1固接柱腳, 節點4是屋頂尖端。F7(member1)是
    node2->node4的45度斜梁。"""
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_node(4, 3, 7)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 4, 's')   # F7: 45度斜梁
    f.add_member(2, 4, 3, 's')
    f.add_member(3, 3, 1, 's')
    f.fix(0)
    f.fix(1)
    return f


# ---- 兩種寫法比對: 手動分解 vs 新介面自動分解 ----
F, angle_deg, a = 10.0, -90.0, 2.0
angle_F7 = math.radians(45)   # F7桿件自己的角度
ang = math.radians(angle_deg)
gx, gy = F * math.cos(ang), F * math.sin(ang)
c, s = math.cos(angle_F7), math.sin(angle_F7)
fx_local_manual = c * gx + s * gy
fy_local_manual = -s * gx + c * gy

f_manual = base_roof()
f_manual.member_point_load(1, a=a, fx=fx_local_manual, fy=fy_local_manual)
r_manual = solve(f_manual)

f_new = base_roof()
f_new.member_point_load(1, a=a, direction='global', F=F, angle_deg=angle_deg)
r_new = solve(f_new)

print("手動分解局部分量: fx={:.4f}, fy={:.4f}".format(fx_local_manual, fy_local_manual))
print()
print("反力比對(node0):")
for label, r, frame in [("手動分解(direction='local')", r_manual, f_manual),
                         ("新介面(direction='global')", r_new, f_new)]:
    ux, uy, rot = frame.dofs_of(0)
    print(f"  {label}: Rx={r.reactions[ux]:.6f}  Ry={r.reactions[uy]:.6f}  M={r.reactions[rot]:.6f}")

diff = max(abs(r_manual.reactions[i] - r_new.reactions[i]) for i in f_manual.dofs_of(0))
print(f"\n兩者差異: {diff:.2e} (應為0, 純浮點精度)")

# ---- 畫圖: 用新介面的版本(更直觀, 不用先手算局部分量) ----
fig = plot_all(f_new, r_new, figsize=(15, 9))
fig.suptitle("member_point_load(direction='global', F=10, angle_deg=-90): "
             "vertical point load on sloped member", fontsize=10)
fig.tight_layout()
fig.savefig(f'{OUT}/point_load_angle_demo.png', dpi=130)
plt.close(fig)

print(f"\n六合一圖已存到: {OUT}/point_load_angle_demo.png")
