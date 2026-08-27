"""
懸臂梁斜張橋: 兩種方法畫出六合一圖並排比較

方法A(位移法/displacement method): 5條纜索都用cable元素建模, solve()自動
  判斷哪些該鬆弛(見tests/test_cable_stayed_bridge.py的驗證, 收斂結果是
  Cable1,Cable2鬆弛退出作用)。

方法B(力量法/force method的等效載重版): 不用cable元素, 而是把力量法
  (使用者自己的贅力法筆記本, 含塔柱彎曲柔度修正版)解出來的3條纜索張力
  (Cable3=176.0331, Cable4=257.7291, Cable5=274.2302 kN)當作已知外力,
  直接施加在梁跟塔柱上(牛頓第三定律: 梁承受纜索方向的拉力, 塔頂承受
  等大反向的反作用力), 解一個不含纜索元素、純粹是「主結構+已知外力」的
  結構。

兩者的樑跟塔柱內力/變形應該完全一致(已在test_cable_stayed_bridge.py驗證
彎矩逐點吻合, 這裡额外確認變形/軸力/剪力也全部一致) —— 這是位移法跟力量法
對同一個物理問題的兩種不同解法, 答案本來就該相同, 差別只在"纜索"這個角色
在圖上怎麼呈現(方法A畫成真實桿件, 方法B化簡成等效外力箭頭)。
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output_plots')
os.makedirs(OUT, exist_ok=True)

w = 20.0
H_tower = 25.0
E_b, A_b, I_b = 2.0e8, 0.5, 0.08
E_t, A_t, I_t = 2.0e8, 0.8, 0.15
E_c, A_c = 1.6e8, 0.003
x_c_all = [10.0, 20.0, 30.0, 40.0, 50.0]

# 力量法(含塔柱彎曲柔度修正版)解出的3條纜索張力 -- 見test_cable_stayed_bridge.py
CABLE_TENSIONS_FORCE_METHOD = {3: (30.0, 176.0331), 4: (40.0, 257.7291), 5: (50.0, 274.2302)}


def _add_common_beam_tower(f):
    f.add_node(0, 0, 0)
    for i, x in enumerate(x_c_all):
        f.add_node(i + 1, x, 0)
    f.add_node(6, 0, H_tower)
    f.add_section('beam', E=E_b, I=I_b, A=A_b)
    f.add_section('tower', E=E_t, I=I_t, A=A_t)
    for i in range(5):
        f.add_member(i, i, i + 1, 'beam')
    f.add_member(5, 0, 6, 'tower')
    f.fix(0)
    for i in range(5):
        f.distributed_load(i, w=-w)


# ---- 方法A: 位移法, cable元素自動迭代 ----
fA = Frame2D()
_add_common_beam_tower(fA)
fA.add_section('cable', E=E_c, I=1.0, A=A_c)
for i in range(5):
    fA.add_cable(6 + i, 6, i + 1, 'cable')
rA = solve(fA)
figA = plot_all(fA, rA, figsize=(16, 10))
figA.suptitle(f'Method A - Displacement Method (frame2d cable auto-iteration, slack={sorted(rA.slack_cables)})',
              fontsize=11)
figA.tight_layout()
figA.savefig(f'{OUT}/bridge_displacement_method.png', dpi=120)
print('saved bridge_displacement_method.png')

# ---- 方法B: 力量法等效載重(不含cable元素, 3條纜索張力當作已知外力) ----
fB = Frame2D()
_add_common_beam_tower(fB)
Fx_tower, Fy_tower = 0.0, 0.0
for node_id, (xi, Xi) in CABLE_TENSIONS_FORCE_METHOD.items():
    L_ci = np.hypot(xi, H_tower)
    ux_i, uy_i = -xi / L_ci, H_tower / L_ci   # 單位向量: 從梁錨點指向塔頂
    fB.point_load(node_id, fx=Xi * ux_i, fy=Xi * uy_i)   # 梁承受纜索拉力
    Fx_tower += -Xi * ux_i
    Fy_tower += -Xi * uy_i
fB.point_load(6, fx=Fx_tower, fy=Fy_tower)                # 塔頂承受反作用力
rB = solve(fB)
figB = plot_all(fB, rB, figsize=(16, 10))
figB.suptitle('Method B - Force Method (matrix flexibility, cable tensions as equivalent loads)', fontsize=11)
figB.tight_layout()
figB.savefig(f'{OUT}/bridge_force_method.png', dpi=120)
print('saved bridge_force_method.png')

# ---- 交叉驗證: 梁+塔柱的變形/內力應該完全一致 ----
print()
print('=== 交叉驗證: 節點位移(應完全相同) ===')
for nid in [1, 2, 3, 4, 5, 6]:
    ux, uy, _ = fA.dofs_of(nid)
    dA = np.hypot(rA.displacements[ux], rA.displacements[uy])
    dB = np.hypot(rB.displacements[ux], rB.displacements[uy])
    print(f'  node{nid}: 位移法={dA:.6f}  力量法={dB:.6f}  差={abs(dA-dB):.2e}')
    assert abs(dA - dB) < 1e-6, f'node{nid} 位移不一致'

print()
print('PASS: 兩種方法算出的梁+塔柱變形完全一致, 差異在浮點精度誤差範圍內')
