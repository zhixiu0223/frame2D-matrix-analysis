"""
驗證案例: frame2d.corotational + frame2d.newton -- 真正的大轉角
co-rotational幾何 + Newton-Raphson平衡疊代求解器。

跟run_pushover()/run_pushover_converged()是完全獨立的第三種求解器,
見frame2d/newton.py開頭docstring的定位說明。

這裡驗證的重點, 依開發過程實際發生的問題排序(不是隨便湊案例):
1. 核心正確性(co-rotational公式本身): 桿件被剛體轉動再多度(不是
   小角度), 只要沒有真的變形, 內力一定是0——這是舊的geometry_update
   做不到、被實際案例找出來的限制, 這裡直接驗證修正了。
2. 跟解析解一致: 彈性懸臂樑, 位移控制/力控制都要對得上F=3EIu/L^3
   (在小位移極限下)。
3. 大位移時的二階效應: 誤差應該隨(u/L)^2縮放(真實的幾何非線性
   效應, 不是bug)。
4. 材料降伏(塑鉸)的路徑相依性: 降伏瞬間彎矩必須連續(不能像
   開發過程中發現的那樣暴跌到接近0), 且收斂後應該符合
   Mp+R*theta_p這個雙折線硬化模型的物理關係。
5. 不收斂時誠實回報, 不會給不可信的答案。
6. release端明確拒絕(還沒支援, 不是給錯誤答案)。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.hinge import HingeState
from frame2d.corotational import corotational_global_force
from frame2d.newton import run_pushover_newton

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


def cantilever_no_hinge():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec')
    f.fix(0)
    return f


# ---- 案例1: 剛體轉動(不管多大)不產生任何內力 ----
print("=== 案例1: 大角度剛體轉動內力應該精確為0 ===")
for angle_deg in [10, 60, 170]:
    theta = np.radians(angle_deg)
    xj_new, yj_new = L * np.cos(theta), L * np.sin(theta)
    # 桿件原始沿+x方向(0,0)-(L,0), 剛體轉動theta角度後j端新位置是
    # (xj_new,yj_new), 對應位移uj=xj_new-L, vj=yj_new-0
    f_elem = corotational_global_force(0, 0, 0, 0, theta, L, 0, xj_new - L, yj_new, theta, E, A, I)
    assert np.max(np.abs(f_elem)) < 1e-8, f"{angle_deg}度剛體轉動應該內力=0, 實際最大值={np.max(np.abs(f_elem))}"
print("PASS: 10/60/170度剛體轉動, 內力都精確為0(co-rotational核心正確性)\n")


# ---- 案例2: 彈性懸臂, 位移控制/力控制都對得上解析解 ----
print("=== 案例2: 彈性懸臂樑對得上F=3EIu/L^3解析解 ===")
f1 = cantilever_no_hinge()
u1, F1, ev1, hs1, conv1 = run_pushover_newton(
    f1, {}, prescribed_dofs=[f1.dofs_of(1)[0]], direction=[1.0],
    target_total=0.0005, d_nominal=0.00005, base_reaction_dofs=[f1.dofs_of(0)[0]],
    control_mode='displacement', tol=1e-9,
)
F_theory = 3 * E * I * 0.0005 / L**3
assert abs(F1[-1] - F_theory) / F_theory < 1e-4, "小位移下應該精確對得上解析解"
assert conv1 is True
print(f"PASS(位移控制): F={F1[-1]:.6f}, 理論={F_theory:.6f}\n")

f2 = cantilever_no_hinge()
F_target = 5.0
u2, F2, ev2, hs2, conv2 = run_pushover_newton(
    f2, {}, prescribed_dofs=[f2.dofs_of(1)[0]], direction=[1.0],
    target_total=F_target, d_nominal=0.5, base_reaction_dofs=[f2.dofs_of(0)[0]],
    control_mode='force', tol=1e-9,
)
u_theory = F_target * L**3 / (3 * E * I)
assert abs(u2[-1] - u_theory) / u_theory < 1e-4, "小力下應該精確對得上解析解"
print(f"PASS(力控制): u={u2[-1]:.6f}, 理論={u_theory:.6f}\n")


# ---- 案例3: 大位移時誤差隨(u/L)^2縮放(真實二階效應, 不是bug) ----
print("=== 案例3: 大位移時的誤差應該隨(u/L)^2縮放 ===")
errs = []
for u_target in [0.05, 0.005]:
    f3 = cantilever_no_hinge()
    u3, F3, _, _, _ = run_pushover_newton(
        f3, {}, prescribed_dofs=[f3.dofs_of(1)[0]], direction=[1.0],
        target_total=u_target, d_nominal=u_target / 10, base_reaction_dofs=[f3.dofs_of(0)[0]],
        control_mode='displacement', tol=1e-9,
    )
    F_theory = 3 * E * I * u_target / L**3
    errs.append(abs(F3[-1] - F_theory) / F_theory)
ratio = errs[0] / errs[1]
assert 80 < ratio < 120, f"位移縮小10倍, 誤差應該縮小約100倍((u/L)^2效應), 實際比值={ratio}"
print(f"PASS: u=0.05時誤差={errs[0]:.3e}, u=0.005時誤差={errs[1]:.3e}, 比值={ratio:.1f}(約100)\n")


# ---- 案例4: 塑鉸降伏瞬間彎矩必須連續(不能暴跌) ----
print("=== 案例4: 塑鉸降伏瞬間彎矩連續, 不會像修正前那樣暴跌到接近0 ===")
f4 = Frame2D()
f4.add_node(0, 0, 0)
f4.add_node(1, 0, L)
f4.add_section('sec', E=E, I=I, A=A)
f4.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
              R_post_yield_i=100.0, R_post_yield_j=100.0)
f4.fix(0)
hs4 = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=100.0, R_post_yield_2=100.0)}
u4, F4, ev4, hsf4, conv4, snaps4 = run_pushover_newton(
    f4, hs4, prescribed_dofs=[f4.dofs_of(1)[0]], direction=[1.0],
    target_total=0.08, d_nominal=0.002, base_reaction_dofs=[f4.dofs_of(0)[0]],
    control_mode='displacement', tol=1e-8, include_snapshots=True,
)
assert len(ev4) == 1, f"應該有1個降伏事件, 實際={len(ev4)}"
yield_step = None
for i, snap in enumerate(snaps4):
    if snap['hinge_states'][0]['yielded'][0]:
        yield_step = i
        break
M1_before = snaps4[yield_step - 1]['member_forces'][0][2]
M1_after = snaps4[yield_step]['member_forces'][0][2]
assert M1_before > 80, f"降伏前一步M1應該接近Mp=100, 實際={M1_before}"
assert M1_after > 80, (
    f"降伏後那一步M1應該保持接近Mp、連續, 不應該暴跌"
    f"(修正前的bug是這裡會掉到<1)——實際M1_before={M1_before}, M1_after={M1_after}"
)
print(f"PASS: 降伏前M1={M1_before:.2f}, 降伏後M1={M1_after:.2f}(連續, 沒有暴跌)\n")


# ---- 案例5: 收斂後的最終彎矩符合Mp+R*theta_p(雙折線硬化模型) ----
print("=== 案例5: 最終彎矩應該符合Mp+R*theta_p ===")
M1_final = snaps4[-1]['member_forces'][0][2]
theta_p_final = hsf4[0].theta_p[0]
M_expected = 100.0 + 100.0 * theta_p_final
rel_err = abs(M1_final - M_expected) / M_expected
assert rel_err < 0.15, (
    f"最終彎矩應該大致符合Mp+R*theta_p(容許step-size相關的imprecision,"
    f"見newton.py已知限制說明——這裡只是要抓修正前那種~99%的災難性錯誤,"
    f"不是要驗證到很高精度), 實際相對誤差={rel_err}"
)
print(f"PASS: 最終M1={M1_final:.3f}, Mp+R*theta_p={M_expected:.3f}, 相對誤差={rel_err:.3f}(<15%)\n")


# ---- 案例6: 不收斂時誠實回報, 不給不可信的答案 ----
print("=== 案例6: 不收斂時應該誠實回報converged=False ===")
f6 = cantilever_no_hinge()
u6, F6, ev6, hs6, conv6 = run_pushover_newton(
    f6, {}, prescribed_dofs=[f6.dofs_of(1)[0]], direction=[1.0],
    target_total=3.9, d_nominal=3.9, base_reaction_dofs=[f6.dofs_of(0)[0]],
    control_mode='displacement', max_iter=3, tol=1e-12,
)
assert conv6 is False, "極端大步長+疊代次數壓低, 應該真的偵測到不收斂"
print(f"PASS: 正確回報converged=False, 在u={u6[-1]}處停止\n")


# ---- 案例7: release端明確拒絕, 不是給錯誤答案 ----
print("=== 案例7: release端應該明確拒絕(還沒支援) ===")
f7 = Frame2D()
f7.add_node(0, 0, 0)
f7.add_node(1, 0, L)
f7.add_section('sec', E=E, I=I, A=A)
f7.add_member(0, node_i=0, node_j=1, section='sec', release_j=True)
f7.fix(0)
try:
    run_pushover_newton(f7, {}, prescribed_dofs=[f7.dofs_of(1)[0]], direction=[1.0],
                         target_total=0.01, d_nominal=0.001, base_reaction_dofs=[f7.dofs_of(0)[0]])
    assert False, "應該要raise ValueError"
except ValueError as e:
    assert 'release' in str(e)
print("PASS: release端正確raise ValueError, 不是給錯誤答案\n")

# ---- 案例8: 重力/桿件內部載重支援, 對得上已驗證過的apply_gravity() ----
print("=== 案例8: 重力(分佈載重)支援, 精確對得上apply_gravity()跟wL^2/12公式 ===")
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import apply_gravity
from frame2d.newton import _gravity_fixed_end_forces


def fixed_fixed_beam_with_udl():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.distributed_load(0, w=10.0, direction='global_y')
    f.fix(0)
    f.fix(1)
    return f


f8a = fixed_fixed_beam_with_udl()
hs8a = initial_hinge_states(f8a)
init_forces8, _ = apply_gravity(f8a, hs8a)
M1_ref = init_forces8[0][2]
M2_ref = init_forces8[0][5]
assert abs(M1_ref - 10.0 * L**2 / 12) < 1e-9, "跟已驗證的apply_gravity()自己先要對得上wL^2/12"

# 直接測_gravity_fixed_end_forces()本身(不透過完整的pushover迴圈,
# 因為這個案例兩端都完全固定、沒有任何自由度可以推, 直接測底層函式
# 比硬湊一個可以推的模型更乾淨): 固定端彎矩(member_fem_local的索引
# 2,5)取負號後應該直接對得上apply_gravity()的end_forces_local慣例
# (見dofmanager.py: end_forces_local = k_local@u_local - f_FE, 變形
# 量在u=0時是0, 所以應該精確等於-f_FE)。
f8b = fixed_fixed_beam_with_udl()
f_ext_gravity8, member_fem_local8 = _gravity_fixed_end_forces(f8b)
M1_newton = -member_fem_local8[0][2]
M2_newton = -member_fem_local8[0][5]
assert abs(M1_newton - M1_ref) < 1e-9 and abs(M2_newton - M2_ref) < 1e-9, (
    f"_gravity_fixed_end_forces()算出的固定端彎矩(取負號後)應該精確"
    f"對得上apply_gravity(), 實際: newton=({M1_newton},{M2_newton}), "
    f"apply_gravity=({M1_ref},{M2_ref})"
)
print(f"PASS: _gravity_fixed_end_forces()精確對得上apply_gravity()跟wL^2/12公式"
      f"(M1={M1_newton:.6f}, M2={M2_newton:.6f})\n")

# ---- 案例9: 重力+側推力同時作用, 完整跑過一次不crash且合理收斂 ----
print("=== 案例9: 重力+側推同時作用的完整案例 ===")
f9 = Frame2D()
f9.add_node(0, 0, 0)
f9.add_node(1, 0, L)
f9.add_section('sec', E=E, I=I, A=A)
f9.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
              R_post_yield_i=100.0, R_post_yield_j=100.0)
f9.distributed_load(0, w=10.0, direction='global_y')
f9.fix(0)
hs9 = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=100.0, R_post_yield_2=100.0)}
u9, F9, ev9, hsf9, conv9 = run_pushover_newton(
    f9, hs9, prescribed_dofs=[f9.dofs_of(1)[0]], direction=[1.0],
    target_total=0.05, d_nominal=0.002, base_reaction_dofs=[f9.dofs_of(0)[0]],
    control_mode='displacement', tol=1e-8,
)
assert conv9 is True, "重力+側推同時作用應該能正常收斂"
print(f"PASS: converged={conv9}, 降伏事件數={len(ev9)}, 最終F={F9[-1]:.4f}\n")


# ---- 案例10: 快照要有完整的Fx/Fy/M(不是只有M1,M2), 否則精確N/V/M/
# 變形圖(member_internal_forces())重建V(x)時會把剪力當成0, 柱子的
# 彎矩圖會被錯誤畫成上下端數值相同的"矩形"——這是使用者用截圖抓出來
# 的真實bug, 見對話紀錄。 ----
print("=== 案例10: 快照的端力向量要有完整Fy(不能只有M1,M2=剪力=0的假象) ===")
from frame2d.newton import _member_end_forces_local

f10 = Frame2D()
f10.add_node(0, 0, 0)
f10.add_node(1, 0, L)
f10.add_section('sec', E=E, I=I, A=A)
f10.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
               R_post_yield_i=1.0, R_post_yield_j=1.0)
f10.fix(0)
hs10 = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u10, F10, ev10, hsf10, conv10, snaps10 = run_pushover_newton(
    f10, hs10, prescribed_dofs=[f10.dofs_of(1)[0]], direction=[1.0],
    target_total=0.01, d_nominal=0.005, base_reaction_dofs=[f10.dofs_of(0)[0]],
    control_mode='displacement', tol=1e-9, include_snapshots=True,
)
Fy1_snap = snaps10[-1]['member_forces'][0][1]
assert abs(Fy1_snap) > 1e-6, (
    f"快照裡的Fy1不應該是0(修正前的bug就是永遠存0, 讓member_internal_"
    f"forces()重建V(x)/M(x)時剪力消失, 柱子彎矩圖變成上下端一樣的"
    f"矩形)——實際Fy1={Fy1_snap}"
)
M1_snap, M2_snap = snaps10[-1]['member_forces'][0][2], snaps10[-1]['member_forces'][0][5]
assert abs(Fy1_snap - (M1_snap + M2_snap) / L) < 1e-6, (
    "快照裡的Fy1應該精確符合(M1+M2)/L這個樑元素平衡關係"
)
print(f"PASS: 快照正確存有Fy1={Fy1_snap:.4f}(不是0), 且符合(M1+M2)/L平衡關係\n")

# ---- 案例11: 重力預載階段要真的解平衡方程式, 不能天真假設u=0就是
# "重力施加前"的狀態——這是使用者用實際截圖(柱子彎矩顯示成0)抓出來
# 的真實bug, 見對話紀錄。用一個不對稱跨度的2層樓構架(跟apply_gravity()
# 對照), 驗證柱子在"第一步"(側推還沒開始, 只有重力)應該有非零彎矩
# (透過剛接節點, 樑的彎曲會傳一部分進柱子), 不是0。 ----
print("=== 案例11: 重力預載要真的疊代求解平衡, 柱子彎矩不應該是0 ===")
from frame2d.dofmanager import initial_hinge_states as _initial_hinge_states


def asym_two_story_frame():
    f = Frame2D()
    coords = {0: (0, 0), 1: (5, 0), 2: (12, 0), 3: (0, 4), 4: (5, 4), 5: (12, 4),
              6: (0, 7.5), 7: (5, 7.5), 8: (12, 7.5)}
    for nid, (x, y) in coords.items():
        f.add_node(nid, x, y)
    f.add_section('sec', E=208e9, I=23500e-8, A=83.37e-4)
    Rmap = {0: 1466.4e3, 1: 1466.4e3, 2: 1466.4e3, 3: 1173.12e3, 4: 837.943e3,
            5: 1675.89e3, 6: 1675.89e3, 7: 1675.89e3, 8: 1173.12e3, 9: 837.943e3}
    conns = [(0, 0, 3), (1, 1, 4), (2, 5, 2), (3, 3, 4), (4, 4, 5),
             (5, 6, 3), (6, 7, 4), (7, 8, 5), (8, 6, 7), (9, 7, 8)]
    for mid, ni, nj in conns:
        f.add_member(mid, node_i=ni, node_j=nj, section='sec', Mp_i=322.5e3, Mp_j=322.5e3,
                     R_post_yield_i=Rmap[mid], R_post_yield_j=Rmap[mid])
    for mid in range(10):
        f.distributed_load(mid, w=15e3, direction='global_y')
    f.fix(0)
    f.fix(1)
    f.fix(2)
    return f


f11a = asym_two_story_frame()
hs11a = _initial_hinge_states(f11a)
init_forces11, _ = apply_gravity(f11a, hs11a)
M1_ref11 = init_forces11[0][2]
assert abs(M1_ref11) > 100, "跟已驗證的apply_gravity()自己先要確認柱子底端彎矩明顯不是0"

f11b = asym_two_story_frame()
hs11b = {mid: HingeState(Mp1=322.5e3, Mp2=322.5e3, R_post_yield_1=1e10, R_post_yield_2=1e10)
         for mid in range(10)}
u11, F11, ev11, hsf11, conv11, u_full11, snaps11 = run_pushover_newton(
    f11b, hs11b, prescribed_dofs=[f11b.dofs_of(7)[0]], direction=[1.0],
    target_total=1e-9, d_nominal=1e-9, base_reaction_dofs=[f11b.dofs_of(0)[0], f11b.dofs_of(1)[0], f11b.dofs_of(2)[0]],
    control_mode='displacement', tol=1e-9, include_final_displacement=True, include_snapshots=True, max_iter=50,
)
M1_newton11 = snaps11[0]['member_forces'][0][2]
assert abs(M1_newton11) > 100, (
    f"修正前的bug: 柱子在重力預載階段會顯示彎矩=0(天真地假設u=0是"
    f"重力施加前的狀態, 沒有真的疊代求解重力平衡)——修正後應該明顯"
    f"不是0, 實際={M1_newton11}"
)
rel_err11 = abs(M1_newton11 - M1_ref11) / abs(M1_ref11)
assert rel_err11 < 0.01, (
    f"newton版本的重力預載柱子底端彎矩應該精確對得上apply_gravity(), "
    f"實際: newton={M1_newton11}, apply_gravity={M1_ref11}, 相對誤差={rel_err11}"
)
print(f"PASS: 重力預載階段柱子底端彎矩={M1_newton11:.2f}(不是0), "
      f"對得上apply_gravity()的{M1_ref11:.2f}(相對誤差{rel_err11:.4f})\n")


# ---- 案例12: frame.point_loads(直接節點力)要真的被套用, 不能靜默
# 忽略——這是補上P-Delta驗證時順便發現的另一個真實缺口(newton.py
# 原本完全沒有處理point_loads, 只處理了distributed_loads/
# member_point_loads)。 ----
print("=== 案例12: frame.point_loads要真的套用, 不能被靜默忽略 ===")
f12 = Frame2D()
f12.add_node(0, 0, 0)
f12.add_node(1, 0, L)
f12.add_section('sec', E=E, I=I, A=A)
f12.add_member(0, node_i=0, node_j=1, section='sec')
f12.fix(0)
f12.point_load(1, fx=100.0, fy=0, m=0)
u12, F12, ev12, hsf12, conv12 = run_pushover_newton(
    f12, {}, prescribed_dofs=[f12.dofs_of(1)[0]], direction=[1.0],
    target_total=1e-9, d_nominal=1e-9, base_reaction_dofs=[f12.dofs_of(0)[0]],
    control_mode='displacement', tol=1e-9,
)
from frame2d.newton import _gravity_fixed_end_forces
f_ext_g, _ = _gravity_fixed_end_forces(f12)
assert abs(f_ext_g[f12.dofs_of(1)[0]] - 100.0) < 1e-9, (
    f"point_loads應該直接加進f_ext_gravity, 實際={f_ext_g[f12.dofs_of(1)[0]]}"
)
print(f"PASS: point_loads正確套用進f_ext_gravity(={f_ext_g[f12.dofs_of(1)[0]]:.2f})\n")


# ---- 案例13: 局部P-Delta(軸力對桿件自身彎曲勁度的修正, 用
# elements.local_geometric_stiffness()的theta子矩陣)——接近挫屈臨界力
# 的受壓柱, 開這個選項應該讓側向勁度明顯變軟(位移明顯變大)。這跟
# co-rotational大轉角是完全獨立的兩件事(見對話紀錄的完整討論): 大
# 轉角處理"桿件整體轉了多少度", 局部P-Delta處理"軸力怎麼影響桿件
# 自己的彎曲勁度"。 ----
print("=== 案例13: 局部P-Delta(軸力對桿件自身彎曲勁度)接近挫屈臨界力時應該明顯軟化 ===")
P_cr = np.pi**2 * E * I / (4 * L**2)   # 懸臂樑歐拉挫屈臨界力
P_axial = 0.81 * P_cr

def cantilever_with_axial_compression():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec')
    f.fix(0)
    f.point_load(1, fx=0, fy=-P_axial, m=0)
    return f


f13a = cantilever_with_axial_compression()
u13a, F13a, ev13a, hsf13a, conv13a = run_pushover_newton(
    f13a, {}, prescribed_dofs=[f13a.dofs_of(1)[0]], direction=[1.0],
    target_total=100.0, d_nominal=10.0, base_reaction_dofs=[f13a.dofs_of(0)[0]],
    control_mode='force', tol=1e-9, use_pdelta=False,
)
f13b = cantilever_with_axial_compression()
u13b, F13b, ev13b, hsf13b, conv13b = run_pushover_newton(
    f13b, {}, prescribed_dofs=[f13b.dofs_of(1)[0]], direction=[1.0],
    target_total=100.0, d_nominal=10.0, base_reaction_dofs=[f13b.dofs_of(0)[0]],
    control_mode='force', tol=1e-9, use_pdelta=True,
)
assert u13b[-1] > u13a[-1] * 1.2, (
    f"接近挫屈臨界力時, 開局部P-Delta應該讓同樣的側推力產生明顯更大的"
    f"位移(勁度變軟), 實際: 沒開={u13a[-1]:.4f}, 開了={u13b[-1]:.4f}"
)
print(f"PASS: 沒開P-Delta u={u13a[-1]:.4f}, 開了u={u13b[-1]:.4f}"
      f"(明顯變大, 符合接近挫屈臨界力時應有的軟化效應)\n")

print("PASS: frame2d.newton所有案例通過")
