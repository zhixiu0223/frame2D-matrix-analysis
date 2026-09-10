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

print("PASS: frame2d.newton所有案例通過")
