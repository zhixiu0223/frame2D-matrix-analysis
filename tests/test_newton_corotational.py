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

print("PASS: frame2d.newton所有案例通過")
