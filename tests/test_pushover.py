"""
驗證案例: frame2d/pushover.py -- 遞增側推(位移控制 + event-to-event)

用最小、最好獨立驗算的案例: 一根固定基底的垂直懸臂柱, 柱底(靠近固定
支承的那一端)設塑鉸容量Mp, 柱頂設塑鉸容量極大(永遠不會降伏), 在柱頂
做水平位移控制側推。這個案例的三個關鍵數字都能獨立算出來跟
run_pushover()的輸出比對, 不用另外架OpenSeesPy模型:

  1. 降伏前的側向勁度: 標準教科書公式3EI/L^3(懸臂梁, 頂端自由轉動、
     只受側向力)
  2. 降伏發生的位置: 底部彎矩M=F*L(靜力學), 應該恰好在F=Mp/L時觸發
  3. 降伏後的側向勁度: 用已經在test_hinge_condensation.py驗證過的
     hinge_bending_stiffness()公式, 把柱頂theta(不受任何彎矩限制的
     自由轉角)靜力凝聚掉, 得到的縮減1自由度勁度

案例A: 上述三個數字, run_pushover()的輸出跟獨立算出來的理論值逐項吻合。
案例B: mechanism_reached應該是False(單一鉸降伏後仍有正的殘餘勁度,
       不是機構), hinge_states[0].yielded應該是[True, False](只有
       柱底那一端降伏, 柱頂那一端Mp設極大, 永遠不會降伏)。
案例C: apply_gravity() -- 先用重力載重跑一次力控制求解, 取得的cum_forces
       應該跟直接呼叫solve_with_hinges()算出來的end_forces_local逐項
       一致(這是_solve_once_dofmanager既有機制的直接重用, 不是另外
       發明一套邏輯)。
案例D: max_steps安全閥 -- d_nominal設成比target_total大很多倍時
       (只需要一步), 不應該誤觸發RuntimeError。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.hinge import HingeState, hinge_bending_stiffness
from frame2d.dofmanager import initial_hinge_states, solve_with_hinges
from frame2d.pushover import run_pushover, apply_gravity

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0
Mp_base = 100.0


def cantilever_column():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)   # 垂直柱: 側推方向(水平)跟軸向(垂直)分開,
                          # 是真實pushover的設定方式, 不是為了測試方便
                          # 特意擺成水平
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec',
                 Mp_i=Mp_base, Mp_j=1e30, R_post_yield_i=1e4, R_post_yield_j=1e4)
    f.fix(0)
    return f


# ---- 案例A: 三個關鍵數字逐項獨立驗算 ----
print("=== 案例A: 彈性勁度/降伏點/降伏後勁度, 三個數字獨立驗算 ===")
f = cantilever_column()
hinge_states = initial_hinge_states(f)
control_dof = f.dofs_of(1)[0]   # node1的ux(水平位移控制)
base_dof = f.dofs_of(0)[0]      # node0的ux反力(底部水平反力=側推的底剪力)

u_hist, F_hist, event_log, hs, mechanism = run_pushover(
    f, hinge_states, prescribed_dofs=[control_dof], direction=[1.0],
    target_total=0.05, d_nominal=0.001, base_reaction_dofs=[base_dof],
)

K_elastic_theory = 3 * E * I / L**3
K_elastic_numeric = (F_hist[1] - F_hist[0]) / (u_hist[1] - u_hist[0])
assert np.isclose(K_elastic_numeric, K_elastic_theory, rtol=1e-6), (
    f"降伏前側向勁度應該是教科書3EI/L^3={K_elastic_theory}, 實際={K_elastic_numeric}"
)
print(f"PASS: 降伏前側向勁度 {K_elastic_numeric:.4f} 跟教科書3EI/L^3={K_elastic_theory:.4f}吻合")

assert len(event_log) == 1, f"這個案例應該恰好只有1個降伏事件, 實際={len(event_log)}"
F_at_yield = event_log[0]['F']
F_yield_theory = Mp_base / L   # 靜力學: 底部彎矩M=F*L, 降伏時M=Mp
assert np.isclose(F_at_yield, F_yield_theory, rtol=1e-9), (
    f"降伏時的底剪力應該是Mp/L={F_yield_theory}, 實際={F_at_yield}"
)
print(f"PASS: 降伏發生在F={F_at_yield:.4f}, 跟靜力學Mp/L={F_yield_theory:.4f}吻合")

hs_post = HingeState(Mp1=Mp_base, Mp2=1e30, R_post_yield_1=1e4, R_post_yield_2=1e4)
hs_post.yielded = [True, False]
Kb = hinge_bending_stiffness(E, I, L, hs_post)   # (v1,theta1,v2,theta2)
K_free = Kb[np.ix_([2, 3], [2, 3])]              # node0固定, 只留v2,theta2
K_post_yield_theory = K_free[0, 0] - K_free[0, 1] * K_free[1, 0] / K_free[1, 1]  # 凝聚掉theta2
K_post_yield_numeric = (F_hist[-1] - F_hist[-2]) / (u_hist[-1] - u_hist[-2])
assert np.isclose(K_post_yield_numeric, K_post_yield_theory, rtol=1e-4), (
    f"降伏後側向勁度應該是{K_post_yield_theory:.4f}(用已驗證過的hinge_bending_stiffness"
    f"凝聚theta2算出來的), 實際={K_post_yield_numeric:.4f}"
)
print(f"PASS: 降伏後側向勁度 {K_post_yield_numeric:.4f} 跟獨立凝聚算出的"
      f"{K_post_yield_theory:.4f}吻合(容許誤差來自步長離散化, rtol=1e-4)\n")


# ---- 案例B: 機構偵測 + 降伏狀態正確性 ----
print("=== 案例B: 單一鉸降伏不是機構, 只有柱底那一端降伏 ===")
assert mechanism is False, "單一鉸降伏後仍有正的殘餘勁度, 不應該判定為機構"
assert hs[0].yielded == [True, False], (
    f"應該只有端0(柱底)降伏, 端1(柱頂, Mp設1e30永遠不會降伏)應該仍是False, "
    f"實際={hs[0].yielded}"
)
print("PASS: mechanism_reached=False, yielded=[True, False]\n")


# ---- 案例C: apply_gravity()直接重用solve_with_hinges(), 不是另一套邏輯 ----
print("=== 案例C: apply_gravity()跟solve_with_hinges()一致 ===")
f2 = cantilever_column()
f2.point_load(1, fx=-5.0)   # 隨便給一個側向力, 純粹測試機制本身
hinge_states2 = initial_hinge_states(f2)
cum_forces, grav_result = apply_gravity(f2, hinge_states2)

direct_result = solve_with_hinges(f2, hinge_states2)
for mid in f2.members:
    assert np.allclose(cum_forces[mid], direct_result.member_results[mid].end_forces_local), (
        f"apply_gravity()回傳的cum_forces應該跟直接呼叫solve_with_hinges()的"
        f"end_forces_local逐項一致(member {mid})"
    )
print("PASS: apply_gravity()的cum_forces跟solve_with_hinges()逐項一致\n")


# ---- 案例D: max_steps安全閥不應該誤觸發 ----
print("=== 案例D: d_nominal >= target_total時(只需要1步)不應該誤觸發max_steps ===")
f3 = cantilever_column()
hinge_states3 = initial_hinge_states(f3)
u3, F3, ev3, hs3, mech3 = run_pushover(
    f3, hinge_states3, prescribed_dofs=[f3.dofs_of(1)[0]], direction=[1.0],
    target_total=0.001, d_nominal=1.0, base_reaction_dofs=[f3.dofs_of(0)[0]],
)
assert len(u3) >= 2, "應該正常跑完, 不應該raise"
print("PASS: 大步長(一步走完)沒有誤觸發max_steps\n")

print("PASS: frame2d/pushover.py 所有案例通過")
