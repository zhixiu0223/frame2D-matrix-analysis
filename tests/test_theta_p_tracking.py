"""
驗證案例: HingeState.theta_p(累積塑性轉角)在run_pushover()裡的追蹤

theta_p是用hinge.beam_internal_rotation()(獨立sympy驗證過, 見
test_hinge_condensation.py)算出的"樑內部端點轉角"跟"外部節點轉角"的差,
不是外部節點轉角本身(見hinge.py跟pushover.py裡對這個區別的說明: 鉸
剛好在完全固定支承上時, 外部節點轉角永遠是0, 不能拿來當塑性轉角)。

案例A: 未降伏的端點, theta_p應該永遠是0(不管推多遠)。
案例B: 已降伏的端點, 側推距離越遠, theta_p應該越大(單調遞增) -- 這是
       "圈圈越推越大"這個視覺化需求實際依賴的性質。
案例C: 對稱門型鋼架(test_pushover_mechanism.py的案例)裡, 兩根柱是
       對稱的, 最終theta_p應該幾乎相等(同一端)。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0
Mp_base = 100.0


def cantilever_column():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec',
                 Mp_i=Mp_base, Mp_j=1e30, R_post_yield_i=1e4, R_post_yield_j=1e4)
    f.fix(0)
    return f


# ---- 案例A: 未降伏端點theta_p永遠是0 ----
print("=== 案例A: 未降伏端點theta_p永遠是0 ===")
f = cantilever_column()
hinge_states = initial_hinge_states(f)
u_hist, F_hist, event_log, hs, mechanism = run_pushover(
    f, hinge_states, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
    target_total=0.05, d_nominal=0.001, base_reaction_dofs=[f.dofs_of(0)[0]],
)
assert hs[0].yielded == [True, False]
assert hs[0].theta_p[1] == 0.0, f"端1(柱頂, 永遠不降伏)theta_p應該是0, 實際={hs[0].theta_p[1]}"
assert hs[0].theta_p[0] > 0.0, f"端0(柱底, 已降伏)theta_p應該大於0, 實際={hs[0].theta_p[0]}"
print(f"PASS: 未降伏端theta_p=0, 已降伏端theta_p={hs[0].theta_p[0]:.6f} > 0\n")


# ---- 案例B: 推越遠, theta_p越大(單調遞增) ----
print("=== 案例B: 側推距離越遠, theta_p單調遞增 ===")
theta_p_values = []
for target in (0.035, 0.04, 0.045, 0.05):
    f_i = cantilever_column()
    hs_i = initial_hinge_states(f_i)
    _, _, _, hs_result, _ = run_pushover(
        f_i, hs_i, prescribed_dofs=[f_i.dofs_of(1)[0]], direction=[1.0],
        target_total=target, d_nominal=0.001, base_reaction_dofs=[f_i.dofs_of(0)[0]],
    )
    theta_p_values.append(hs_result[0].theta_p[0])

for i in range(1, len(theta_p_values)):
    assert theta_p_values[i] > theta_p_values[i - 1], (
        f"theta_p應該隨側推距離單調遞增, 實際={theta_p_values}"
    )
print(f"PASS: theta_p隨側推距離單調遞增 {[f'{v:.6f}' for v in theta_p_values]}\n")


# ---- 案例C: 對稱門型鋼架, 兩根柱最終theta_p應該幾乎相等 ----
print("=== 案例C: 對稱結構下兩根柱的theta_p應該幾乎相等 ===")
Mp, R_py, H, Lspan = 100.0, 1.0, 4.0, 6.0
f2 = Frame2D()
f2.add_node(0, 0, 0); f2.add_node(1, Lspan, 0)
f2.add_node(2, 0, H); f2.add_node(3, Lspan, H)
f2.add_section('col', E=E, I=I, A=A)
f2.add_section('beam', E=E, I=I, A=A)
f2.add_member(0, node_i=0, node_j=2, section='col', Mp_i=Mp, Mp_j=Mp, R_post_yield_i=R_py, R_post_yield_j=R_py)
f2.add_member(1, node_i=1, node_j=3, section='col', Mp_i=Mp, Mp_j=Mp, R_post_yield_i=R_py, R_post_yield_j=R_py)
f2.add_member(2, node_i=2, node_j=3, section='beam')
f2.fix(0); f2.fix(1)
hinge_states2 = initial_hinge_states(f2)
_, _, _, hs2, mech2 = run_pushover(
    f2, hinge_states2, prescribed_dofs=[f2.dofs_of(2)[0], f2.dofs_of(3)[0]], direction=[1.0, 1.0],
    target_total=0.045, d_nominal=0.0005, base_reaction_dofs=[f2.dofs_of(0)[0], f2.dofs_of(1)[0]],
    mechanism_ratio_limit=1e-6,
)
for end_idx in (0, 1):
    tp0, tp1 = hs2[0].theta_p[end_idx], hs2[1].theta_p[end_idx]
    assert np.isclose(tp0, tp1, rtol=1e-6), (
        f"對稱結構下兩柱端{end_idx}的theta_p應該幾乎相等, 實際={tp0} vs {tp1}"
    )
print(f"PASS: 兩柱theta_p對稱一致 (端0: {hs2[0].theta_p[0]:.6f}≈{hs2[1].theta_p[0]:.6f}, "
      f"端1: {hs2[0].theta_p[1]:.6f}≈{hs2[1].theta_p[1]:.6f})\n")

print("PASS: theta_p追蹤所有案例通過")
