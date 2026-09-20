"""
驗證支座沉陷(Support的非0數值)、溫度效應(thermal_load)在Newton/
corotational_oneshot這兩個求解器裡也正確生效——這兩個之前完全沒有
讀取支承的指定位移值(只當作"是否固定", 不管數值), 等同靜默忽略沉陷,
這次補上。

驗證策略: 拿已經驗證過的線性解(frame2d.dofmanager.solve, 見
test_thermal_load.py/既有測試)當比對基準, 確認Newton疊代到的反力
在合理誤差內(co-rotational跟線性理論本來就會有微小的幾何非線性
差異, 不是bug, 見對話紀錄裡其他跨求解器比對案例)。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.hinge import HingeState
from frame2d.newton import (run_pushover_newton, run_pushover_corotational_oneshot,
                             _assemble_global, _member_theta_def, _member_moments)

E, I, A, L = 200e9, 8e-5, 1e-2, 6.0


def three_node_settlement_model():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_node(2, 2 * L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.add_member(1, node_i=1, node_j=2, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.support(0, ux=0.0, uy=0.0, rot=0.0)
    f.support(1, ux=0.0, uy=-0.02, rot=None)   # 支承1沉陷2cm
    f.support(2, ux=0.0, uy=0.0, rot=None)
    return f


def _newton_reaction_moment(f, hs, u_full):
    """用正確的member_ref重新組裝, 取node0的彎矩反力(跟run_pushover_
    newton()內部收斂後更新member_ref的流程一致, 不能用member_ref全0
    去重算, 那樣member_ref會跟u_full已經有的變形量不一致, 得到錯誤的
    反力數字。"""
    member_ref = {mid: (0.0, 0.0, 0.0, 0.0) for mid in f.members}
    n_dof = 3 * len(f.nodes)
    for mid in f.members:
        th1d, th2d = _member_theta_def(f, mid, u_full)
        M1_def, M2_def = _member_moments(f, mid, u_full, hs, member_ref, use_pdelta=False)
        member_ref[mid] = (M1_def, M2_def, th1d, th2d)
    f_int, _ = _assemble_global(f, hs, u_full, n_dof, member_ref, False)
    return f_int[f.dofs_of(0)[2]]


print("=== 案例1: 支座沉陷在Newton裡正確生效(位移精確等於指定值,")
print("    反力在合理誤差內對得上線性解) ===")
f_lin = three_node_settlement_model()
r_lin = solve(f_lin)
M0_theory = r_lin.reactions[f_lin.dofs_of(0)[2]]

f_nt = three_node_settlement_model()
hs = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0),
      1: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u, F_, ev, hsf, conv, u_full = run_pushover_newton(
    f_nt, hs, prescribed_dofs=[f_nt.dofs_of(1)[1]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f_nt.dofs_of(0)[0]], tol=1e-6, include_final_displacement=True,
)
uy1 = u_full[f_nt.dofs_of(1)[1]]
assert abs(uy1 - (-0.02)) < 1e-12, f"支承1的位移應該精確等於指定的沉陷值-0.02, 實際={uy1}"
M0_newton = _newton_reaction_moment(f_nt, hs, u_full)
rel_err = abs(M0_newton - M0_theory) / abs(M0_theory)
assert rel_err < 1e-3, f"Newton的反力應該在合理誤差內對得上線性解{M0_theory}, 實際={M0_newton}, 相對誤差={rel_err}"
assert conv is True, "有沉陷的模型不應該讓Newton疊代變得不收斂"
print(f"PASS: uy1={uy1}(精確等於-0.02), Newton M0={M0_newton:.4f}(線性解{M0_theory:.4f}), 相對誤差={rel_err:.2e}\n")


print("=== 案例2: 支座沉陷在corotational單步裡也正確生效(位移精確等於指定值) ===")
f_os = three_node_settlement_model()
hs2 = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0),
       1: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u2, F2, ev2, hsf2, conv2, u_full2 = run_pushover_corotational_oneshot(
    f_os, hs2, prescribed_dofs=[f_os.dofs_of(1)[1]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f_os.dofs_of(0)[0]], include_final_displacement=True,
)
uy1_os = u_full2[f_os.dofs_of(1)[1]]
assert abs(uy1_os - (-0.02)) < 1e-12, f"支承1的位移應該精確等於指定的沉陷值-0.02, 實際={uy1_os}"
M0_os = _newton_reaction_moment(f_os, hs2, u_full2)
rel_err_os = abs(M0_os - M0_theory) / abs(M0_theory)
assert rel_err_os < 1e-3, f"corotational單步的反力應該在合理誤差內對得上線性解, 相對誤差={rel_err_os}"
print(f"PASS: uy1={uy1_os}(精確等於-0.02), M0={M0_os:.4f}(線性解{M0_theory:.4f}), 相對誤差={rel_err_os:.2e}\n")


print("=== 案例3: 沒有沉陷時(所有支承都是0或None), 完全不影響既有結果 ===")


def three_node_no_settlement():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_node(2, 2 * L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.add_member(1, node_i=1, node_j=2, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.support(0, ux=0.0, uy=0.0, rot=0.0)
    f.support(1, ux=0.0, uy=0.0, rot=None)
    f.support(2, ux=0.0, uy=0.0, rot=None)
    f.point_load(1, fx=0, fy=-50e3, m=0)
    return f


f3a = three_node_no_settlement()
hs3a = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0),
        1: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
_, _, _, _, _, u_full3a = run_pushover_newton(
    f3a, hs3a, prescribed_dofs=[f3a.dofs_of(1)[1]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f3a.dofs_of(0)[0]], tol=1e-6, include_final_displacement=True,
)

f3b = three_node_no_settlement()
hs3b = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0),
        1: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
_, _, _, _, _, u_full3b = run_pushover_newton(
    f3b, hs3b, prescribed_dofs=[f3b.dofs_of(1)[1]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f3b.dofs_of(0)[0]], tol=1e-6, include_final_displacement=True,
)
assert np.allclose(u_full3a, u_full3b, atol=1e-12), "沒有沉陷時, 結果應該完全一致(逐位元/浮點層級)"
print("PASS: 沒有沉陷的模型, 結果不受這次修改影響\n")

print("PASS: frame2d 支座沉陷(Newton/corotational)所有案例通過")


print("=== 案例4: 溫度效應(thermal_load)在Newton裡正確生效(懸臂樑獨立驗證,")
print("    自由端應該精確位移alpha*delta_T*L, 沒有任何內力) ===")
alpha, delta_T = 1.2e-5, 30.0
theory_ux1 = alpha * delta_T * L


def cantilever_thermal():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A, alpha=alpha)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.thermal_load(0, delta_T=delta_T)
    f.fix(0)
    return f


f4 = cantilever_thermal()
hs4 = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u4, F4, ev4, hsf4, conv4, u_full4 = run_pushover_newton(
    f4, hs4, prescribed_dofs=[f4.dofs_of(1)[1]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f4.dofs_of(0)[0]], tol=1e-6, include_final_displacement=True,
)
ux1_4 = u_full4[f4.dofs_of(1)[0]]
assert abs(ux1_4 - theory_ux1) < 1e-6, f"Newton裡懸臂樑軸向熱伸長應該精確等於{theory_ux1}, 實際={ux1_4}"
print(f"PASS: Newton自由端熱伸長={ux1_4:.8f}(理論{theory_ux1:.8f})\n")

print("=== 案例5: 溫度效應在co-rotational單步裡也正確生效 ===")
f5 = cantilever_thermal()
hs5 = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u5, F5, ev5, hsf5, conv5, u_full5 = run_pushover_corotational_oneshot(
    f5, hs5, prescribed_dofs=[f5.dofs_of(1)[1]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f5.dofs_of(0)[0]], include_final_displacement=True,
)
ux1_5 = u_full5[f5.dofs_of(1)[0]]
assert abs(ux1_5 - theory_ux1) < 1e-6, f"co-rotational單步裡懸臂樑軸向熱伸長應該精確等於{theory_ux1}, 實際={ux1_5}"
print(f"PASS: co-rotational單步自由端熱伸長={ux1_5:.8f}(理論{theory_ux1:.8f})\n")

print("PASS: frame2d 溫度效應(Newton/corotational)所有案例通過")
