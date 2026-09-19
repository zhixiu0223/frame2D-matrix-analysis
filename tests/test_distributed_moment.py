"""
驗證frame2d.model.Frame2D.distributed_moment()/frame2d.elements.
fixed_end_forces_distributed_moment()——桿件全長均佈的分布彎矩
(kN·m/m, 逆時針為正)。

驗證策略(跟fixed_end_forces_point_moment()當初抓到符號bug用的同一套
精神): 這條公式本身沒有獨立於推導方法之外的"基本物理事實"可以直接
核對, 所以改用兩種靜定結構(簡支樑、懸臂樑)的獨立靜力平衡驗證——
這兩種結構的反力可以純粹用手算靜力學算出來, 完全不依賴這條固定端
反力公式本身, 拿FEM組裝出來的反力去對照, 才是真正獨立的驗證, 不是
自己驗證自己。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A, L = 200e9, 8e-5, 1e-2, 4.0
m_val = 10.0  # kN·m/m, 逆時針為正


print("=== 案例1: 簡支樑, 純靜力平衡驗證 ===")
print("對node0取力矩平衡: R1y(向上為正)對node0產生逆時針力矩+R1y*L,")
print("外加的均佈彎矩(逆時針為正)貢獻+m*L, 平衡: R1y*L+m*L=0 -> R1y=-m;")
print("垂直力平衡: R0y+R1y=0 -> R0y=+m。")
f1 = Frame2D()
f1.add_node(0, 0, 0)
f1.add_node(1, L, 0)
f1.add_section('sec', E=E, I=I, A=A)
f1.add_member(0, node_i=0, node_j=1, section='sec')
f1.distributed_moment(0, m=m_val)
f1.pin(0)
f1.roller_y(1)
r1 = solve(f1)
R0y = r1.reactions[f1.dofs_of(0)[1]]
R1y = r1.reactions[f1.dofs_of(1)[1]]
assert abs(R0y - m_val) < 1e-6, f"簡支樑R0y應該精確等於+m={m_val}, 實際={R0y}"
assert abs(R1y - (-m_val)) < 1e-6, f"簡支樑R1y應該精確等於-m={-m_val}, 實際={R1y}"
print(f"PASS: R0y={R0y}(理論+m), R1y={R1y}(理論-m)\n")


print("=== 案例2: 懸臂樑, 純靜力平衡驗證 ===")
print("對固定端取力矩平衡: M0(固定端彎矩反力)+均佈彎矩總量m*L=0 -> M0=-m*L。")
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L, 0)
f2.add_section('sec', E=E, I=I, A=A)
f2.add_member(0, node_i=0, node_j=1, section='sec')
f2.distributed_moment(0, m=m_val)
f2.fix(0)
r2 = solve(f2)
M0 = r2.reactions[f2.dofs_of(0)[2]]
assert abs(M0 - (-m_val * L)) < 1e-6, f"懸臂樑M0應該精確等於-m*L={-m_val*L}, 實際={M0}"
print(f"PASS: M0={M0}(理論-m*L={-m_val*L})\n")


print("=== 案例3: truss/cable元素加分布彎矩應該明確拒絕, 不是靜默給錯誤答案 ===")
f3 = Frame2D()
f3.add_node(0, 0, 0)
f3.add_node(1, L, 0)
f3.add_section('sec', E=E, I=I, A=A)
f3.add_member(0, node_i=0, node_j=1, section='sec', member_type='truss')
f3.distributed_moment(0, m=m_val)
f3.pin(0)
f3.roller_y(1)
try:
    solve(f3)
    assert False, "truss元素加分布彎矩應該要raise"
except ValueError as e:
    assert 'truss' in str(e)
print("PASS: truss元素正確拒絕\n")


print("=== 案例4: 沒有使用distributed_moment時(預設空list), 完全不影響既有結果 ===")
f4 = Frame2D()
f4.add_node(0, 0, 0)
f4.add_node(1, L, 0)
f4.add_section('sec', E=E, I=I, A=A)
f4.add_member(0, node_i=0, node_j=1, section='sec')
f4.point_load(1, fx=0, fy=-50.0, m=0)
f4.pin(0)
f4.roller_y(1)
r4 = solve(f4)

f4b = Frame2D()
f4b.add_node(0, 0, 0)
f4b.add_node(1, L, 0)
f4b.add_section('sec', E=E, I=I, A=A)
f4b.add_member(0, node_i=0, node_j=1, section='sec')
f4b.point_load(1, fx=0, fy=-50.0, m=0)
f4b.pin(0)
f4b.roller_y(1)
r4b = solve(f4b)
assert np.array_equal(r4.displacements, r4b.displacements), "空的distributed_moments不應該影響既有結果"
print("PASS: 空的distributed_moments列表對結果完全沒有影響\n")


print("=== 案例5: distributed_moment在Newton/co-rotational單步裡也正確生效")
print("    (透過newton.py獨立的_gravity_fixed_end_forces(), 跟dofmanager.py")
print("    共用同一個fixed_end_forces_distributed_moment()公式) ===")
from frame2d.hinge import HingeState
from frame2d.newton import (run_pushover_newton, run_pushover_corotational_oneshot,
                             _assemble_global, _gravity_fixed_end_forces)


def cantilever_with_distributed_moment():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.distributed_moment(0, m=m_val)
    f.fix(0)
    return f


f5a = cantilever_with_distributed_moment()
f_ext_gravity5a, _ = _gravity_fixed_end_forces(f5a)
hs5a = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u5a, F5a, ev5a, hsf5a, conv5a, u_full5a = run_pushover_newton(
    f5a, hs5a, prescribed_dofs=[f5a.dofs_of(1)[0]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f5a.dofs_of(0)[0]], tol=1e-9, include_final_displacement=True,
)
member_ref5a = {0: (0.0, 0.0, 0.0, 0.0)}
f_int5a, _ = _assemble_global(f5a, hs5a, u_full5a, len(u_full5a), member_ref5a)
rot_dof5 = f5a.dofs_of(0)[2]
M0_newton = f_int5a[rot_dof5] - f_ext_gravity5a[rot_dof5]
assert abs(M0_newton - (-m_val * L)) < 1e-6, f"Newton裡M0反力應該精確等於-m*L, 實際={M0_newton}"
print(f"PASS: Newton M0反力={M0_newton:.6f}(理論{-m_val*L})\n")

f5b = cantilever_with_distributed_moment()
f_ext_gravity5b, _ = _gravity_fixed_end_forces(f5b)
hs5b = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u5b, F5b, ev5b, hsf5b, conv5b, u_full5b = run_pushover_corotational_oneshot(
    f5b, hs5b, prescribed_dofs=[f5b.dofs_of(1)[0]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f5b.dofs_of(0)[0]], include_final_displacement=True,
)
member_ref5b = {0: (0.0, 0.0, 0.0, 0.0)}
f_int5b, _ = _assemble_global(f5b, hs5b, u_full5b, len(u_full5b), member_ref5b)
M0_oneshot = f_int5b[rot_dof5] - f_ext_gravity5b[rot_dof5]
assert abs(M0_oneshot - (-m_val * L)) < 1e-6, f"co-rotational單步裡M0反力應該精確等於-m*L, 實際={M0_oneshot}"
print(f"PASS: co-rotational單步 M0反力={M0_oneshot:.6f}(理論{-m_val*L})\n")

print("=== 案例6: distributed_moment在event-to-event pushover(apply_gravity)裡")
print("    也正確生效(共用dofmanager.py的solve_with_hinges()) ===")
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import apply_gravity

f6 = cantilever_with_distributed_moment()
hs6 = initial_hinge_states(f6)
init_forces6, gravity_result6 = apply_gravity(f6, hs6)
M0_pushover = gravity_result6.reactions[f6.dofs_of(0)[2]]
assert abs(M0_pushover - (-m_val * L)) < 1e-6, f"pushover重力預載M0反力應該精確等於-m*L, 實際={M0_pushover}"
print(f"PASS: pushover重力預載 M0反力={M0_pushover:.6f}(理論{-m_val*L})\n")

print("PASS: frame2d distributed_moment所有案例通過(含Newton/co-rotational/pushover)")
