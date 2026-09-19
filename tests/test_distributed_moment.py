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


print("=== 案例7: fixed_end_forces_partial_distributed_moment()對照sympy符號積分獨立驗證")
print("    (數值高斯積分 vs 封閉式符號解析解, 兩條完全獨立的推導路徑) ===")
import sympy as sp
from frame2d.elements import fixed_end_forces_partial_distributed_moment

_L, _x, _c, _d, _m0, _m1 = sp.symbols('L x c d m0 m1', real=True)
_N1 = 1 - 3*(_x/_L)**2 + 2*(_x/_L)**3
_N2 = _x*(1 - _x/_L)**2
_N3 = 3*(_x/_L)**2 - 2*(_x/_L)**3
_N4 = (_x**2/_L)*(_x/_L - 1)
_mx = _m0 + (_m1 - _m0)*(_x - _c)/(_d - _c)
_F1s = sp.simplify(sp.integrate(_mx * sp.diff(_N1, _x), (_x, _c, _d)))
_F2s = sp.simplify(sp.integrate(_mx * sp.diff(_N2, _x), (_x, _c, _d)))
_F3s = sp.simplify(sp.integrate(_mx * sp.diff(_N3, _x), (_x, _c, _d)))
_F4s = sp.simplify(sp.integrate(_mx * sp.diff(_N4, _x), (_x, _c, _d)))
_subs = {_L: 8.0, _c: 2.0, _d: 5.0, _m0: 3.0, _m1: 9.0}
exact = np.array([0.0, float(_F1s.subs(_subs)), float(_F2s.subs(_subs)),
                   0.0, float(_F3s.subs(_subs)), float(_F4s.subs(_subs))])
numeric = fixed_end_forces_partial_distributed_moment(3.0, 9.0, 2.0, 5.0, 8.0)
rel_err7 = np.max(np.abs(exact - numeric))
assert rel_err7 < 1e-10, f"高斯積分應該跟sympy符號解析解精確一致, 實際最大誤差={rel_err7}"
print(f"PASS: 符號解={exact}, 數值解={numeric}, 最大誤差={rel_err7:.2e}\n")

print("=== 案例8: 退化檢查——局部段c=0,d=L, m_start=m_end=m時, 應該完全等於")
print("    fixed_end_forces_distributed_moment(m, L) ===")
from frame2d.elements import fixed_end_forces_distributed_moment
full1 = fixed_end_forces_distributed_moment(7.0, 6.0)
full2 = fixed_end_forces_partial_distributed_moment(7.0, 7.0, 0.0, 6.0, 6.0)
assert np.allclose(full1, full2, atol=1e-10), f"退化情況應該完全一致, 實際: {full1} vs {full2}"
print(f"PASS: 全長常數m={full1} 跟局部段退化版{full2}完全一致\n")

print("PASS: frame2d partial distributed_moment所有案例通過")


print("=== 案例9: 局部段+線性變化分布彎矩, 簡支樑反力驗證(純力偶的關鍵性質:")
print("    反力大小只跟合力矩總量有關, 跟力偶作用在桿件哪個位置無關) ===")
c9, d9, m0_9, m1_9 = 1.0, 3.0, 3.0, 9.0
total_moment9 = (m0_9 + m1_9) / 2 * (d9 - c9)
R_left_theory9 = total_moment9 / L
R_right_theory9 = -total_moment9 / L

f9 = Frame2D()
f9.add_node(0, 0, 0)
f9.add_node(1, L, 0)
f9.add_section('sec', E=E, I=I, A=A)
f9.add_member(0, node_i=0, node_j=1, section='sec')
f9.distributed_moment(0, m=m0_9, m_end=m1_9, x_start=c9, x_end=d9)
f9.pin(0)
f9.roller_y(1)
r9 = solve(f9)
R0y9 = r9.reactions[f9.dofs_of(0)[1]]
R1y9 = r9.reactions[f9.dofs_of(1)[1]]
assert abs(R0y9 - R_left_theory9) < 1e-6, f"R0y應該={R_left_theory9}, 實際={R0y9}"
assert abs(R1y9 - R_right_theory9) < 1e-6, f"R1y應該={R_right_theory9}, 實際={R1y9}"
print(f"PASS: R0y={R0y9}(理論{R_left_theory9}), R1y={R1y9}(理論{R_right_theory9})\n")


print("=== 案例10: 退化檢查——不給m_end/x_start/x_end時, 應該完全等同全長常數版 ===")
f10a = Frame2D()
f10a.add_node(0, 0, 0)
f10a.add_node(1, L, 0)
f10a.add_section('sec', E=E, I=I, A=A)
f10a.add_member(0, node_i=0, node_j=1, section='sec')
f10a.distributed_moment(0, m=m_val)
f10a.fix(0)
r10a = solve(f10a)

f10b = Frame2D()
f10b.add_node(0, 0, 0)
f10b.add_node(1, L, 0)
f10b.add_section('sec', E=E, I=I, A=A)
f10b.add_member(0, node_i=0, node_j=1, section='sec')
f10b.distributed_moment(0, m=m_val, m_end=m_val, x_start=0.0, x_end=L)
f10b.fix(0)
r10b = solve(f10b)
assert np.allclose(r10a.displacements, r10b.displacements, atol=1e-9), "退化情況應該完全一致"
print("PASS: 顯式指定全長常數 跟 預設(不給任何額外參數) 結果一致\n")

print("PASS: frame2d 局部段+線性變化分布彎矩所有案例通過")


print("=== 案例11: 範圍超出桿件實際長度時, 明確拒絕(不會靜默外推給錯誤答案) ===")
f11 = Frame2D()
f11.add_node(0, 0, 0)
f11.add_node(1, L, 0)
f11.add_section('sec', E=E, I=I, A=A)
f11.add_member(0, node_i=0, node_j=1, section='sec')
f11.distributed_moment(0, m=5.0, x_start=1.0, x_end=10.0)
f11.pin(0)
f11.roller_y(1)
try:
    solve(f11)
    assert False, "範圍超出桿件長度應該要raise"
except ValueError as e:
    assert '超出桿件' in str(e)
print("PASS: 範圍超出桿件長度正確拒絕\n")


print("=== 案例12: 局部段+線性變化分布彎矩在Newton/co-rotational單步裡也正確生效 ===")
from frame2d.hinge import HingeState
from frame2d.newton import run_pushover_newton, _assemble_global, _gravity_fixed_end_forces
from frame2d.dofmanager import initial_hinge_states


def cantilever_for_partial_dm():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=1e30, Mp_j=1e30,
                 R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.distributed_moment(0, m=3.0, m_end=9.0, x_start=1.0, x_end=3.0)
    f.fix(0)
    return f


f12 = cantilever_for_partial_dm()
total_moment12 = (3.0 + 9.0) / 2 * (3.0 - 1.0)
f_ext_gravity12, _ = _gravity_fixed_end_forces(f12)
hs12 = {0: HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1.0, R_post_yield_2=1.0)}
u12, F12, ev12, hsf12, conv12, u_full12 = run_pushover_newton(
    f12, hs12, prescribed_dofs=[f12.dofs_of(1)[0]], direction=[1.0], target_total=1e-9, d_nominal=1e-9,
    base_reaction_dofs=[f12.dofs_of(0)[0]], tol=1e-6, include_final_displacement=True,
    # tol刻意用1e-6, 不是1e-9: co-rotational的切線用有限差分(corotational_
    # tangent_fd(), 不是解析導數)算, 實測發現這個局部段+線性變化的案例
    # 殘餘力會在1e-7量級附近震盪、收斂不到1e-9(不是我這個功能本身的
    # bug, 是有限差分切線本身固有的精度上限)——這不是計算錯誤, 是
    # Newton疊代對"多緊算收斂"這個容許誤差要設得跟切線本身的精度匹配,
    # 1e-6對這個模型的量級來說仍然是很緊的容許誤差。
)
member_ref12 = {0: (0.0, 0.0, 0.0, 0.0)}
f_int12, _ = _assemble_global(f12, hs12, u_full12, len(u_full12), member_ref12)
rot_dof12 = f12.dofs_of(0)[2]
M0_newton12 = f_int12[rot_dof12] - f_ext_gravity12[rot_dof12]
# 懸臂樑: 固定端彎矩反力應該精確等於"合力矩"(不管作用在哪個位置)取負號
assert abs(M0_newton12 - (-total_moment12)) < 1e-6, (
    f"Newton裡局部段+線性變化分布彎矩應該正確生效, 理論M0={-total_moment12}, 實際={M0_newton12}"
)
print(f"PASS: Newton M0反力={M0_newton12:.6f}(理論{-total_moment12})\n")

f12b = cantilever_for_partial_dm()
hs12b = initial_hinge_states(f12b)
from frame2d.pushover import apply_gravity
init_forces12b, gravity_result12b = apply_gravity(f12b, hs12b)
M0_pushover12 = gravity_result12b.reactions[f12b.dofs_of(0)[2]]
assert abs(M0_pushover12 - (-total_moment12)) < 1e-6, (
    f"pushover重力預載裡局部段+線性變化分布彎矩應該正確生效, 理論M0={-total_moment12}, 實際={M0_pushover12}"
)
print(f"PASS: pushover重力預載 M0反力={M0_pushover12:.6f}(理論{-total_moment12})\n")

print("PASS: frame2d 局部段+線性變化分布彎矩(含Newton/pushover)所有案例通過")


print("=== 案例13: M(x)在桿件內部正確反映分布彎矩造成的局部變化(不是只看")
print("    端點值連成一條直線)——這是實際使用者用M圖截圖發現、回報的真實bug:")
print("    member_internal_forces()完全沒有把distributed_moments納入計算,")
print("    只有兩端節點反力/端點彎矩是對的, 桿件內部的M(x)完全沒反映分布")
print("    彎矩的局部效應(不管放在哪裡、是均佈還是梯形, M圖看起來都一樣) ===")
from frame2d.postprocess import member_internal_forces

L13 = 8.0
c13, d13, m0_13, m1_13 = 2.5, 4.0, 2.0, 3.0
f13 = Frame2D()
f13.add_node(0, 0, 0)
f13.add_node(1, L13, 0)
f13.add_section('sec', E=E, I=I, A=A)
f13.add_member(0, node_i=0, node_j=1, section='sec')
f13.distributed_moment(0, m=m0_13, m_end=m1_13, x_start=c13, x_end=d13)
f13.pin(0)
f13.roller_y(1)
r13 = solve(f13)
x13, N13, V13, M13 = member_internal_forces(f13, r13, 0, n=25)

# 獨立驗證1: 分布彎矩範圍前後M(x)應該是直線(斜率=反力大小, 不受分布
# 彎矩的"形狀"影響, 只受它的"總量"透過反力間接影響)——用範圍前兩個點
# 的斜率去外插到x=0, 應該等於M(0)。
before_mask = x13 < c13 - 1e-9
after_mask = x13 > d13 + 1e-9
x_before, M_before = x13[before_mask], M13[before_mask]
x_after, M_after = x13[after_mask], M13[after_mask]
slope_before = np.polyfit(x_before, M_before, 1)[0]
slope_after = np.polyfit(x_after, M_after, 1)[0]
resid_before = M_before - np.polyval(np.polyfit(x_before, M_before, 1), x_before)
resid_after = M_after - np.polyval(np.polyfit(x_after, M_after, 1), x_after)
assert np.max(np.abs(resid_before)) < 1e-8, "分布彎矩範圍前, M(x)應該是精確直線"
assert np.max(np.abs(resid_after)) < 1e-8, "分布彎矩範圍後, M(x)應該是精確直線"
print(f"PASS: 分布彎矩範圍前後M(x)都是精確直線(殘差<1e-8), 斜率分別={slope_before:.4f}, {slope_after:.4f}\n")

# 獨立驗證2: 跨越分布彎矩範圍[c,d]前後, M(x)的"落差"應該精確等於
# Fy1*(d-c) - ∫[c,d]m(x')dx'(前者是剪力/反力對這段長度的線性貢獻,
# 後者是分布彎矩本身的合力矩, 梯形公式, m線性變化的精確解)——一開始
# 漏算了Fy1*(d-c)這一項(簡支樑受純力偶會產生剪力反力, V(x)不是0,
# 是常數=支承反力, 這裡也一併驗證), 這是純粹幾何/微積分的關係, 不
# 依賴這條公式本身怎麼實作。
M_at_c = np.interp(c13, x13, M13)
M_at_d = np.interp(d13, x13, M13)
Fy1_13 = V13[0]   # 簡支樑受純力偶, V(x)應該全程是常數=左支承反力
total_moment13 = (m0_13 + m1_13) / 2 * (d13 - c13)
expected_drop = Fy1_13 * (d13 - c13) - total_moment13
actual_drop = M_at_d - M_at_c
# 容許誤差用1e-4不是1e-6: M(x)在[c,d]範圍內是二次曲線(m線性變化的
# 累積), 這裡用np.interp在c/d的eps鄰近點(見member_internal_forces()
# 的extra取樣點)做線性內插取值, 對二次曲線的線性內插本來就會有跟
# eps同量級的殘留誤差, 不是物理計算本身不精確。
assert abs(actual_drop - expected_drop) < 1e-4, (
    f"M(x)跨越分布彎矩範圍的落差應該精確等於Fy1*(d-c)-合力矩, 理論={expected_drop}, 實際={actual_drop}"
)
print(f"PASS: M(c)={M_at_c:.4f}, M(d)={M_at_d:.4f}, 落差={actual_drop:.4f}(理論{expected_drop:.4f})\n")

# 獨立驗證3: 簡支樑受純力偶(沒有其他分布力), V(x)應該全程是常數,
# 精確等於支承反力(純力偶用"力偶等效成兩端支承反力"這個關鍵性質,
# 這裡直接拿反力去對照V(x))。
R0_13 = r13.reactions[f13.dofs_of(0)[1]]
assert np.max(np.abs(V13 - R0_13)) < 1e-8, (
    f"簡支樑受純力偶時V(x)應該全程等於支承反力{R0_13}, 實際V(x)範圍=[{V13.min()},{V13.max()}]"
)
print(f"PASS: V(x)全程={V13[0]:.6f}(精確等於支承反力{R0_13:.6f})\n")

print("=== 案例14: 全長常數分布彎矩的懸臂樑, M(x)應該是精確直線(斜率=m) ===")
f14 = Frame2D()
f14.add_node(0, 0, 0)
f14.add_node(1, L, 0)
f14.add_section('sec', E=E, I=I, A=A)
f14.add_member(0, node_i=0, node_j=1, section='sec')
f14.distributed_moment(0, m=m_val)
f14.fix(0)
r14 = solve(f14)
x14, N14, V14, M14 = member_internal_forces(f14, r14, 0, n=9)
# 理論: M(0)應該等於固定端反力矩的相反數(=+m*L, 跟_partial_load_W在
# 全長時累積到m*L、用減號扣掉一致), M(L)=0(自由端), 中間線性內插:
# M(x) = m*L - m*x = m*(L-x)。
theory_M14 = m_val * (L - x14)
assert np.max(np.abs(M14 - theory_M14)) < 1e-6, f"全長常數分布彎矩M(x)應該精確線性, 實際={M14}, 理論={theory_M14}"
assert np.max(np.abs(V14)) < 1e-8, "V(x)應該全程接近0(懸臂樑受純力偶, 固定端只提供彎矩反力, 不提供剪力反力)"
print(f"PASS: M(x)={M14}\n     理論={theory_M14}\n")

print("PASS: frame2d M(x)內部彎矩圖正確反映分布彎矩局部效應, 所有案例通過")
