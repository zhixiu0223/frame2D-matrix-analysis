"""
驗證frame2d.model.Frame2D.thermal_load()——桿件的溫度效應載重, 分成
軸向(均勻溫度變化)+彎曲(截面深度方向溫度梯度)兩個獨立分量。

驗證策略: 用懸臂樑(一端自由)——這種邊界條件下, 桿件可以完全自由地
伸縮/彎曲, 是靜定結構, 純粹靠幾何(自由熱應變/熱曲率的積分)就能
獨立算出自由端的位移/轉角, 完全不依賴這條固定端反力公式本身怎麼
實作, 而且理論上桿件內部不會有任何應力/內力(沒有東西在阻止它自由
變形)——這是最直接、最不依賴實作細節的驗證方式。額外用固定-固定樑
交叉驗證M(x)的量值。
"""
import numpy as np
from frame2d import Frame2D, solve

E, A, I, L = 200e9, 1e-2, 8e-5, 4.0
alpha = 1.2e-5   # 典型鋼材熱膨脹係數, 1/°C
depth = 0.3      # 截面深度, m


print("=== 案例1: 懸臂樑軸向均勻升溫, 應該可以自由伸長, 沒有任何應力 ===")
delta_T = 30.0
f1 = Frame2D()
f1.add_node(0, 0, 0)
f1.add_node(1, L, 0)
f1.add_section('sec', E=E, I=I, A=A, alpha=alpha)
f1.add_member(0, node_i=0, node_j=1, section='sec')
f1.thermal_load(0, delta_T=delta_T)
f1.fix(0)
r1 = solve(f1)
ux1 = r1.displacements[f1.dofs_of(1)[0]]
theory_ux1 = alpha * delta_T * L
assert abs(ux1 - theory_ux1) < 1e-9, f"自由端軸向熱伸長應該精確等於alpha*delta_T*L={theory_ux1}, 實際={ux1}"
N1 = r1.member_results[0].end_forces_local[3]
assert abs(N1) < 1e-6, f"靜定自由伸縮不應該有任何軸力, 實際Fx2={N1}"
print(f"PASS: 自由端伸長={ux1:.8f}(理論{theory_ux1:.8f}), 軸力={N1:.2e}(應接近0)\n")


print("=== 案例2: 懸臂樑受截面溫度梯度, 應該可以自由彎曲, 沒有任何彎矩 ===")
dT_top, dT_bottom = 40.0, 10.0
kappa_theory = alpha * (dT_top - dT_bottom) / depth
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L, 0)
f2.add_section('sec', E=E, I=I, A=A, alpha=alpha, depth=depth)
f2.add_member(0, node_i=0, node_j=1, section='sec')
f2.thermal_load(0, delta_T_top=dT_top, delta_T_bottom=dT_bottom)
f2.fix(0)
r2 = solve(f2)
rot1 = r2.displacements[f2.dofs_of(1)[2]]
uy1 = r2.displacements[f2.dofs_of(1)[1]]
theory_rot1 = kappa_theory * L
theory_uy1 = kappa_theory * L**2 / 2
assert abs(rot1 - theory_rot1) < 1e-9, f"自由端轉角應該精確等於kappa*L={theory_rot1}, 實際={rot1}"
assert abs(uy1 - theory_uy1) < 1e-9, f"自由端撓度應該精確等於kappa*L^2/2={theory_uy1}, 實際={uy1}"
M2 = r2.member_results[0].end_forces_local[5]
assert abs(M2) < 1e-6, f"靜定自由彎曲不應該有任何彎矩, 實際M2={M2}"
print(f"PASS: 自由端轉角={rot1:.8f}(理論{theory_rot1:.8f}), 撓度={uy1:.8f}(理論{theory_uy1:.8f}), M={M2:.2e}(應接近0)\n")


print("=== 案例3: 固定-固定樑受溫度梯度, M(x)應該全程是常數=EI*kappa_thermal ===")
from frame2d.postprocess import member_internal_forces

f3 = Frame2D()
f3.add_node(0, 0, 0)
f3.add_node(1, L, 0)
f3.add_section('sec', E=E, I=I, A=A, alpha=alpha, depth=depth)
f3.add_member(0, node_i=0, node_j=1, section='sec')
f3.thermal_load(0, delta_T_top=dT_top, delta_T_bottom=dT_bottom)
f3.fix(0)
f3.fix(1)
r3 = solve(f3)
x3, N3, V3, M3 = member_internal_forces(f3, r3, 0, n=9)
M_theory3 = E * I * kappa_theory
assert np.max(np.abs(M3 - (-M_theory3))) < 1e-3, f"M(x)應該全程等於-EI*kappa={-M_theory3}, 實際={M3}"
assert np.max(np.abs(V3)) < 1e-6, "均勻溫度梯度不應該對V(x)有任何貢獻(對稱載重, 沒有淨側向力)"
print(f"PASS: M(x)={M3}\n     理論(全程常數)={-M_theory3}\n")


print("=== 案例4: 固定-固定樑軸向均勻升溫, 軸力應該全程是常數=-EA*alpha*delta_T(升溫受束制=壓力) ===")
f4 = Frame2D()
f4.add_node(0, 0, 0)
f4.add_node(1, L, 0)
f4.add_section('sec', E=E, I=I, A=A, alpha=alpha)
f4.add_member(0, node_i=0, node_j=1, section='sec')
f4.thermal_load(0, delta_T=delta_T)
f4.fix(0)
f4.fix(1)
r4 = solve(f4)
x4, N4, V4, M4 = member_internal_forces(f4, r4, 0, n=5)
N_theory4 = -E * A * alpha * delta_T   # N=-Fx1慣例(拉力為正), 升溫受束制應該是壓力(負值)
assert np.max(np.abs(N4 - N_theory4)) < 1e-3, f"N(x)應該全程等於-EA*alpha*delta_T={N_theory4}, 實際={N4}"
print(f"PASS: N(x)={N4}(理論{N_theory4}, 負值=壓力, 升溫受束制符合預期)\n")


print("=== 案例5: 沒有設定alpha的斷面用delta_T, 應該明確拒絕(不會靜默忽略) ===")
f5 = Frame2D()
f5.add_node(0, 0, 0)
f5.add_node(1, L, 0)
f5.add_section('sec', E=E, I=I, A=A)   # 沒有alpha
f5.add_member(0, node_i=0, node_j=1, section='sec')
f5.thermal_load(0, delta_T=delta_T)
f5.fix(0)
try:
    solve(f5)
    assert False, "沒設定alpha卻用delta_T應該要raise"
except ValueError as e:
    assert 'alpha' in str(e)
print("PASS: 沒設定alpha正確拒絕\n")


print("=== 案例6: 有溫度梯度但斷面沒設定depth, 應該明確拒絕 ===")
f6 = Frame2D()
f6.add_node(0, 0, 0)
f6.add_node(1, L, 0)
f6.add_section('sec', E=E, I=I, A=A, alpha=alpha)   # 沒有depth
f6.add_member(0, node_i=0, node_j=1, section='sec')
f6.thermal_load(0, delta_T_top=dT_top, delta_T_bottom=dT_bottom)
f6.fix(0)
try:
    solve(f6)
    assert False, "有溫度梯度但沒設定depth應該要raise"
except ValueError as e:
    assert 'depth' in str(e)
print("PASS: 沒設定depth正確拒絕\n")


print("=== 案例7: truss元素承受溫度梯度應該明確拒絕(沒有彎曲勁度), 但均勻溫度變化可以 ===")
f7 = Frame2D()
f7.add_node(0, 0, 0)
f7.add_node(1, L, 0)
f7.add_section('sec', E=E, I=I, A=A, alpha=alpha, depth=depth)
f7.add_member(0, node_i=0, node_j=1, section='sec', member_type='truss')
f7.thermal_load(0, delta_T_top=dT_top, delta_T_bottom=dT_bottom)
f7.pin(0)
f7.roller_y(1)
try:
    solve(f7)
    assert False, "truss元素承受溫度梯度應該要raise"
except ValueError as e:
    assert 'truss' in str(e)
print("PASS: truss元素承受溫度梯度正確拒絕\n")

f7b = Frame2D()
f7b.add_node(0, 0, 0)
f7b.add_node(1, L, 0)
f7b.add_section('sec', E=E, I=I, A=A, alpha=alpha)
f7b.add_member(0, node_i=0, node_j=1, section='sec', member_type='truss')
f7b.thermal_load(0, delta_T=delta_T)
f7b.pin(0)
f7b.roller_y(1)
r7b = solve(f7b)  # 不應該raise, truss可以承受均勻軸向熱效應
print("PASS: truss元素承受均勻溫度變化正常求解(不受彎曲熱效應限制)\n")


print("=== 案例8: thermal_load()完全沒給任何參數時應該明確拒絕 ===")
f8 = Frame2D()
f8.add_node(0, 0, 0)
f8.add_node(1, L, 0)
f8.add_section('sec', E=E, I=I, A=A)
f8.add_member(0, node_i=0, node_j=1, section='sec')
try:
    f8.thermal_load(0)
    assert False, "沒給任何參數應該要raise"
except ValueError as e:
    assert 'thermal_load' in str(e)
print("PASS: 沒給任何參數正確拒絕\n")


print("=== 案例9: 沒有使用thermal_load時(預設空list), 完全不影響既有結果 ===")
f9a = Frame2D()
f9a.add_node(0, 0, 0)
f9a.add_node(1, L, 0)
f9a.add_section('sec', E=E, I=I, A=A)
f9a.add_member(0, node_i=0, node_j=1, section='sec')
f9a.point_load(1, fx=0, fy=-50.0, m=0)
f9a.fix(0)
r9a = solve(f9a)

f9b = Frame2D()
f9b.add_node(0, 0, 0)
f9b.add_node(1, L, 0)
f9b.add_section('sec', E=E, I=I, A=A)
f9b.add_member(0, node_i=0, node_j=1, section='sec')
f9b.point_load(1, fx=0, fy=-50.0, m=0)
f9b.fix(0)
r9b = solve(f9b)
assert np.array_equal(r9a.displacements, r9b.displacements), "空的thermal_loads不應該影響既有結果"
print("PASS: 空的thermal_loads列表對結果完全沒有影響\n")

print("PASS: frame2d thermal_load所有案例通過")
