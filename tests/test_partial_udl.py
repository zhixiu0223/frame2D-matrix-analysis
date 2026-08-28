"""
驗證案例10: 桿件局部段(不是整根桿件)均佈載重 (Load System v2 Phase 2)

驗證策略跟Phase1一致: 優先用簡支梁的決定性反力(純靜力學, 完全不依賴
fixed_end_forces_partial_udl公式本身對不對)當基準。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.elements import fixed_end_forces_udl, fixed_end_forces_partial_udl

E, I, A = 200e6, 8e-5, 1e-2


def check(label, fem, exact, tol=1e-8):
    rel = abs(fem - exact) / max(abs(exact), 1e-12)
    print(f"  {label}: FEM={fem:.8f}  解析解={exact:.8f}  rel_err={rel:.2e}")
    assert rel < tol, f"{label} 不吻合"


# ---- 案例A: 退化情況(c=0,d=L)應該精確等於既有的全長UDL公式 ----
print("=== 案例A: 局部段公式在c=0,d=L退化情況下應精確等於既有全長UDL公式 ===")
L = 8.0
for w1, w2 in [(-15.0, -15.0), (-10.0, -22.0)]:
    old = fixed_end_forces_udl(w1, w2, L)
    new = fixed_end_forces_partial_udl(w1, w2, 0.0, L, L)
    diff = np.max(np.abs(old - new))
    print(f"  w=({w1},{w2}): max_diff={diff:.2e}")
    assert diff < 1e-9
print("PASS: 退化情況精確吻合\n")


# ---- 案例B: 簡支梁, 局部段均佈載重(不對稱, 不含端點), 對照純靜力學反力 ----
print("=== 案例B: 簡支梁局部段均佈載重, 對照純靜力學決定性反力 ===")
L, w, cB, dB = 12.0, -20.0, 3.0, 7.0   # 載重只在[3,7]之間, 大小20kN/m向下
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.distributed_load(0, w=w, x_start=cB, x_end=dB)
r = solve(f)

# 純靜力學: 總力 = w*(d-c), 作用在形心(c+d)/2, 對節點1取矩算R0
W_total = w * (dB - cB)
centroid = (cB + dB) / 2
R0_exact = -W_total * (L - centroid) / L   # 注意W_total是負值(向下), 反力應為正(向上)
R1_exact = -W_total * centroid / L
check("R0(pin)", r.reactions[f.dofs_of(0)[1]], R0_exact)
check("R1(roller)", r.reactions[f.dofs_of(1)[1]], R1_exact)
print("PASS: 簡支梁局部段均佈載重反力精確吻合純靜力學決定性解\n")


# ---- 案例C: 簡支梁, 局部段"線性變化"(梯形)載重, 對照純靜力學反力 ----
print("=== 案例C: 簡支梁局部段線性變化(梯形)載重, 對照純靜力學決定性反力 ===")
L, w1v, w2v, cC, dC = 10.0, -10.0, -30.0, 2.0, 6.0
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's')
f2.pin(0)
f2.roller_y(1)
f2.distributed_load(0, w=w1v, w_end=w2v, x_start=cC, x_end=dC)
r2 = solve(f2)

# 梯形載重拆成均佈w1v(整段)+三角形(0到w2v-w1v), 分開算合力跟形心再疊加
span = dC - cC
W_uniform = w1v * span
centroid_uniform = (cC + dC) / 2
dw = w2v - w1v
W_tri = 0.5 * dw * span
centroid_tri = cC + 2.0 / 3.0 * span   # 三角形(0在c, 峰值在d)形心距c為2/3*span
W_total = W_uniform + W_tri
moment_about_0 = W_uniform * centroid_uniform + W_tri * centroid_tri
centroid_total = moment_about_0 / W_total

R0_exact = -W_total * (L - centroid_total) / L
R1_exact = -W_total * centroid_total / L
check("R0(pin)", r2.reactions[f2.dofs_of(0)[1]], R0_exact, tol=1e-6)
check("R1(roller)", r2.reactions[f2.dofs_of(1)[1]], R1_exact, tol=1e-6)
print("PASS: 簡支梁局部段梯形載重反力精確吻合純靜力學決定性解\n")


# ---- 案例D: 診斷圖 V(x)/M(x) 在載重範圍外應為常數(局部段外沒有載重) ----
print("=== 案例D: 診斷圖確認載重範圍外V(x)為常數(不受局部段外的載重影響) ===")
from frame2d.postprocess import member_internal_forces
xB, NB, VB, MB = member_internal_forces(f, r, 0, n=41)
# x<cB: V應該等於R0(常數, 還沒進入載重段); x>dB: V應該等於-R1(常數, 離開載重段後不再變化)
margin = 0.05
V_before = VB[xB < cB - margin]
V_after = VB[xB > dB + margin]
assert np.allclose(V_before, V_before[0], atol=1e-8), "載重段前V(x)應為常數"
assert np.allclose(V_after, V_after[0], atol=1e-8), "載重段後V(x)應為常數"
print(f"  載重段前V(x)={V_before[0]:.4f} (應等於R0={r.reactions[f.dofs_of(0)[1]]:.4f})")
print(f"  載重段後V(x)={V_after[0]:.4f}")
assert abs(V_before[0] - r.reactions[f.dofs_of(0)[1]]) < 1e-8
print("PASS: 局部段外V(x)確實維持常數, 沒有被載重段外的區域錯誤影響\n")

print("PASS: Phase 2(局部段均佈載重) 全數驗證通過")
