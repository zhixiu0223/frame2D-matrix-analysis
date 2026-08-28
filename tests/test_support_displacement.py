"""
驗證案例11: 支承強制位移 (沉陷分析, Load System v2 Phase 3)

驗證策略:
案例A: 簡支梁(靜定結構)支承沉陷 -> 靜定結構的支承沉陷只會造成剛體轉動,
       不會產生任何內力/反力(沒有多餘拘束可以「抵抗」沉陷)。這是最乾淨的
       驗證, 完全不依賴任何公式, 純粹是靜定結構的基本性質。
案例B: 固接+滾支承的梁(一次靜不定)支承沉陷 -> 對照傾角變位法的經典公式
       M_A = -3EIΔ/L² (使用者自己sd_framework用的同一種方法, 標準教科書結果)。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 1e-2


def check(label, fem, exact, tol=1e-8, atol=1e-9):
    """exact接近0時改用絕對誤差比較(相對誤差在分母趨近0時會失真)"""
    if abs(exact) < atol:
        err = abs(fem - exact)
        print(f"  {label}: FEM={fem:.10f}  解析解={exact:.10f}  abs_err={err:.2e}")
        assert err < atol, f"{label} 不吻合"
    else:
        rel = abs(fem - exact) / abs(exact)
        print(f"  {label}: FEM={fem:.8f}  解析解={exact:.8f}  rel_err={rel:.2e}")
        assert rel < tol, f"{label} 不吻合"


# ---- 案例A: 簡支梁支承沉陷, 靜定結構應該零反力零內力(純剛體轉動) ----
print("=== 案例A: 簡支梁(靜定)支承沉陷, 不應產生任何反力/內力 ===")
L = 10.0
delta = -0.02   # 節點1沉陷2cm(向下)
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.support(1, uy=delta)   # 滾支承, 但這次uy不是0而是強制沉陷delta
r = solve(f)

check("node1實際位移", r.displacements[f.dofs_of(1)[1]], delta, tol=1e-10)
check("node0反力Rx", r.reactions[f.dofs_of(0)[0]], 0.0, tol=1e-9)
check("node0反力Ry", r.reactions[f.dofs_of(0)[1]], 0.0, tol=1e-9)
check("node1反力Ry", r.reactions[f.dofs_of(1)[1]], 0.0, tol=1e-9)
ef = r.member_results[0].end_forces_local
check("桿端彎矩M1", ef[2], 0.0, tol=1e-9)
check("桿端彎矩M2", ef[5], 0.0, tol=1e-9)
print("PASS: 靜定結構沉陷不產生任何反力/內力(純剛體轉動), 完全符合結構學基本性質\n")


# ---- 案例B: 固接+滾支承梁(一次靜不定)支承沉陷, 對照傾角變位法經典公式 ----
print("=== 案例B: 一次靜不定梁支承沉陷, 對照傾角變位法 M_A=-3EIΔ/L² ===")
L = 8.0
delta = -0.015   # 滾支承端沉陷1.5cm
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's')
f2.fix(0)                    # A端固接
f2.support(1, uy=delta)      # B端滾支承, 沉陷delta (ux仍自由, 不設定=None)
r2 = solve(f2)

# 傾角變位法經典公式(見docstring): M_AB = -3*EI*delta/L^2
M_AB_exact = -3 * E * I * delta / L**2
ef2 = r2.member_results[0].end_forces_local
M_AB_fem = ef2[2]
check("M_AB(固接端彎矩)", M_AB_fem, M_AB_exact, tol=1e-6)

# B端(滾支承, 只擋uy)彎矩應為0(沒有東西可以在那裡產生力矩)
check("M_BA(滾支承端彎矩)", ef2[5], 0.0, tol=1e-9)

# 交叉確認: 反力平衡(沒有外加載重, 只有沉陷, 反力總和應為0)
Ry_A = r2.reactions[f2.dofs_of(0)[1]]
Ry_B = r2.reactions[f2.dofs_of(1)[1]]
check("垂直反力平衡(RyA+RyB=0)", Ry_A + Ry_B, 0.0, tol=1e-9)
print("PASS: 一次靜不定梁沉陷分析吻合傾角變位法經典公式\n")

print("PASS: Phase 3(支承強制位移/沉陷分析) 全數驗證通過")
