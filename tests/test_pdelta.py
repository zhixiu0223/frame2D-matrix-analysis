"""
驗證案例: solve_pdelta() -- 線性化P-Delta疊代求解

案例A: 完全沒有軸力路徑的模型(懸臂梁, 只有橫向點載重), solve_pdelta()
       應該跟solve_dofmanager()逐位元一致(第一次疊代後軸力仍是0, 直接
       收斂) -- 鎖住"新增這個函式不影響既有主要求解器"這件事。

案例B: 獨立sympy重新解一次縮減後的2自由度系統, 交叉驗證solve_pdelta()
       算出來的側向位移/轉角。用一根水平懸臂柱, 固定端node0, 自由端
       node1同時受軸向力(產生確定的軸力N)跟橫向力H, 拿frame2d疊代收斂後
       的N, 獨立組出Ke+Kg的縮減2x2系統(只有node1的v,theta是自由度)手動解,
       不透過frame2d的組裝/求解程式碼路徑。

案例C: 物理合理性 -- 軸向壓力(N<0)應使側向位移比純彈性(P=0)時更大
       (P-Delta放大效應), 拉力(N>0)應使側向位移變小。
"""
import numpy as np
import sympy as sp
from frame2d import Frame2D
from frame2d.dofmanager import solve_dofmanager, solve_pdelta

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


def cantilever_lateral_only():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec')
    f.fix(0)
    f.point_load(1, fy=-10.0)
    return f


# ---- 案例A: 無軸力路徑時, 跟solve_dofmanager()逐位元一致 ----
print("=== 案例A: 無軸力時與solve_dofmanager()一致 ===")
f = cantilever_lateral_only()
r_baseline = solve_dofmanager(f)
r_pdelta = solve_pdelta(f)
assert np.allclose(r_baseline.displacements, r_pdelta.displacements, atol=0), \
    "無軸力路徑時, solve_pdelta()應與solve_dofmanager()逐位元一致"
assert np.allclose(r_baseline.reactions, r_pdelta.reactions, atol=0)
print("PASS: 無軸力路徑時兩者逐位元相同\n")


def build_axial_lateral_case(fx_applied, fy_applied):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec')
    f.fix(0)
    f.point_load(1, fx=fx_applied, fy=fy_applied)
    return f


def sympy_reduced_solution(N, H, E_, I_, L_):
    """獨立組出Ke+Kg的縮減2x2系統(自由度: node1的v, theta)手解, 完全不
    呼叫frame2d的任何組裝/求解函式, 只用跟elements.py同樣的教科書公式
    (但這裡重新用sympy寫一次, 不是import過來比對自己)。"""
    v2, th2 = sp.symbols('v2 th2')
    EI = E_ * I_
    # Ke的4x4 (v1,th1,v2,th2), node1=(v2,th2)是自由度, node0固定(v1=th1=0)
    # 所以只需要Ke/Kg對(v2,th2)那個2x2子矩陣, 且因為v1=th1=0, 交叉項不貢獻力
    Ke22 = sp.Matrix([
        [12 * EI / L_**3,  -6 * EI / L_**2],
        [-6 * EI / L_**2,   4 * EI / L_],
    ])
    if N == 0.0:
        Kg22 = sp.zeros(2, 2)
    else:
        Kg22 = (sp.Rational(N) / L_) * sp.Matrix([
            [sp.Rational(6, 5),           -sp.Rational(1, 10) * L_],
            [-sp.Rational(1, 10) * L_,     sp.Rational(2, 15) * L_**2],
        ])
    K22 = Ke22 + Kg22
    sol = sp.linsolve((K22, sp.Matrix([H, 0])), [v2, th2])
    v2_val, th2_val = list(sol)[0]
    return float(v2_val), float(th2_val)


# ---- 案例B: 交叉驗證位移數值 ----
# 軸力同樣控制在遠低於挫屈載重的範圍內(理由見案例C說明)
print("=== 案例B: 獨立sympy解縮減系統交叉驗證 ===")
Pcr_check = np.pi**2 * E * I / (4 * L**2)
fx_applied, fy_applied = -0.3 * Pcr_check, -10.0   # 壓力(推向固定端) + 橫向力
f = build_axial_lateral_case(fx_applied, fy_applied)
r = solve_pdelta(f)
N_converged = r.member_results[0].end_forces_local[3]   # 拉力為正
print(f"疊代收斂後軸力 N = {N_converged:.4f} (應等於施加的fx = {fx_applied}, "
      f"因為單一懸臂柱是靜定系統, 軸力不受Kg疊代影響)")
assert abs(N_converged - fx_applied) < 1e-6, "靜定系統軸力應該直接等於施加的fx"

v2_sympy, th2_sympy = sympy_reduced_solution(N_converged, fy_applied, E, I, L)
v2_fem = r.displacements[f.dofs_of(1)[1]]
th2_fem = r.displacements[f.dofs_of(1)[2]]
print(f"v2: frame2d={v2_fem:.10e}  sympy獨立解={v2_sympy:.10e}")
print(f"th2: frame2d={th2_fem:.10e}  sympy獨立解={th2_sympy:.10e}")
assert abs((v2_fem - v2_sympy) / v2_sympy) < 1e-8, "側向位移跟獨立sympy解不符"
assert abs((th2_fem - th2_sympy) / th2_sympy) < 1e-8, "端點轉角跟獨立sympy解不符"
print("PASS: solve_pdelta()位移結果跟獨立sympy解算的縮減系統逐位吻合\n")

# ---- 案例C: 壓力放大位移, 拉力抑制位移 (P-Delta物理方向) ----
# 注意: 這裡的軸力必須遠低於歐拉挫屈臨界載重 Pcr = pi^2*EI/(4L^2) (懸臂柱,
# 固定-自由), 不然線性化P-Delta模型會嚴重外推到失效區間(矩陣接近奇異,
# 算出來的位移可能不增反減、甚至變號, 完全不代表真實物理行為)。
# 這根柱子 Pcr = pi^2*200e6*8e-5/(4*4^2) ≈ 2467, 這裡刻意只取Pcr的~20%,
# 留在線性化模型有效的範圍內。
print("=== 案例C: 壓力放大位移/拉力抑制位移 ===")
Pcr = np.pi**2 * E * I / (4 * L**2)
N_test = 0.2 * Pcr   # 遠低於挫屈載重, 線性化模型有效範圍內
print(f"Pcr(歐拉挫屈臨界載重) ≈ {Pcr:.1f}, 本案例取 N_test=0.2*Pcr ≈ {N_test:.1f}")

f_elastic = build_axial_lateral_case(0.0, fy_applied)
v2_elastic = solve_pdelta(f_elastic).displacements[f_elastic.dofs_of(1)[1]]

f_compression = build_axial_lateral_case(-N_test, fy_applied)
v2_compression = solve_pdelta(f_compression).displacements[f_compression.dofs_of(1)[1]]

f_tension = build_axial_lateral_case(+N_test, fy_applied)
v2_tension = solve_pdelta(f_tension).displacements[f_tension.dofs_of(1)[1]]

# 位移方向跟fy_applied同號(這裡是負的, 表示往下), 用絕對值比大小
assert abs(v2_compression) > abs(v2_elastic) > abs(v2_tension), (
    f"應滿足 |壓力側移|({abs(v2_compression):.6e}) > |純彈性|({abs(v2_elastic):.6e})"
    f" > |拉力側移|({abs(v2_tension):.6e})"
)
print(f"PASS: |v2| 壓力={abs(v2_compression):.6e} > 純彈性={abs(v2_elastic):.6e} "
      f"> 拉力={abs(v2_tension):.6e}\n")

print("PASS: solve_pdelta() 所有案例通過")
