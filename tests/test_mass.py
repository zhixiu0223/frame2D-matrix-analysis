"""
驗證案例: frame2d.mass -- 質量矩陣 (動力分析 D1)

層1 獨立推導: 局部6x6一致質量矩陣, 用sympy從Hermite形函數直接積分
    m_ij = ∫ m̄ N_i N_j dx, 不看mass.py裡寫死的156/22L/54/-13L係數。
層2 不依賴係數的解析不變量:
    (a) 總質量: rᵀ M r = 總質量 (集中/一致, 含斜桿、桁架、節點質量、release)
    (b) 繞原點剛體轉動的極慣量: 一致質量必須精確等於
        Σ ρAL(|c|² + L²/12)  (c = 桿件形心到原點的向量, 均勻桿的極慣量),
        集中質量則是 Σ ρAL(|c|² + L²/4)
        —— 剛體轉動位移場(線性)剛好在形函數空間內, 所以一致質量是精確的,
        不是近似。這一條同時驗證了旋轉矩陣T的用法。
層3 結構性質: 對稱、半正定、集中質量是對角、DOF編號跟assemble_K一致。
層4 跟D0的K接起來: 頂端質量懸臂柱, 靜力凝縮無質量DOF後 ω² = 3EI/(mL³);
    同一個物理結構分別用SI(Pa,N,kg)與kN-m-ton兩套單位算, ω必須完全相同
    (單位一致性檢核)。
層5 明確拒絕: cable / equal_dof / 非零指定位移 / 非法參數。
"""
import numpy as np
import sympy as sp

from frame2d import Frame2D
from frame2d.assembly import assemble_K
from frame2d.mass import (assemble_M, consistent_mass_local, influence_vector,
                          total_mass)

E, I, A = 200e6, 8e-5, 1e-2


def rel(a, b):
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


# ------------------------------------------------------------------
# 層1: sympy獨立推導局部一致質量矩陣
# ------------------------------------------------------------------
print("=== 層1: sympy從形函數積分推導局部一致質量矩陣 ===")
x, Ls, mbar = sp.symbols('x L mbar', positive=True)
xi = x / Ls
N_bend = [1 - 3 * xi**2 + 2 * xi**3,
          Ls * (xi - 2 * xi**2 + xi**3),
          3 * xi**2 - 2 * xi**3,
          Ls * (-xi**2 + xi**3)]
N_axial = [1 - xi, xi]

L_val, mbar_val = 3.7, 2.5
m_sym = sp.zeros(6, 6)
bend_idx, axial_idx = [1, 2, 4, 5], [0, 3]
for a, ia in enumerate(bend_idx):
    for b, ib in enumerate(bend_idx):
        m_sym[ia, ib] = sp.integrate(mbar * N_bend[a] * N_bend[b], (x, 0, Ls))
for a, ia in enumerate(axial_idx):
    for b, ib in enumerate(axial_idx):
        m_sym[ia, ib] = sp.integrate(mbar * N_axial[a] * N_axial[b], (x, 0, Ls))
m_ref = np.array(m_sym.subs({Ls: L_val, mbar: mbar_val}).evalf().tolist(), dtype=float)
m_code = consistent_mass_local(mbar_val, L_val)
d = rel(m_code, m_ref)
print(f"  局部6x6一致質量矩陣 vs sympy積分: 最大相對差 = {d:.2e}")
assert d < 1e-13, "一致質量矩陣係數跟形函數積分不一致"


# ------------------------------------------------------------------
# 層2/3: 不變量與結構性質(一個含斜桿、桁架、節點質量的模型)
# ------------------------------------------------------------------
def mixed_model(release=False):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0).add_node(4, 3, 7)
    f.add_section('beam', E=E, I=I, A=A, rho=7.85)
    f.add_section('bar', E=E, I=I, A=2e-3, rho=7.85)
    f.add_member(0, 0, 1, 'beam').add_member(1, 1, 2, 'beam').add_member(2, 3, 2, 'beam')
    f.add_member(3, 1, 4, 'beam')                 # 斜桿(frame)
    f.add_truss(4, 4, 2, 'bar')                   # 斜桿(truss)
    if release:
        f.members[1].release_j = True
    f.fix(0).fix(3)
    f.add_mass(1, mx=3.0, my=1.0, Iz=0.5)
    f.add_mass(4, mx=2.0, my=2.5, Iz=0.2)
    f.add_mass(4, mx=1.0)                          # 同節點再加一次, 要加總
    return f


def rigid_rotation_vector(f, n_total):
    r = np.zeros(n_total)
    for nid, node in f.nodes.items():
        ux, uy, rot = f.dofs_of(nid)
        r[ux], r[uy], r[rot] = -node.y, node.x, 1.0
    return r


def expected_polar(f, kind):
    tot = 0.0
    for m in f.members.values():
        sec = f.sections[m.section]
        pi, pj = f.nodes[m.node_i], f.nodes[m.node_j]
        L = np.hypot(pj.x - pi.x, pj.y - pi.y)
        cx, cy = (pi.x + pj.x) / 2, (pi.y + pj.y) / 2
        c2 = cx * cx + cy * cy
        extra = L * L / 12 if kind == 'consistent' else L * L / 4
        tot += sec.rho * sec.A * L * (c2 + extra)
    for nm in f.node_masses:
        n = f.nodes[nm.node]
        tot += nm.mx * n.y**2 + nm.my * n.x**2 + nm.Iz
    return tot


print("=== 層2/3: 總質量與極慣量不變量、結構性質 ===")
f = mixed_model()
K_asm = assemble_K(f)
for kind in ('lumped', 'consistent'):
    M = assemble_M(f, kind)
    assert M.shape == K_asm.K.shape, "M與K的維度(DOF編號)不一致"
    for direction in ('x', 'y'):
        r = influence_vector(f, direction)
        got, want = float(r @ M @ r), total_mass(f, direction)
        print(f"  [{kind}] 方向{direction}: rᵀMr = {got:.10f}, 總質量 = {want:.10f}")
        assert abs(got - want) < 1e-10 * want, f"{kind}/{direction}: 總質量不守恆"
    rr = rigid_rotation_vector(f, M.shape[0])
    got, want = float(rr @ M @ rr), expected_polar(f, kind)
    print(f"  [{kind}] 繞原點剛體轉動: rᵀMr = {got:.10f}, 解析極慣量 = {want:.10f}")
    assert abs(got - want) < 1e-10 * want, f"{kind}: 極慣量不符解析解"
    assert np.allclose(M, M.T, rtol=0, atol=1e-12 * np.max(np.abs(M))), f"{kind}: M不對稱"
    ev = np.linalg.eigvalsh(M)
    assert ev[0] > -1e-12 * ev[-1], f"{kind}: M不是半正定"
    if kind == 'lumped':
        assert np.count_nonzero(M - np.diag(np.diag(M))) == 0, "集中質量矩陣必須是對角的"

# 節點質量落在正確的DOF、且同節點多次add_mass會加總
Ml = assemble_M(mixed_model(), 'lumped')
f4 = mixed_model()
ux4, uy4, rot4 = f4.dofs_of(4)
tributary = 7.85 * 2e-3 * np.hypot(3, 3) / 2 + 7.85 * 1e-2 * np.hypot(3, 3) / 2   # truss(4-2)與frame(1-4)各分一半
assert abs(Ml[ux4, ux4] - (3.0 + tributary)) < 1e-12, "節點4的ux質量 = 節點質量(2+1) + 兩根桿件各一半"
assert abs(Ml[rot4, rot4] - 0.2) < 1e-15, "節點4的轉角質量只來自Iz"
print("  節點質量加總與DOF位置正確")

print("=== 層2b: release端(額外DOF) ===")
fr = mixed_model(release=True)
Kr = assemble_K(fr)
assert Kr.n_extra_dof == 1
for kind in ('lumped', 'consistent'):
    M = assemble_M(fr, kind)
    assert M.shape == Kr.K.shape
    r = influence_vector(fr, 'x')
    assert abs(r @ M @ r - total_mass(fr, 'x')) < 1e-10 * total_mass(fr, 'x')
    extra = Kr.n_node_dof
    row = np.abs(M[extra, :]).max()
    print(f"  [{kind}] 額外DOF那一列最大質量項 = {row:.3e}")
    if kind == 'lumped':
        assert row == 0.0, "集中質量下release專屬轉角DOF必須無質量(D2要靜力凝縮)"
    else:
        assert row > 0.0, "一致質量下release專屬DOF有自己的慣性項"

# ------------------------------------------------------------------
# 層4: 跟D0的K接起來 -- 頂端質量懸臂柱 + 單位一致性
# ------------------------------------------------------------------
def tip_mass_cantilever(E_, A_, I_, m_tip, L_):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, L_, 0)
    f.add_section('s', E=E_, I=I_, A=A_)          # 桿件本身沒有質量(rho=None)
    f.add_member(0, 0, 1, 's')
    f.fix(0)
    f.add_mass(1, my=m_tip)
    return f


def sdof_omega(f):
    """靜力凝縮無質量DOF後的頻率。動力DOF = M對角>0的自由DOF。"""
    asm = assemble_K(f)
    M = assemble_M(f, 'lumped')
    fixed = np.zeros(asm.n_dof, dtype=bool)
    for s in f.supports:
        for dof, val in zip(f.dofs_of(s.node), (s.ux, s.uy, s.rot)):
            if val is not None:
                fixed[dof] = True
    free = np.where(~fixed)[0]
    dyn = np.array([i for i in free if M[i, i] > 0])
    mless = np.array([i for i in free if M[i, i] == 0])
    Kdd = asm.K[np.ix_(dyn, dyn)]
    Kdm = asm.K[np.ix_(dyn, mless)]
    Kmm = asm.K[np.ix_(mless, mless)]
    Keff = Kdd - Kdm @ np.linalg.solve(Kmm, Kdm.T)
    Md = M[np.ix_(dyn, dyn)]
    return float(np.sqrt(np.linalg.eigvals(np.linalg.solve(Md, Keff)).real.min()))


print("=== 層4: 頂端質量懸臂柱 ω² = 3EI/(mL³), 兩套單位 ===")
L_c = 3.0
# kN, m, ton, s
w_kn = sdof_omega(tip_mass_cantilever(200e6, 1e-2, 8e-5, 2.0, L_c))
# SI: Pa, m, N, kg (同一個物理結構: E 200 GPa, 質量2 ton = 2000 kg)
w_si = sdof_omega(tip_mass_cantilever(200e9, 1e-2, 8e-5, 2000.0, L_c))
w_exact = np.sqrt(3 * 200e6 * 8e-5 / (2.0 * L_c**3))
print(f"  kN-m-ton: ω = {w_kn:.10f} rad/s; SI: ω = {w_si:.10f} rad/s; 解析解 = {w_exact:.10f}")
assert abs(w_kn - w_exact) / w_exact < 1e-10, "kN-m-ton: 跟解析解不一致"
assert abs(w_si - w_exact) / w_exact < 1e-10, "SI: 跟解析解不一致"
assert abs(w_si - w_kn) / w_kn < 1e-10, "兩套單位算同一個結構, ω必須一致"
# 反面檢核: 故意把ton當kg(SI力+kg質量以外的錯誤組合)必須讓ω差√1000倍
w_bad = sdof_omega(tip_mass_cantilever(200e9, 1e-2, 8e-5, 2.0, L_c))
print(f"  (反面) SI力單位配上『2』(當成kg而不是2000kg): ω = {w_bad:.4f}, 是正確值的 {w_bad / w_exact:.3f} 倍 (預期√1000={np.sqrt(1000):.3f})")
assert abs(w_bad / w_exact - np.sqrt(1000)) < 1e-9

# ------------------------------------------------------------------
# 層5: 明確拒絕與參數檢查
# ------------------------------------------------------------------
print("=== 層5: 明確拒絕 ===")


def expect_error(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:38]}...)")
        return
    raise AssertionError(f"{label}: 應該raise ValueError")


f = mixed_model(); f.add_section('c', E=E, I=I, A=1e-3, rho=1.0); f.add_cable(9, 0, 2, 'c')
expect_error("含cable", lambda: assemble_M(f))
f = mixed_model(); f.equal_dof(1, 2, ux=True)
expect_error("含equal_dof", lambda: assemble_M(f))
f = mixed_model(); f.support(2, uy=-0.01)
expect_error("非零指定位移", lambda: assemble_M(f))
expect_error("非法kind", lambda: assemble_M(mixed_model(), 'diagonal'))
expect_error("負密度", lambda: Frame2D().add_section('s', 1, 1, 1, rho=-1.0))
expect_error("負質量", lambda: Frame2D().add_mass(0, mx=-1.0))
f = mixed_model(); f.add_mass(99, mx=1.0)
expect_error("質量指定到不存在的節點", lambda: assemble_M(f))
expect_error("非法方向", lambda: influence_vector(mixed_model(), 'z'))

# rho=None: 桿件沒有分佈質量, 只剩節點質量
f = Frame2D(); f.add_node(0, 0, 0).add_node(1, 2, 0)
f.add_section('s', E=E, I=I, A=A); f.add_member(0, 0, 1, 's'); f.fix(0); f.add_mass(1, mx=5.0)
M0 = assemble_M(f, 'consistent')
assert np.count_nonzero(M0) == 1 and abs(M0[f.dofs_of(1)[0], f.dofs_of(1)[0]] - 5.0) < 1e-15
print("  rho=None: 桿件不貢獻質量, 只剩節點質量 OK")

print("\n全部通過: 質量矩陣與sympy積分、總質量、極慣量、SDOF解析解、單位一致性都吻合。")
