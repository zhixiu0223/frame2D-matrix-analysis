"""
驗證案例: frame2d.modal.eigen() -- 模態分析 (動力分析 D2), 只依賴 numpy

層1 解析解(頻率)
  (a) 頂端質量懸臂柱: ω² = 3EI/(mL³), 質量歸一化後 Γ=√m、有效質量=m、累積質量比=1
  (b) 兩質量懸臂(無質量梁)的精確柔度解: 頻率與參與係數; 以及軸向極剛(A預設1e8)+兩方向
      質量、ω²跨度 4e12 的梯度式結構仍精確到 1e-16 (對角縮放的價值)
  (c) 均質懸臂梁彎曲頻率 (β_nL)²√(EI/ρAL⁴): 一致質量從上方 O(h⁴) 收斂、
      集中質量從下方 O(h²) 收斂
  (d) 簡支梁 (nπ)²√(EI/ρAL⁴)
層2 模態性質解析解
  懸臂梁第1模態有效質量比 (∫φ dx)²/(L∫φ² dx) = 0.61308 (φ為精確模態形狀)
層3 結構性質(不依賴外部參考)
  正交性 φᵀMφ = I、φᵀKφ = diag(ω²); 取全部模態時有效質量完備 ΣM* = 可動質量;
  n_modes 截斷、正負號慣例
層4 獨立路徑交叉驗證
  (a) 靜力凝縮: 集中質量(轉角無質量, 走凝縮) vs 給每個轉角加極小轉動慣量ε(質量正定,
      不需凝縮), ε→0 必須收斂到同一組頻率
  (b) release端專屬DOF: 「固定+兩端release」與「鉸支承+滾支承」是同一個物理結構,
      特徵值必須完全相同(集中與一致質量都要)
  (c) 單位一致性: 同一結構用 SI(Pa, kg) 與 kN·m·ton 算, ω 必須一致
層5 明確拒絕: 支承不足(剛體模態)、機構、沒有質量、質量全在支承上
"""
import numpy as np

from frame2d import Frame2D
from frame2d.assembly import assemble_K
from frame2d.mass import assemble_M
from frame2d.modal import eigen

E, I, A, RHO, L = 200e6, 8e-5, 1e-2, 7.85, 3.0
BETA_L = [1.8751040687, 4.6940911330, 7.8547574382]     # 懸臂梁 cosh·cos = -1 的根


def rel(a, b):
    return np.abs(np.asarray(a) / np.asarray(b) - 1.0)


def cantilever(n, kind_ux_restrained=True, E_=E, I_=I, A_=A, rho=RHO):
    f = Frame2D()
    for i in range(n + 1):
        f.add_node(i, L * i / n, 0)
    f.add_section('s', E=E_, I=I_, A=A_, rho=rho)
    for i in range(n):
        f.add_member(i, i, i + 1, 's')
    f.fix(0)
    if kind_ux_restrained:                       # 拿掉軸向模態, 只留彎曲
        for i in range(1, n + 1):
            f.support(i, ux=0.0)
    return f


def simply_supported(n):
    f = Frame2D()
    for i in range(n + 1):
        f.add_node(i, L * i / n, 0)
    f.add_section('s', E=E, I=I, A=A, rho=RHO)
    for i in range(n):
        f.add_member(i, i, i + 1, 's')
    f.pin(0)
    f.roller_y(n)
    for i in range(1, n):
        f.support(i, ux=0.0)
    return f


# ------------------------------------------------------------------
print("=== 層1a: 頂端質量懸臂柱 ===")
m_tip = 2.0
f = Frame2D()
f.add_node(0, 0, 0).add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.fix(0)
f.add_mass(1, my=m_tip)
md = eigen(f)
w_exact = np.sqrt(3 * E * I / (m_tip * L**3))
print(f"  ω = {md.omega[0]:.10f}, 解析解 {w_exact:.10f}, Γy = {md.gamma['y'][0]:.10f} (√m = {np.sqrt(m_tip):.10f})")
assert md.n_modes == 1 and rel(md.omega[0], w_exact) < 1e-12
assert abs(abs(md.gamma['y'][0]) - np.sqrt(m_tip)) < 1e-12
assert abs(md.eff_mass['y'][0] - m_tip) < 1e-12 and abs(md.cum_ratio['y'][0] - 1.0) < 1e-12
assert abs(md.period[0] - 2 * np.pi / w_exact) < 1e-12 and abs(md.frequency[0] - w_exact / (2 * np.pi)) < 1e-12

# ------------------------------------------------------------------
print("=== 層1b: 兩質量懸臂(無質量梁, 精確解析) ===")
# 無質量懸臂梁, 質量 m1 在 L/2、m2 在 L(只考慮橫向)。載重作用在節點時Hermite三次
# 形函數給出精確節點撓度, 所以有限元素結果應該是「精確」的(不是近似)。
m1, m2 = 3.0, 2.0
f = Frame2D()
f.add_node(0, 0, 0).add_node(1, L / 2, 0).add_node(2, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's').add_member(1, 1, 2, 's')
f.fix(0)
f.support(1, ux=0.0).support(2, ux=0.0)
f.add_mass(1, my=m1).add_mass(2, my=m2)
md = eigen(f)
# 獨立手算: 柔度矩陣(梁理論, 位置 a=L/2 與 L), 完全不看frame2d
Fm = (L**3 / (E * I)) * np.array([[1 / 24, 5 / 48], [5 / 48, 1 / 3]])
Mm = np.diag([m1, m2])
lam_ana, V = np.linalg.eig(Fm @ Mm)                   # F M φ = (1/ω²) φ
order = np.argsort(-lam_ana.real)
w_ana = 1.0 / np.sqrt(lam_ana.real[order])
V = V.real[:, order]
Vn = V / np.sqrt(np.einsum('ij,ik,kj->j', V, Mm, V))  # 對M正規化
gam_ana = np.abs(Vn.T @ Mm @ np.ones(2))
print(f"  ω 模型 {md.omega} 解析 {w_ana}")
print(f"  |Γy| 模型 {np.abs(md.gamma['y'])} 解析 {gam_ana}")
assert md.n_modes == 2 and np.all(rel(md.omega, w_ana) < 1e-11), "兩質量懸臂頻率不符精確解"
assert np.allclose(np.abs(md.gamma['y']), gam_ana, rtol=1e-10)
assert abs(md.eff_mass['y'].sum() - (m1 + m2)) < 1e-10 and md.modes_needed(0.9, 'y') in (1, 2)
assert md.modes_needed(2.0, 'y') is None

print("=== 層1b2: 軸向極剛(A預設1e8) + 兩方向質量: 跨度十幾個數量級仍然精確 ===")
f = Frame2D()
f.add_node(0, 0, 0).add_node(1, L, 0)
f.add_section('s', E=E, I=I)                            # A預設1e8, 沒有分佈質量
f.add_member(0, 0, 1, 's')
f.fix(0)
f.add_mass(1, mx=m_tip, my=m_tip)
md = eigen(f)
w_bend = np.sqrt(3 * E * I / (m_tip * L**3))
w_axial = np.sqrt(E * 1e8 / L / m_tip)
print(f"  ω² 跨度 {md.omega[-1]**2 / md.omega[0]**2:.1e}; 彎曲 ω 相對誤差 {rel(md.omega[0], w_bend):.2e}, 軸向 {rel(md.omega[1], w_axial):.2e}")
assert rel(md.omega[0], w_bend) < 1e-9 and rel(md.omega[1], w_axial) < 1e-9, "梯度式(軸向極剛)結構精度不足"

# ------------------------------------------------------------------
print("=== 層1c: 均質懸臂梁彎曲頻率與收斂階數 ===")
exact = np.array([b**2 * np.sqrt(E * I / (RHO * A * L**4)) for b in BETA_L])
errs = {}
for kind in ('lumped', 'consistent'):
    errs[kind] = {}
    for n in (4, 8, 16, 32):
        md = eigen(cantilever(n), n_modes=3, mass=kind)
        errs[kind][n] = md.omega / exact - 1.0
        print(f"  {kind:>10} n={n:<2}: 相對誤差 {errs[kind][n]}")
for n in (4, 8, 16, 32):
    assert np.all(errs['consistent'][n] > 0), "一致質量應從上方收斂(頻率偏高)"
    assert np.all(errs['lumped'][n] < 0), "集中質量應從下方收斂(頻率偏低)"
for n in (4, 8, 16):
    r_c = errs['consistent'][n][0] / errs['consistent'][2 * n][0]
    r_l = errs['lumped'][n][0] / errs['lumped'][2 * n][0]
    print(f"  第1模態誤差比 (n={n}→{2 * n}): 一致 {r_c:.2f} (理論 16), 集中 {r_l:.2f} (理論 4)")
    assert 12 < r_c < 20, f"一致質量收斂階數不是 O(h⁴): 比值 {r_c}"
    assert 3.5 < r_l < 4.5, f"集中質量收斂階數不是 O(h²): 比值 {r_l}"
assert errs['consistent'][16][0] < 2e-7 and errs['consistent'][32][0] < 1e-8

print("=== 層1d: 簡支梁 ===")
exact_ss = np.array([(k * np.pi)**2 * np.sqrt(E * I / (RHO * A * L**4)) for k in (1, 2, 3)])
for kind in ('lumped', 'consistent'):
    md = eigen(simply_supported(16), n_modes=3, mass=kind)
    e = md.omega / exact_ss - 1.0
    print(f"  {kind:>10}: 相對誤差 {e}")
    assert np.all(np.abs(e) < 2e-4) and np.abs(e[0]) < 2e-6
    assert (np.all(e > 0) if kind == 'consistent' else np.all(e < 0))

# ------------------------------------------------------------------
print("=== 層2: 懸臂梁第1模態有效質量比 ===")
b = BETA_L[0]
sig = (np.cosh(b) + np.cos(b)) / (np.sinh(b) + np.sin(b))
xs = np.linspace(0, 1, 400001)
phi_ex = np.cosh(b * xs) - np.cos(b * xs) - sig * (np.sinh(b * xs) - np.sin(b * xs))
trapz = getattr(np, 'trapezoid', None) or np.trapz
ratio_exact = trapz(phi_ex, xs)**2 / trapz(phi_ex**2, xs)
md = eigen(cantilever(16), n_modes=1, mass='consistent')
ratio_fe = md.eff_mass['y'][0] / md.mass_total['y']
print(f"  精確 {ratio_exact:.8f} (文獻約0.613), 有限元素(n=16) {ratio_fe:.8f}, 相對差 {abs(ratio_fe / ratio_exact - 1):.2e}")
assert abs(ratio_exact - 0.6131) < 1e-3
assert abs(ratio_fe / ratio_exact - 1) < 1e-6

# ------------------------------------------------------------------
def mixed_frame():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4.5).add_node(3, 6.5, 0).add_node(4, 3, 7.5)
    f.add_section('c', E=E, I=I, A=A, rho=RHO)
    f.add_section('b', E=E, I=2 * I, A=2 * A, rho=RHO)
    f.add_section('t', E=E, I=I, A=1e-3, rho=RHO)
    f.add_member(0, 0, 1, 'c').add_member(1, 1, 2, 'b').add_member(2, 3, 2, 'c').add_member(3, 1, 4, 'b')
    f.add_truss(4, 4, 2, 't')
    f.fix(0).fix(3)
    f.add_mass(1, mx=3.0, my=1.0, Iz=0.5).add_mass(2, mx=2.0, my=2.5).add_mass(4, mx=1.5, my=1.0, Iz=0.1)
    return f


print("=== 層3: 正交性、完備性、截斷、正負號 ===")
for kind in ('lumped', 'consistent'):
    f = mixed_frame()
    md = eigen(f, mass=kind)
    K = assemble_K(f).K
    M = assemble_M(f, kind)
    P = md.phi
    orth_m = np.max(np.abs(P.T @ M @ P - np.eye(md.n_modes)))
    orth_k = np.max(np.abs(P.T @ K @ P - np.diag(md.omega**2))) / np.max(md.omega**2)
    print(f"  [{kind}] 模態數 {md.n_modes}: max|φᵀMφ-I| = {orth_m:.2e}, max|φᵀKφ-diag(ω²)|/ω²max = {orth_k:.2e}")
    assert orth_m < 1e-9 and orth_k < 1e-9
    # 獨立算「可動質量」: 總質量 - 直接掛在支承DOF上的質量 (Schur補: M_FF - M_FD M_DD⁻¹ M_DF)
    fixed = np.zeros(M.shape[0], dtype=bool)
    for sp_ in f.supports:
        for dof, val in zip(f.dofs_of(sp_.node), (sp_.ux, sp_.uy, sp_.rot)):
            if val is not None:
                fixed[dof] = True
    D = np.where((~fixed) & (np.abs(M).max(axis=1) > 0))[0]
    F = np.where(fixed)[0]
    Mtilde = M[np.ix_(F, F)] - M[np.ix_(F, D)] @ np.linalg.solve(M[np.ix_(D, D)], M[np.ix_(D, F)])
    for d in ('x', 'y'):
        r = np.zeros(M.shape[0])
        for nid in f.nodes:
            r[f.dofs_of(nid)[0 if d == 'x' else 1]] = 1.0
        base_mass = float(r[F] @ Mtilde @ r[F])
        want_free = md.mass_total[d] - base_mass
        s = md.eff_mass[d].sum()
        print(f"  [{kind}] 方向{d}: ΣM* = {s:.10f}, 可動質量 = {md.mass_free[d]:.10f}, "
              f"總質量-支承質量 = {want_free:.10f}, 累積比末項 = {md.cum_ratio[d][-1]:.10f}")
        assert abs(s - md.mass_free[d]) < 1e-9 * md.mass_free[d], f"{kind}/{d}: 取全部模態時有效質量必須完備"
        assert abs(md.mass_free[d] - want_free) < 1e-9 * want_free, f"{kind}/{d}: 可動質量與獨立Schur補算法不一致"
        assert abs(md.cum_ratio[d][-1] - 1.0) < 1e-9
        assert md.mass_free[d] <= md.mass_total[d] + 1e-12
        assert abs(md.cum_ratio_total[d][-1] - md.mass_free[d] / md.mass_total[d]) < 1e-12
        if kind == 'lumped':                                    # 集中質量: 扣掉的正好是支承DOF上的節點質量
            fixed_mass = sum(M[f.dofs_of(sp_.node)[0 if d == 'x' else 1], f.dofs_of(sp_.node)[0 if d == 'x' else 1]]
                             for sp_ in f.supports)
            assert abs(md.mass_total[d] - md.mass_free[d] - fixed_mass) < 1e-9
    assert np.all(np.diff(md.omega) >= -1e-12), "頻率必須由小到大"
    # 截斷
    md3 = eigen(f, n_modes=3, mass=kind)
    assert md3.n_modes == 3 and np.allclose(md3.omega, md.omega[:3], rtol=1e-12)
    # 正負號慣例與可重現性(一致質量下所有自由DOF都是動力DOF, 全向量最大分量必為正)
    md_again = eigen(f, mass=kind)
    assert np.array_equal(md.phi, md_again.phi), "模態形狀必須可重現(含正負號)"
    if kind == 'consistent':
        for j in range(md.n_modes):
            assert md.phi[np.argmax(np.abs(md.phi[:, j])), j] > 0, "正負號慣例: 絕對值最大分量為正"
    txt = md.table()
    assert 'T (s)' in txt and len(txt.splitlines()) == md.n_modes + 1
print("  模態表輸出:"); print("\n".join("    " + l for l in eigen(mixed_frame(), n_modes=3).table().splitlines()))

# ------------------------------------------------------------------
print("=== 層4a: 靜力凝縮 vs 極小轉動慣量(不需凝縮的獨立路徑) ===")
f0 = mixed_frame()
md0 = eigen(f0, n_modes=6, mass='lumped')                    # 走凝縮
for eps in (1e-4, 1e-6, 1e-8):
    f1 = mixed_frame()
    for nid in f1.nodes:
        f1.add_mass(nid, Iz=eps)                             # 每個轉角加極小質量 -> M正定, 不凝縮
    md1 = eigen(f1, n_modes=6, mass='lumped')
    e = np.max(rel(md1.omega, md0.omega))
    print(f"  ε = {eps:.0e}: 最大頻率相對差 = {e:.2e}")
assert e < 1e-6, "ε→0 時應收斂到靜力凝縮的結果"

print("=== 層4b: release端專屬DOF ===")
def chain(n, mode):
    f = Frame2D()
    for i in range(n + 1):
        f.add_node(i, L * i / n, 0)
    f.add_section('s', E=E, I=I, A=A, rho=RHO)
    for i in range(n):
        f.add_member(i, i, i + 1, 's')
    for i in range(1, n):
        f.support(i, ux=0.0)
    if mode == 'pinned':
        f.pin(0)
        f.roller_y(n)
    else:                                                    # 固定 + 兩端release: 物理上還是簡支
        f.fix(0)
        f.support(n, ux=None, uy=0.0, rot=0.0)
        f.members[0].release_i = True
        f.members[n - 1].release_j = True
    return f
for kind in ('lumped', 'consistent'):
    a = eigen(chain(8, 'pinned'), n_modes=5, mass=kind)
    bb = eigen(chain(8, 'released'), n_modes=5, mass=kind)
    e = np.max(rel(bb.omega, a.omega))
    print(f"  [{kind}] 鉸支承 vs 固定+release: 最大頻率相對差 = {e:.2e}")
    assert e < 1e-10

print("=== 層4c: 單位一致性 ===")
def frame_in_units(E_, rho_, mass_scale, force_scale):
    f = mixed_frame()
    for s in f.sections.values():
        s.E = E_
        s.rho = rho_
    for nm in f.node_masses:
        nm.mx *= mass_scale; nm.my *= mass_scale; nm.Iz *= mass_scale
    return f
a = eigen(frame_in_units(E, RHO, 1.0, 1.0), n_modes=5)                # kN, m, ton
b_ = eigen(frame_in_units(E * 1e3, RHO * 1e3, 1e3, 1e3), n_modes=5)   # N, m, kg
print(f"  kN·m·ton vs SI: 最大 ω 相對差 = {np.max(rel(a.omega, b_.omega)):.2e}")
assert np.max(rel(a.omega, b_.omega)) < 1e-12

# ------------------------------------------------------------------
print("=== 層5: 明確拒絕 ===")
def expect_error(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該raise ValueError")

f = Frame2D(); f.add_node(0, 0, 0).add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A, rho=RHO); f.add_member(0, 0, 1, 's')
expect_error("完全沒有支承(剛體模態)", lambda: eigen(f))
f = Frame2D(); f.add_node(0, 0, 0).add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A, rho=RHO); f.add_member(0, 0, 1, 's'); f.pin(0)
expect_error("鉸支承懸臂(轉動機構)", lambda: eigen(f))
f = Frame2D(); f.add_node(0, 0, 0).add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A); f.add_member(0, 0, 1, 's'); f.fix(0)
expect_error("沒有任何質量", lambda: eigen(f))
f.add_mass(0, mx=5.0, my=5.0)                                  # 質量只在固定端
expect_error("質量全部在支承上", lambda: eigen(f))
expect_error("非法mass種類", lambda: eigen(mixed_frame(), mass='bogus'))
expect_error("n_modes<1", lambda: eigen(mixed_frame(), n_modes=0))

print("\n全部通過: 頻率、模態性質、正交性、完備性、收斂階數、凝縮、release、單位一致性都吻合。")
