"""
驗證案例: frame2d.spectrum -- 反應譜分析 (動力分析 D3), 只依賴 numpy

層1 SDOF 解析解: 頂端質量懸臂柱(單一模態)。常數譜 S_a≡A0 時:
    位移 = S_d = A0/ω², 基底剪力 = m·A0(標準 SDOF 公式), 反力方向與量值(R=Ku-F 慣例,
    跟 SolveResult.reactions 同一套, 見 dofmanager.solve_dofmanager 的 `R = K @ u - F`)。
層2 兩質量懸臂(無質量梁, 精確解析, 沿用 test_modal.py 的獨立柔度矩陣參考): 每個模態的
    訊號位移、SRSS/CQC 組合位移, 跟手算的柔度矩陣解對照。
層3 反力/基底剪力的簽名恆等式(不依賴 combine, 對每個模態分別成立, 是全域動量平衡的結果):
    Σ_dof(modal_reactions[dof][i]·r[dof]) == -modal_base_shear[i]   (精確相等)
    在多支承、斜桿、桁架混合的門型剛架上驗證, 比 SDOF 更有說服力。
層4 獨立 K 解交叉驗證: 對其中一個模態, 直接解 K u = f_i*(自由DOF線性方程式, 不用
    modal.py/spectrum.py 共用的任何捷徑), 確認等於 modal_disp[i](驗證「D_i 不用解方程式」
    這個捷徑本身是對的, 不是自我一致的巧合)。
層5 CQC 性質: ρ 對角=1、對稱; 頻率分離良好時 CQC → SRSS。
層6 n_modes 截斷、明確拒絕(方向/combine/damping/n_modes/spectrum回傳負值或nan)。
層7 對 OpenSeesPy `responseSpectrumAnalysis`(選用, 沒裝就 SKIPPED): 單一模態(集中質量、
    一致質量)的位移與反力, 到機器精度。**API 的多模態 SRSS/CQC 組合語法不明確**(嘗試多種
    `-mode` 傳法都沒有得到看起來正確的組合結果, 見 test_spectrum_vs_openseespy.py 開頭說明),
    所以多模態組合只用層2/3/4 的獨立解析/數值交叉驗證, 不對 OpenSeesPy 比對。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.assembly import assemble_K
from frame2d.mass import assemble_M
from frame2d.modal import eigen
from frame2d.spectrum import combine_cqc, combine_srss, cqc_rho, response_spectrum

E, I, A, RHO, L = 200e6, 8e-5, 1e-2, 7.85, 3.0


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: SDOF 解析解(頂端質量懸臂柱) ===")
m_tip = 2.0
f = Frame2D()
f.add_node(0, 0, 0).add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.fix(0)
f.add_mass(1, my=m_tip)
md = eigen(f)
w = md.omega[0]
A0 = 5.0
res = response_spectrum(md, lambda T: A0, direction='y', combine='SRSS')
tip = f.dofs_of(1)[1]
Sd_exact = A0 / w**2
print(f"  位移 = {res.displacements[tip]:.10f}, 解析 S_d = {Sd_exact:.10f}")
assert rel(res.displacements[tip], Sd_exact) < 1e-12
print(f"  基底剪力 = {res.base_shear:.10f}, 解析 m·A0 = {m_tip * A0:.10f}")
assert rel(res.base_shear, m_tip * A0) < 1e-12
# 反力慣例 R = K@D - F: 跟 dofmanager.solve_dofmanager 同一套, 這裡驗證大小相等、方向相反
dof_y0 = f.dofs_of(0)[1]
assert rel(res.modal_reactions[dof_y0][0], -res.modal_base_shear[0]) < 1e-12
print(f"  底部反力(uy) = {res.reactions[dof_y0]:.10f}, |應等於基底剪力| = {res.base_shear:.10f}")
assert rel(res.reactions[dof_y0], res.base_shear) < 1e-12

# ------------------------------------------------------------------
print("=== 層2: 兩質量懸臂(精確解析, 沿用 test_modal.py 的獨立柔度矩陣參考) ===")
m1, m2 = 3.0, 2.0
f2 = Frame2D()
f2.add_node(0, 0, 0).add_node(1, L / 2, 0).add_node(2, L, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's').add_member(1, 1, 2, 's')
f2.fix(0)
f2.support(1, ux=0.0).support(2, ux=0.0)
f2.add_mass(1, my=m1).add_mass(2, my=m2)
md2 = eigen(f2)
Fm = (L**3 / (E * I)) * np.array([[1 / 24, 5 / 48], [5 / 48, 1 / 3]])
Mm = np.diag([m1, m2])
lam_ana, V = np.linalg.eig(Fm @ Mm)
order = np.argsort(-lam_ana.real)
w_ana = 1.0 / np.sqrt(lam_ana.real[order])
V = V.real[:, order]
Vn = V / np.sqrt(np.einsum('ij,ik,kj->j', V, Mm, V))          # 對 M 正規化
gam_ana = Vn.T @ Mm @ np.ones(2)                                # 帶正負號的 Γ

A0 = 7.0
spectrum_const = lambda T: A0
res2 = response_spectrum(md2, spectrum_const, direction='y', combine='SRSS')
dof1, dof2 = f2.dofs_of(1)[1], f2.dofs_of(2)[1]
for i in range(2):
    Sd_i = A0 / w_ana[i]**2
    u1_ana = gam_ana[i] * Sd_i * Vn[0, i]
    u2_ana = gam_ana[i] * Sd_i * Vn[1, i]
    # frame2d 的正負號慣例跟手算的特徵向量可能整體差一個負號, 用比值的絕對值比對星等
    ratio1 = res2.modal_displacements[i, dof1] / u1_ana
    ratio2 = res2.modal_displacements[i, dof2] / u2_ana
    print(f"  模態{i + 1}: ω 相對誤差 {rel(md2.omega[i], w_ana[i]):.2e}, 位移比值(應為±1) {ratio1:.10f} / {ratio2:.10f}")
    assert rel(md2.omega[i], w_ana[i]) < 1e-10
    assert abs(abs(ratio1) - 1) < 1e-9 and abs(abs(ratio2) - 1) < 1e-9 and abs(ratio1 - ratio2) < 1e-9

u1_srss_ana = np.sqrt(sum((gam_ana[i] * (A0 / w_ana[i]**2) * Vn[0, i])**2 for i in range(2)))
u2_srss_ana = np.sqrt(sum((gam_ana[i] * (A0 / w_ana[i]**2) * Vn[1, i])**2 for i in range(2)))
print(f"  SRSS 位移: frame2d {res2.displacements[dof1]:.10f}/{res2.displacements[dof2]:.10f}, "
      f"解析 {u1_srss_ana:.10f}/{u2_srss_ana:.10f}")
assert rel(res2.displacements[dof1], u1_srss_ana) < 1e-10
assert rel(res2.displacements[dof2], u2_srss_ana) < 1e-10

res2c = response_spectrum(md2, spectrum_const, direction='y', combine='CQC', damping=0.05)
rho_ana = cqc_rho(w_ana, 0.05)
u1_cqc_ana = np.sqrt(sum(rho_ana[i, j] * (gam_ana[i] * (A0 / w_ana[i]**2) * Vn[0, i])
                         * (gam_ana[j] * (A0 / w_ana[j]**2) * Vn[0, j]) for i in range(2) for j in range(2)))
print(f"  CQC 位移(dof1): frame2d {res2c.displacements[dof1]:.10f}, 解析 {u1_cqc_ana:.10f}")
assert rel(res2c.displacements[dof1], u1_cqc_ana) < 1e-10

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


print("=== 層3: 反力/基底剪力的簽名恆等式(混合門型剛架, 多支承) ===")
for kind in ('lumped', 'consistent'):
    fm = mixed_frame()
    mdm = eigen(fm, n_modes=5, mass=kind)
    for direction in ('x', 'y'):
        res3 = response_spectrum(mdm, lambda T: 3.0 + 2.0 / T, direction=direction, combine='SRSS')
        li = 0 if direction == 'x' else 1
        r_full = np.zeros(res3.modal_displacements.shape[1])
        for nid in fm.nodes:
            r_full[fm.dofs_of(nid)[li]] = 1.0
        worst = 0.0
        for i in range(mdm.n_modes if False else 5):
            s = sum(res3.modal_reactions[d][i] * r_full[d] for d in res3.modal_reactions)
            worst = max(worst, rel(s, -res3.modal_base_shear[i]))
        print(f"  [{kind}/{direction}] 5 個模態的 Σ(反力·r) vs -V_i, 最大相對誤差 = {worst:.2e}")
        assert worst < 1e-9

# ------------------------------------------------------------------
print("=== 層4: 獨立 K 線性方程式交叉驗證(不用 K@D=f 的捷徑, 直接解自由DOF方程式) ===")
fm = mixed_frame()
mdm = eigen(fm, n_modes=3, mass='consistent')
res4 = response_spectrum(mdm, lambda T: 4.0, direction='x', combine='SRSS')
asm = assemble_K(fm)
M = assemble_M(fm, 'consistent')
fixed = np.zeros(asm.n_dof, dtype=bool)
for s in fm.supports:
    for dof, val in zip(fm.dofs_of(s.node), (s.ux, s.uy, s.rot)):
        if val is not None:
            fixed[dof] = True
free = np.where(~fixed)[0]
for i in range(mdm.n_modes):
    fi = mdm.gamma['x'][i] * 4.0 * (M @ mdm.phi[:, i])
    u_free = np.linalg.solve(asm.K[np.ix_(free, free)], fi[free])
    u_full = np.zeros(asm.n_dof)
    u_full[free] = u_free
    err = np.max(np.abs(u_full - res4.modal_displacements[i])) / max(np.max(np.abs(u_full)), 1e-300)
    print(f"  模態{i + 1}: 獨立解 K u=f_i 與 D_i 的最大相對差 = {err:.2e}")
    assert err < 1e-9, "K D_i = f_i 這個捷徑跟直接解線性方程式不一致"

# ------------------------------------------------------------------
print("=== 層5: CQC 性質 ===")
om_close = np.array([10.0, 10.5, 30.0])
rho = cqc_rho(om_close, 0.05)
assert np.allclose(np.diag(rho), 1.0) and np.allclose(rho, rho.T)
assert rho[0, 1] > 0.3 and rho[0, 2] < 0.01, "頻率相近該有明顯相關, 頻率差很多該幾乎不相關"
print(f"  ρ(10,10.5)={rho[0, 1]:.4f}(密集模態, 應該明顯>0), ρ(10,30)={rho[0, 2]:.4f}(應該≈0)")
# 對照 Der Kiureghian & Rosenblueth 公式手算出的精確數值(不是只檢查「大概對」的區間), 這樣
# 分母漏乘一項這種分子分母都不為0、比值仍落在(0,1)區間的錯誤才抓得到
rho_hand = 0.8074520303382796
assert abs(rho[0, 1] - rho_hand) < 1e-12, f"ρ(10,10.5,ζ=0.05) 應該精確等於手算值 {rho_hand}, 實際 {rho[0, 1]}"
assert np.all(np.isfinite(rho)), "ρ 矩陣不能有 nan/inf(對角線 r=1 時分母不能是0)"

om_wide = np.array([10.0, 100.0, 1000.0])
vals = np.array([1.0, -0.6, 0.3])
s_srss = combine_srss(vals)
s_cqc = combine_cqc(vals, om_wide, 0.05)
print(f"  頻率分離良好時: SRSS = {s_srss:.6f}, CQC = {s_cqc:.6f}, 相對差 = {rel(s_cqc, s_srss):.2e}")
assert rel(s_cqc, s_srss) < 1e-3, "頻率分離良好時 CQC 應該收斂到 SRSS"

# ------------------------------------------------------------------
print("=== 層6: n_modes 截斷、明確拒絕 ===")
fm = mixed_frame()
mdm = eigen(fm, n_modes=5, mass='lumped')
r_all = response_spectrum(mdm, lambda T: 3.0, direction='x', n_modes=5)
r_1 = response_spectrum(mdm, lambda T: 3.0, direction='x', n_modes=1)
print(f"  累積質量比(總): n_modes=1 -> {r_1.cum_ratio_total:.4f}, n_modes=5 -> {r_all.cum_ratio_total:.4f}")
assert r_1.cum_ratio_total <= r_all.cum_ratio_total + 1e-12


def expect_error(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該raise ValueError")


expect_error("非法方向", lambda: response_spectrum(mdm, lambda T: 1.0, direction='z'))
expect_error("非法combine", lambda: response_spectrum(mdm, lambda T: 1.0, combine='average'))
expect_error("damping超出範圍", lambda: response_spectrum(mdm, lambda T: 1.0, damping=1.5))
expect_error("n_modes超出範圍", lambda: response_spectrum(mdm, lambda T: 1.0, n_modes=99))
expect_error("n_modes<1", lambda: response_spectrum(mdm, lambda T: 1.0, n_modes=0))
expect_error("spectrum回傳負值", lambda: response_spectrum(mdm, lambda T: -1.0))
expect_error("spectrum回傳nan", lambda: response_spectrum(mdm, lambda T: float('nan')))

print("\n全部通過: SDOF解析解、兩質量懸臂柔度矩陣參考、反力恆等式、獨立K解、CQC性質都吻合。")
