"""
驗證案例: frame2d.newmark -- 線性時程分析 (Newmark-β 直接積分) (動力分析 D4), 只依賴 numpy

層1 SDOF自由振動(無阻尼): 解析解 u(t)=u0·cos(ωt); Δt 收斂階數應該是 O(Δt²)(平均加速度法
    的標準結果, 見 Chopra)。
層2 SDOF自由振動(有阻尼): 解析解 u(t)=e^{-ζωt}(u0·cos(ω_d t)+...); 用對數遞減率從數值結果
    反推 ζ, 跟給定值比對, 順便驗證 damping_matrix 參數本身有作用。
層3 SDOF諧和強迫振動(有阻尼): 穩態振幅對動力放大係數解析公式, 相位落後角對解析公式。
層4 SDOF無阻尼階躍載重: 峰值 = 2倍靜位移, 峰值時刻 = 半週期(古典結果)。
層5 兩質量懸臂(無質量梁, 精確解析, 沿用 test_modal.py 的獨立柔度矩陣參考) 對脈衝載重的反應:
    用獨立的 Duhamel 積分(數值求積, 不呼叫 frame2d.modal/newmark 的任何邏輯)算出參考解,
    對照 Newmark 直接積分。
層6 無質量DOF回填的靜力平衡殘差檢查: 每個時間點的 u_m(t) 都要讓 K_md u_d + K_mm u_m − F_m
    精確為0(獨立於 Newmark 遞迴本身的殘差檢查)。
層7 明確拒絕。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.excitation import force_series_from_pattern, harmonic, pulse, step
from frame2d.newmark import newmark_integrate

E, I, A, L = 200e6, 8e-5, 1e-2, 3.0
M_TIP = 2.0


def cantilever_tip_mass(m=M_TIP):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 's')
    f.fix(0)
    f.add_mass(1, my=m)
    return f


KE = 3 * E * I / L**3
W = (KE / M_TIP) ** 0.5
T_NAT = 2 * np.pi / W


def damping_matrix_at(frame, node, direction, value):
    """在指定 DOF 的對角線放一個阻尼值, 其餘是0(SDOF阻尼測試用)。"""
    from frame2d.assembly import assemble_K
    n = assemble_K(frame).n_dof
    C = np.zeros((n, n))
    li = {'x': 0, 'y': 1, 'rot': 2}[direction]
    dof = frame.dofs_of(node)[li]
    C[dof, dof] = value
    return C


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: SDOF自由振動(無阻尼), 解析解與 Δt 收斂階數 ===")
u0 = 0.01
errs = {}
for k, dt in enumerate((T_NAT / 20, T_NAT / 40, T_NAT / 80, T_NAT / 160)):
    f = cantilever_tip_mass()
    nsteps = int(round(10 * T_NAT / dt))
    init = np.zeros(6)
    init[f.dofs_of(1)[1]] = u0
    res = newmark_integrate(f, dt, nsteps, initial_disp=init)
    tip = f.dofs_of(1)[1]
    u_exact = u0 * np.cos(W * res.t)
    err = float(np.max(np.abs(res.u[:, tip] - u_exact)))
    errs[dt] = err
    print(f"  Δt=T/{T_NAT / dt:.0f}: 最大誤差 = {err:.4e}")
dts = sorted(errs)
ratios = [errs[dts[i + 1]] / errs[dts[i]] for i in range(len(dts) - 1)]
print(f"  Δt減半時的誤差比值: {ratios} (理論值4.0, 平均加速度法是O(Δt²))")
assert all(3.7 < r < 4.3 for r in ratios), "Newmark(平均加速度法)的Δt收斂階數應該接近O(Δt²)"
assert errs[dts[0]] < 1e-4, "最小Δt時誤差應該很小"

# ------------------------------------------------------------------
print("=== 層2: SDOF自由振動(有阻尼), 對數遞減率反推 ζ ===")
zeta_true = 0.04
c = 2 * zeta_true * W * M_TIP           # C = 2ζωm
f = cantilever_tip_mass()
C = damping_matrix_at(f, 1, 'y', c)
dt = T_NAT / 100
nsteps = int(round(15 * T_NAT / dt))
init = np.zeros(6)
init[f.dofs_of(1)[1]] = u0
res = newmark_integrate(f, dt, nsteps, damping_matrix=C, initial_disp=init)
tip = f.dofs_of(1)[1]
wd = W * np.sqrt(1 - zeta_true**2)
u_exact = u0 * np.exp(-zeta_true * W * res.t) * (np.cos(wd * res.t) + zeta_true / np.sqrt(1 - zeta_true**2) * np.sin(wd * res.t))
err = float(np.max(np.abs(res.u[:, tip] - u_exact)))
print(f"  對解析解的最大誤差 = {err:.4e}(相對振幅 {err / u0:.2e})")
assert err / u0 < 5e-3

peaks_idx = [i for i in range(1, len(res.t) - 1)
            if res.u[i, tip] > res.u[i - 1, tip] and res.u[i, tip] > res.u[i + 1, tip] and res.u[i, tip] > 0]
p1, p2 = res.u[peaks_idx[0], tip], res.u[peaks_idx[3], tip]      # 隔3個週期的兩個峰值
delta = np.log(p1 / p2) / 3
zeta_fit = delta / np.sqrt(4 * np.pi**2 + delta**2)
print(f"  對數遞減率反推 ζ = {zeta_fit:.5f}, 給定 ζ = {zeta_true}")
assert abs(zeta_fit - zeta_true) < 1e-3

# ------------------------------------------------------------------
print("=== 層3: SDOF諧和強迫振動, 穩態振幅與相位對動力放大係數解析公式 ===")
zeta = 0.05
c = 2 * zeta * W * M_TIP
omega_ratio = 0.6
omega_f = omega_ratio * W
F0 = 10.0
f = cantilever_tip_mass()
C = damping_matrix_at(f, 1, 'y', c)
p = np.zeros(6)
p[f.dofs_of(1)[1]] = 1.0
force = force_series_from_pattern(p, harmonic(F0, omega_f))
dt = T_NAT / 200
n_cycles_force = 60
nsteps = int(round(n_cycles_force * (2 * np.pi / omega_f) / dt))
res = newmark_integrate(f, dt, nsteps, force=force, damping_matrix=C)
tip = f.dofs_of(1)[1]
u_static = F0 / KE
Rd = 1.0 / np.sqrt((1 - omega_ratio**2)**2 + (2 * zeta * omega_ratio)**2)
phase = np.arctan2(2 * zeta * omega_ratio, 1 - omega_ratio**2)
tail = res.t > res.t[-1] - 5 * (2 * np.pi / omega_f)            # 只看最後幾圈(穩態)
amp_num = (res.u[tail, tip].max() - res.u[tail, tip].min()) / 2
u_exact_tail = Rd * u_static * np.sin(omega_f * res.t[tail] - phase)
print(f"  數值穩態振幅 = {amp_num:.6f}, 解析 Rd·u_static = {Rd * u_static:.6f}, 相對差 = {rel(amp_num, Rd * u_static):.2e}")
assert rel(amp_num, Rd * u_static) < 1e-3
err_tail = float(np.max(np.abs(res.u[tail, tip] - u_exact_tail)))
print(f"  穩態波形(含相位)對解析解的最大誤差 = {err_tail:.4e}(相對振幅 {err_tail / (Rd * u_static):.2e})")
assert err_tail / (Rd * u_static) < 2e-2

# ------------------------------------------------------------------
print("=== 層4: SDOF無阻尼階躍載重, 峰值=2倍靜位移, 峰值時刻=半週期 ===")
F0 = 8.0
f = cantilever_tip_mass()
p = np.zeros(6)
p[f.dofs_of(1)[1]] = 1.0
force = force_series_from_pattern(p, step(F0, t_start=0.0))
dt = T_NAT / 400
nsteps = int(round(1.2 * T_NAT / dt))
res = newmark_integrate(f, dt, nsteps, force=force)
tip = f.dofs_of(1)[1]
u_static = F0 / KE
i_peak = int(np.argmax(res.u[:, tip]))
print(f"  數值峰值 = {res.u[i_peak, tip]:.6f}(解析 2·u_static = {2 * u_static:.6f}), "
      f"峰值時刻 = {res.t[i_peak]:.6f}s(解析 T/2 = {T_NAT / 2:.6f}s)")
assert rel(res.u[i_peak, tip], 2 * u_static) < 2e-3
assert rel(res.t[i_peak], T_NAT / 2) < 1e-2

# ------------------------------------------------------------------
print("=== 層5: 兩質量懸臂對脈衝載重的反應, 對獨立 Duhamel 積分(數值求積) ===")
m1, m2 = 3.0, 2.0
f2 = Frame2D()
f2.add_node(0, 0, 0).add_node(1, L / 2, 0).add_node(2, L, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's').add_member(1, 1, 2, 's')
f2.fix(0)
f2.support(1, ux=0.0).support(2, ux=0.0)
f2.add_mass(1, my=m1).add_mass(2, my=m2)

Fm = (L**3 / (E * I)) * np.array([[1 / 24, 5 / 48], [5 / 48, 1 / 3]])
Mm = np.diag([m1, m2])
lam_ana, V = np.linalg.eig(Fm @ Mm)
order = np.argsort(-lam_ana.real)
w_ana = 1.0 / np.sqrt(lam_ana.real[order])
V = V.real[:, order]
Vn = V / np.sqrt(np.einsum('ij,ik,kj->j', V, Mm, V))              # 對 M 正規化 (2 modes)
r_load = np.array([0.0, 1.0])                                      # 只在節點2(質量m2)施力
# 這是「直接施加的點力」, 不是地震(慣性)力: 模態力的振幅是 φᵢᵀp(直接對力向量取內積),
# **不是** φᵢᵀMp(那是地震參與係數 Γ 的定義, 這裡不適用——地震的等效力本身就是 -Mr·ag(t),
# 所以才會多一個 M; 直接點力 p(t) 不需要這個 M)。這裡一開始就是抄了 Γ 的公式, 導致算出來的
# 參考解剛好差 2 倍(用「單位點力靜位移 Fm@r_load」反推才抓到, 見這裡的獨立驗證)。
gamma_ana = Vn.T @ r_load

t_start, duration, amp = 0.02, 0.03, 12.0


def p_t(t):
    return amp if t_start <= t < t_start + duration else 0.0


T_END = 0.5
n_quad = 20000
tq = np.linspace(0, T_END, n_quad)
pq = np.array([p_t(tt) for tt in tq])


def duhamel(node_weight, wn):
    """無阻尼 Duhamel 積分(對每個模態座標), node_weight是 gamma_ana[mode]。"""
    integrand_cos = pq * np.cos(wn * tq)
    integrand_sin = pq * np.sin(wn * tq)
    A_t = np.concatenate([[0], np.cumsum((integrand_cos[1:] + integrand_cos[:-1]) / 2 * np.diff(tq))])
    B_t = np.concatenate([[0], np.cumsum((integrand_sin[1:] + integrand_sin[:-1]) / 2 * np.diff(tq))])
    q_t = (node_weight / wn) * (A_t * np.sin(wn * tq) - B_t * np.cos(wn * tq))
    return q_t


q1 = duhamel(gamma_ana[0], w_ana[0])
q2 = duhamel(gamma_ana[1], w_ana[1])
u2_ana = Vn[1, 0] * q1 + Vn[1, 1] * q2                              # 節點2的位移(獨立柔度矩陣座標)

f2b = Frame2D()
f2b.add_node(0, 0, 0).add_node(1, L / 2, 0).add_node(2, L, 0)
f2b.add_section('s', E=E, I=I, A=A)
f2b.add_member(0, 0, 1, 's').add_member(1, 1, 2, 's')
f2b.fix(0)
f2b.support(1, ux=0.0).support(2, ux=0.0)
f2b.add_mass(1, my=m1).add_mass(2, my=m2)
p2 = np.zeros(9)
p2[f2b.dofs_of(2)[1]] = 1.0
force2 = force_series_from_pattern(p2, pulse(amp, t_start, duration))
dt2 = T_END / n_quad
nsteps2 = n_quad - 1
res2 = newmark_integrate(f2b, dt2, nsteps2, force=force2)
u2_num = res2.u[:, f2b.dofs_of(2)[1]]
err = float(np.max(np.abs(u2_num - u2_ana)))
scale = float(np.max(np.abs(u2_ana)))
print(f"  節點2位移: 最大誤差 = {err:.4e}(相對峰值 {err / scale:.2e}), 峰值 數值={u2_num.max():.6f} 解析={u2_ana.max():.6f}")
assert err / scale < 5e-3

# ------------------------------------------------------------------
print("=== 層6: 無質量DOF回填的靜力平衡殘差(獨立檢查) ===")
from frame2d.assembly import assemble_K
asm = assemble_K(f2b)
K = asm.K
fixed = np.zeros(asm.n_dof, dtype=bool)
for s in f2b.supports:
    for dof, val in zip(f2b.dofs_of(s.node), (s.ux, s.uy, s.rot)):
        if val is not None:
            fixed[dof] = True
from frame2d.mass import assemble_M
M2 = assemble_M(f2b, 'lumped')
mless = np.where((~fixed) & (np.abs(M2).max(axis=1) == 0))[0]
print(f"  無質量DOF: {mless.tolist()}")
resid = (res2.u @ K.T)[:, mless]                                    # K @ u(t), 限制在無質量DOF那幾列
maxresid = float(np.max(np.abs(resid)))
scale2 = float(np.max(np.abs(res2.u @ K.T)))
print(f"  K@u(t) 在無質量DOF上的最大殘差 = {maxresid:.4e}(相對尺度 {scale2:.4e}, F_m(t)應該恆為0)")
assert maxresid < 1e-8 * scale2, "無質量DOF的回填不滿足靜力平衡(K@u在該處應該恆為0, 因為 F_m=0)"

# ------------------------------------------------------------------
print("=== 層7: 明確拒絕 ===")


def expect_error(label, fn, exc=ValueError):
    try:
        fn()
    except exc as e:
        print(f"  {label}: {exc.__name__} OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該raise {exc.__name__}")


f = cantilever_tip_mass()
expect_error("dt<=0", lambda: newmark_integrate(f, 0.0, 10))
expect_error("n_steps<1", lambda: newmark_integrate(f, 0.01, 0))
expect_error("beta<=0", lambda: newmark_integrate(f, 0.01, 10, beta=0.0))
expect_error("gamma<0.5", lambda: newmark_integrate(f, 0.01, 10, gamma=0.3))
expect_error("force(t)長度不對", lambda: newmark_integrate(f, 0.01, 5, force=lambda t: np.zeros(3)))


def bad_force(t):
    v = np.zeros(6)
    v[f.dofs_of(1)[2]] = 1.0                # 轉角DOF(無質量)上施力
    return v


expect_error("力施加在無質量DOF上", lambda: newmark_integrate(f, 0.01, 5, force=bad_force))
bad_C = np.eye(6)                            # 在無質量DOF上也有阻尼
expect_error("阻尼矩陣在無質量DOF上非零", lambda: newmark_integrate(f, 0.01, 5, damping_matrix=bad_C))
expect_error("initial_disp長度不對", lambda: newmark_integrate(f, 0.01, 5, initial_disp=np.zeros(3)))
bad_init = np.zeros(6)
bad_init[f.dofs_of(0)[1]] = 0.01             # 支承(固定)DOF上有非零初始位移
expect_error("固定DOF上有非零初始位移", lambda: newmark_integrate(f, 0.01, 5, initial_disp=bad_init))
f_nomass = Frame2D()
f_nomass.add_node(0, 0, 0).add_node(1, L, 0)
f_nomass.add_section('s', E=E, I=I, A=A)
f_nomass.add_member(0, 0, 1, 's')
f_nomass.fix(0)
expect_error("沒有任何質量", lambda: newmark_integrate(f_nomass, 0.01, 5))
# 注意: 跟 solve()/eigen() 不同, Newmark 積分不要求結構有支承——沒有支承的自由體(剛體模態)
# 對時程分析是良態的問題(K_hat 裡的質量項讓有效勁度矩陣仍然正定, 剛體模態只是自由飄移/轉動,
# 不是數值奇異), 所以這裡驗證的是「不會誤報」, 不是「應該報錯」
f_free = Frame2D()
f_free.add_node(0, 0, 0).add_node(1, L, 0)
f_free.add_section('s', E=E, I=I, A=A, rho=7.85)
f_free.add_member(0, 0, 1, 's')
res_free = newmark_integrate(f_free, 0.01, 5)
print(f"  沒有支承的自由體(剛體模態): 不報錯, 成功積分 {res_free.n_steps} 步 OK")

f_mech = Frame2D()                            # 真正的機構: release端沒有任何勁度支撐, K_mm奇異
f_mech.add_node(0, 0, 0).add_node(1, L, 0)
f_mech.add_section('s', E=E, I=I, A=A, rho=7.85)
f_mech.add_member(0, 0, 1, 's', release_i=True, release_j=True)
f_mech.fix(0)
f_mech.pin(1)                                  # 兩端都release, 桿件本身變成無彎矩勁度的機構
expect_error("release端形成機構(K_mm奇異)", lambda: newmark_integrate(f_mech, 0.01, 5))

print("\n全部通過: SDOF解析解(自由振動/阻尼/諧和/階躍)、獨立Duhamel積分、無質量DOF回填殘差都吻合。")
