"""
驗證案例: frame2d.damping -- Rayleigh 阻尼 (動力分析 D5), 只依賴 numpy

層1 rayleigh_coefficients(): 解出的 (α,β) 代回 rayleigh_damping_ratio() 精確重現給定的
    ζᵢ、ζⱼ(包括 ζᵢ≠ζⱼ 的一般情況); ζⱼ=None 化簡成經典公式
    α=2ζωᵢωⱼ/(ωᵢ+ωⱼ), β=2ζ/(ωᵢ+ωⱼ)。
層2 SDOF: rayleigh_damping_matrix() 精確等於古典粘滯阻尼 c=2ζωm(代數恆等式, 不是數值近似)。
層3 SDOF damped自由振動: 對解析解(跟 test_newmark.py 層2 同一條公式, 這裡改用
    rayleigh_damping_matrix() 產生的 C, 驗證兩條路徑一致)。
層4 MDOF(6元素懸臂梁)模態投影檢核: 用 modal.eigen() 獨立算出的模態形狀, 把 Rayleigh 阻尼矩陣
    投影到**每一個**模態(不只是拿來當控制頻率的那兩個), ζₙ_投影=φₙᵀCφₙ/(2ωₙ) 對照
    rayleigh_damping_ratio() 的解析公式, 逐一比對。
層5 SDOF對諧和地震輸入的穩態反應: 跟 test_newmark.py 層3 同一套動力放大係數公式, 只是這次
    的力是 -m·ag(t) 不是直接施力。
層6 絕對加速度殘差: m·a_abs+c·v_rel+k·u_rel 應該精確為0(這是運動方程式本身的恆等式, 獨立於
    怎麼算出 a_abs, 拿來驗證 absolute_acceleration() 的重建公式)。
層7 明確拒絕。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.damping import rayleigh_coefficients, rayleigh_damping_matrix, rayleigh_damping_ratio
from frame2d.excitation import absolute_acceleration, ground_motion_force, harmonic
from frame2d.modal import eigen
from frame2d.newmark import newmark_integrate

E, I, A, RHO, L = 200e6, 8e-5, 1e-2, 7.85, 3.0
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


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: rayleigh_coefficients() 解析驗證 ===")
wi, wj = 10.0, 55.0
for zi, zj in ((0.05, 0.05), (0.03, 0.08), (0.10, 0.02)):
    alpha, beta = rayleigh_coefficients(wi, wj, zi, zj if zj != zi else None)
    got_i = rayleigh_damping_ratio(alpha, beta, wi)
    got_j = rayleigh_damping_ratio(alpha, beta, wj)
    print(f"  ζᵢ={zi}, ζⱼ={zj}: 反推 ζᵢ={got_i:.6f}, ζⱼ={got_j:.6f}")
    assert abs(got_i - zi) < 1e-12 and abs(got_j - zj) < 1e-12
alpha, beta = rayleigh_coefficients(wi, wj, 0.05)
alpha_classic = 2 * 0.05 * wi * wj / (wi + wj)
beta_classic = 2 * 0.05 / (wi + wj)
print(f"  ζⱼ=None(同一ζ): α={alpha:.6f}(經典公式{alpha_classic:.6f}), β={beta:.6f}(經典公式{beta_classic:.6f})")
assert abs(alpha - alpha_classic) < 1e-12 and abs(beta - beta_classic) < 1e-12

# ------------------------------------------------------------------
print("=== 層2: SDOF, C 精確等於古典 c=2ζωm ===")
zeta = 0.04
f = cantilever_tip_mass()
C, alpha, beta = rayleigh_damping_matrix(f, W, W * 1.37, zeta, zeta)   # 兩個"控制頻率"對SDOF沒差, 只要都=ζ
tip = f.dofs_of(1)[1]
c_classic = 2 * zeta * W * M_TIP
print(f"  C[tip,tip] = {C[tip, tip]:.8f}, 古典 2ζωm = {c_classic:.8f}")
assert rel(C[tip, tip], c_classic) < 1e-10
assert np.count_nonzero(C) == 1, "SDOF只有一個動態DOF, C矩陣應該只有一個非零項"

# ------------------------------------------------------------------
print("=== 層3: SDOF damped自由振動對解析解(用 rayleigh_damping_matrix 產生的 C) ===")
u0 = 0.01
dt = T_NAT / 100
nsteps = int(round(15 * T_NAT / dt))
init = np.zeros(6)
init[tip] = u0
res = newmark_integrate(f, dt, nsteps, damping_matrix=C, initial_disp=init)
wd = W * np.sqrt(1 - zeta**2)
u_exact = u0 * np.exp(-zeta * W * res.t) * (np.cos(wd * res.t) + zeta / np.sqrt(1 - zeta**2) * np.sin(wd * res.t))
err = float(np.max(np.abs(res.u[:, tip] - u_exact)))
print(f"  最大誤差 = {err:.4e}(相對振幅 {err / u0:.2e})")
assert err / u0 < 5e-3


# ------------------------------------------------------------------
def cantilever_mdof(n=6):
    f = Frame2D()
    for i in range(n + 1):
        f.add_node(i, L * i / n, 0)
    f.add_section('s', E=E, I=I, A=A, rho=RHO)
    for i in range(n):
        f.add_member(i, i, i + 1, 's')
    f.fix(0)
    for i in range(1, n + 1):
        f.support(i, ux=0.0)
    return f


print("=== 層4: MDOF模態投影檢核(獨立於Newmark積分本身) ===")
f6 = cantilever_mdof()
md = eigen(f6, mass='lumped')          # 全部模態
i1, i3 = 0, 2                          # 拿第1、第3模態當控制頻率
w1, w3 = md.omega[i1], md.omega[i3]
zeta_target = 0.05
C, alpha, beta = rayleigh_damping_matrix(f6, w1, w3, zeta_target, mass_kind='lumped')
print(f"  控制頻率: ω1={w1:.3f}, ω3={w3:.3f}, α={alpha:.6f}, β={beta:.6f}")
worst = 0.0
for n in range(md.n_modes):
    phi = md.phi[:, n]
    zeta_proj = float(phi @ C @ phi) / (2 * md.omega[n])
    zeta_ana = rayleigh_damping_ratio(alpha, beta, md.omega[n])
    d = rel(zeta_proj, zeta_ana)
    worst = max(worst, d)
    mark = " <- 控制頻率" if n in (i1, i3) else ""
    print(f"  模態{n + 1}(ω={md.omega[n]:8.2f}): 投影ζ={zeta_proj:.6f}, 解析ζ(ω)={zeta_ana:.6f}, 相對差={d:.2e}{mark}")
assert worst < 1e-9, "模態投影的阻尼比應該跟解析公式 rayleigh_damping_ratio 逐一相符"
assert abs(rayleigh_damping_ratio(alpha, beta, w1) - zeta_target) < 1e-12
assert abs(rayleigh_damping_ratio(alpha, beta, w3) - zeta_target) < 1e-12
assert rayleigh_damping_ratio(alpha, beta, (w1 + w3) / 2) < zeta_target, "兩個控制頻率之間的阻尼比應該比目標值低(Rayleigh阻尼的已知特性)"

# ------------------------------------------------------------------
print("=== 層5: SDOF對諧和地震輸入的穩態反應, 對動力放大係數解析公式 ===")
zeta = 0.05
f = cantilever_tip_mass()
C, _, _ = rayleigh_damping_matrix(f, W, W * 1.001, zeta)
A0 = 3.0                                # 地面加速度振幅
omega_ratio = 0.65
omega_f = omega_ratio * W
ag = harmonic(A0, omega_f)
force = ground_motion_force(f, 'y', ag)
dt = T_NAT / 200
n_cycles = 60
nsteps = int(round(n_cycles * (2 * np.pi / omega_f) / dt))
res = newmark_integrate(f, dt, nsteps, force=force, damping_matrix=C)
u_static_equiv = M_TIP * A0 / KE        # 等效靜位移: F0=m·A0, u_static=F0/k
Rd = 1.0 / np.sqrt((1 - omega_ratio**2)**2 + (2 * zeta * omega_ratio)**2)
tail = res.t > res.t[-1] - 5 * (2 * np.pi / omega_f)
amp_num = (res.u[tail, tip].max() - res.u[tail, tip].min()) / 2
print(f"  數值穩態振幅 = {amp_num:.6f}, 解析 Rd·u_static = {Rd * u_static_equiv:.6f}, "
      f"相對差 = {rel(amp_num, Rd * u_static_equiv):.2e}")
assert rel(amp_num, Rd * u_static_equiv) < 1e-3

# ------------------------------------------------------------------
print("=== 層6: 絕對加速度殘差 m·a_abs+c·v_rel+k·u_rel ≈ 0(運動方程式本身的恆等式) ===")
a_abs = absolute_acceleration(res, 'y', ag)
resid = M_TIP * a_abs[:, tip] + C[tip, tip] * res.v[:, tip] + KE * res.u[:, tip]
scale = float(np.max(np.abs(M_TIP * a_abs[:, tip])))
print(f"  最大殘差 = {np.max(np.abs(resid)):.4e}(相對尺度 {scale:.4e})")
assert np.max(np.abs(resid)) < 1e-8 * scale
rot = f.dofs_of(1)[2]                   # 轉角DOF: 影響向量r在這裡是0, 不該被地面加速度影響
diff_rot = float(np.max(np.abs(a_abs[:, rot] - res.a[:, rot])))
print(f"  轉角DOF(r=0): a_abs 應該精確等於 a_rel, 最大差 = {diff_rot:.4e}")
assert diff_rot == 0.0, "地面加速度不該影響跟激振方向無關的分量(這裡用來抓「絕對加速度重建" \
                        "誤用純量廣播、沒有乘上影響向量r」這類錯誤)"

# ------------------------------------------------------------------
print("=== 層7: 明確拒絕 ===")


def expect_error(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該raise ValueError")


expect_error("wi<=0", lambda: rayleigh_coefficients(-1.0, 10.0, 0.05))
expect_error("wj<=0", lambda: rayleigh_coefficients(10.0, 0.0, 0.05))
expect_error("wi==wj", lambda: rayleigh_coefficients(10.0, 10.0, 0.05))
expect_error("zeta為負", lambda: rayleigh_coefficients(10.0, 20.0, -0.01))
expect_error("rayleigh_damping_ratio: omega<=0", lambda: rayleigh_damping_ratio(1.0, 0.001, -5.0))
expect_error("ground_motion_force: 非法方向", lambda: ground_motion_force(cantilever_tip_mass(), 'z', harmonic(1.0, 1.0)))
expect_error("absolute_acceleration: 非法方向", lambda: absolute_acceleration(res, 'z', ag))

print("\n全部通過: rayleigh_coefficients解析驗證、SDOF代數恆等式、MDOF模態投影、諧和地震動力放大、絕對加速度殘差都吻合。")
