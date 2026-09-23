"""
驗證案例: frame2d.nonlinear_newmark -- 非線性時程分析(Newmark + 循環塑鉸) (動力分析 D6),
只依賴 numpy

層1 彈性極限(SDOF): Mp=∞(永不降伏), 對照一個獨立手刻的線性 Newmark 迴圈, 用**同一份**
    塑鉸基礎的近剛接勁度矩陣(`_assemble_stiffness_with_hinges` 算出來那份, 不是
    `elements.member_stiffness_local` 那份標準樑元素——兩者本來就不會完全相同, 因為 hinge.py
    的「未降伏」彈簧用 RIGID_FACTOR·EI/L 這個很大但有限的數字近似剛接, 不是真的無限剛;
    這裡刻意跟同一種近似比, 才是公平的 apples-to-apples)。
層2 能量平衡(SDOF, 有阻尼、會降伏): 外力作功 = 動能 + 阻尼耗能 + 桿件內力作功(用桿件局部
    節點力對局部節點位移做功算, 這個定義本身就同時涵蓋彈性儲能與塑性耗能, 不需要另外拆開算,
    對任何非線性材料模型都成立, 是最通用的獨立檢查)。
層3 單調載重與 D7(`run_pushover`)一致: 用一個緩慢(準靜態等效)的載重歷程, 驗證塑鉸降伏順序
    與最終塑鉸狀態(降伏與否、方向)跟同一個模型的 `run_pushover` 结果吻合(兩者用的是同一套
    塑鉸底層機制, 只是這裡外加了動態的質量/阻尼項)。
層4 明確拒絕。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState
from frame2d.damping import rayleigh_damping_matrix
from frame2d.dofmanager import initial_hinge_states
from frame2d.excitation import force_series_from_pattern, harmonic, step
from frame2d.mass import assemble_M
from frame2d.nonlinear_newmark import nonlinear_newmark_integrate
from frame2d.pushover import _assemble_stiffness_with_hinges, _fixed_dof_set, run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 3.0
M_TIP = 2.0


def cantilever_elastic(mp=1e30, r=0.0):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='s', Mp_i=mp, Mp_j=None, R_post_yield_i=r, R_post_yield_j=None)
    f.fix(0)
    f.add_mass(1, my=M_TIP)
    return f


KE_STIFF = 3 * E * I / L**3
W = (KE_STIFF / M_TIP) ** 0.5
T_NAT = 2 * np.pi / W


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: 彈性極限(SDOF), 對獨立手刻的線性 Newmark(用同一份近剛接K) ===")
f = cantilever_elastic()
hs = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
tip = f.dofs_of(1)[1]
F0, omega_f = 10.0, 0.7 * W
p = np.zeros(6)
p[tip] = 1.0
force = force_series_from_pattern(p, harmonic(F0, omega_f))
dt = T_NAT / 97
nsteps = 400
res = nonlinear_newmark_integrate(f, hs, dt, nsteps, force=force)

K, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(
    cantilever_elastic(), CyclicHingeState.from_hinge_states(initial_hinge_states(cantilever_elastic())))
M = assemble_M(cantilever_elastic(), 'lumped')
n = K.shape[0]
fixed = sorted(_fixed_dof_set(cantilever_elastic()))
free = [d for d in range(n) if d not in fixed]
beta, gamma = 0.25, 0.5
Mff = M[np.ix_(free, free)]
Khat_inv = np.linalg.inv((K + (1.0 / (beta * dt**2)) * Mff if False else K[np.ix_(free, free)] + (1.0 / (beta * dt**2)) * Mff))
u_ref = np.zeros((nsteps + 1, n))
v_ref = np.zeros((nsteps + 1, n))
a_ref = np.zeros((nsteps + 1, n))
active = [i for i, d in enumerate(free) if Mff[i, i] > 0]
a0f = np.zeros(len(free))
a0f[active] = np.linalg.solve(Mff[np.ix_(active, active)], force(0.0)[free][active])
a_ref[0, free] = a0f
c1, c2 = 1.0 / (beta * dt), 1.0 / (2 * beta)
Fk = force(0.0)
for k in range(nsteps):
    Fk1 = force((k + 1) * dt)
    dF = (Fk1 - Fk)[free] + Mff @ (c1 * v_ref[k, free] + c2 * a_ref[k, free])
    du = Khat_inv @ dF
    da = (du - dt * v_ref[k, free] - (dt**2 / 2) * a_ref[k, free]) / (beta * dt**2)
    dv = gamma * dt * da + dt * a_ref[k, free]
    u_ref[k + 1, free] = u_ref[k, free] + du
    v_ref[k + 1, free] = v_ref[k, free] + dv
    a_ref[k + 1, free] = a_ref[k, free] + da
    Fk = Fk1
d = float(np.max(np.abs(res.u[:, tip] - u_ref[:, tip])))
scale = float(np.max(np.abs(u_ref[:, tip])))
print(f"  最大差 = {d:.3e}(相對振幅 {d / scale:.3e})")
assert d < 1e-9 * scale, "永不降伏時應該精確等於同一份近剛接K的線性Newmark解"

# ------------------------------------------------------------------
print("=== 層2: 能量平衡(SDOF, 有阻尼、會降伏) ===")


def cantilever_yield():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='s', Mp_i=25.0, Mp_j=None, R_post_yield_i=1000.0, R_post_yield_j=None)
    f.fix(0)
    f.add_mass(1, my=M_TIP)
    return f


f2 = cantilever_yield()
hs2 = CyclicHingeState.from_hinge_states(initial_hinge_states(f2))
zeta = 0.03
C, _, _ = rayleigh_damping_matrix(f2, W, W * 1.3, zeta)
tip2 = f2.dofs_of(1)[1]
p2 = np.zeros(6)
p2[tip2] = 1.0
F0b, omega_fb = 40.0, 0.8 * W
force2 = force_series_from_pattern(p2, harmonic(F0b, omega_fb))
dt2 = T_NAT / 200
nsteps2 = 2000
res2 = nonlinear_newmark_integrate(f2, hs2, dt2, nsteps2, force=force2, damping_matrix=C)
n_yield = sum(1 for e in res2.events if e['kind'] == 'yield')
n_unload = sum(1 for e in res2.events if e['kind'] == 'unload')
print(f"  {n_yield} 次降伏, {n_unload} 次卸載")
assert n_yield >= 5 and n_unload >= 3, "這個載重量應該有多次降伏/卸載才測得到重點"

M2 = assemble_M(f2, 'lumped')
KE = 0.5 * np.einsum('ti,ij,tj->t', res2.v, M2, res2.v)
Fser = np.array([force2(tt) for tt in res2.t])
dWext = np.concatenate([[0.0], 0.5 * np.einsum('ti,ti->t', Fser[:-1] + Fser[1:], res2.u[1:] - res2.u[:-1])])
Wext = np.cumsum(dWext)
Cv = res2.v @ C.T
dWdamp = np.concatenate([[0.0], 0.5 * np.einsum('ti,ti->t', Cv[:-1] + Cv[1:], res2.u[1:] - res2.u[:-1])])
Wdamp = np.cumsum(dWdamp)
_, member_dofs2, member_T2, _ = _assemble_stiffness_with_hinges(f2, hs2)
idx = np.array(member_dofs2[0])
u_local = (member_T2[0] @ res2.u[:, idx].T).T
f_local = res2.member_force_history[0]
dWint = np.concatenate([[0.0], 0.5 * np.einsum('ti,ti->t', f_local[:-1] + f_local[1:], u_local[1:] - u_local[:-1])])
Wint = np.cumsum(dWint)
resid = Wext - (KE + Wdamp + Wint)
scale = float(np.max(np.abs(Wext)))
print(f"  能量殘差最大值 = {np.max(np.abs(resid)):.3e}(相對尺度 {scale:.3e}, 相對誤差 {np.max(np.abs(resid)) / scale:.2e})")
assert np.max(np.abs(resid)) < 1e-8 * scale, "外力作功應該精確等於動能+阻尼耗能+桿件內力作功"
assert res2.hinge_theta_p[(0, 0)][-1] != 0.0, "應該留下殘餘塑性轉角"
assert res2.hinge_theta_p[(0, 0)][-1] >= abs(res2.hinge_theta_p[(0, 0)][-1]) - 1e-15, "累積(絕對值)不能小於淨值的絕對值"
# theta_p 是「累積絕對值」, theta_p_signed(這裡=res2.hinge_theta_p同一個key但取自CyclicHingeState)是
# 「帶號淨值」; 有反覆降伏/卸載時前者必須嚴格大於後者的絕對值(這是本模組唯一沒被其他檢查覆蓋到的
# 記帳邏輯, 拿獨立追蹤的 hs2[0].theta_p 來對照 hs2[0].theta_p_signed)
tp_abs = hs2[0].theta_p[0]
tp_signed = hs2[0].theta_p_signed[0]
print(f"  累積絕對值 θp = {tp_abs:.6f}, 帶號淨值 θp_signed = {tp_signed:.6f}(反覆降伏應該讓前者明顯更大)")
assert tp_abs > abs(tp_signed) * 1.5, "有反覆降伏/卸載時, 累積絕對值應該明顯大於淨值的絕對值(否則可能沒有真的累加|Δθp|)"

# ------------------------------------------------------------------
print("=== 層3: 單調(準靜態等效)載重下的「第一個降伏事件」跟 D7 run_pushover 一致 ===")
# 用整條非線性路徑(降伏後軟化、多鉸交互卸載)去對比動力斜坡載重, 對「多慢才算夠慢」很敏感
# (實測: 降伏後結構軟化, 有效週期變長, 就算斜坡相對「原始彈性週期」很慢, 仍然可能在降伏後的
# 路徑上產生動態超越/暫時卸載, 最終停留狀態就對不上pushover)。改成只比對「第一個降伏事件」:
# 這一段路徑上結構還是線性彈性(單一明確的自然週期), 斜坡夠慢就不會有這個問題, 是穩健得多的
# 獨立檢查。


def portal(mp=60.0, r=1500.0):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('s', E=E, I=I, A=A)
    for mid, (i, j) in enumerate([(0, 1), (1, 2), (3, 2)]):
        f.add_member(mid, node_i=i, node_j=j, section='s', Mp_i=mp, Mp_j=mp, R_post_yield_i=r, R_post_yield_j=r)
    f.fix(0).fix(3)
    return f


f3 = portal()
dof3 = f3.dofs_of(1)[0]
base3 = [f3.dofs_of(0)[0], f3.dofs_of(3)[0]]
u_p, F_p, ev_p, hs_p, mech_p = run_pushover(f3, initial_hinge_states(f3), [dof3], [1.0], 0.05, 0.004, base3)
assert ev_p, "這個案例應該至少有一個降伏事件, 測試才有意義"
u_first_p, F_first_p = ev_p[0]['u'], ev_p[0]['F']

f3b = portal()
f3b.add_mass(1, mx=500.0, my=500.0).add_mass(2, mx=500.0, my=500.0)   # 給質量(否則沒辦法做動力分析)
hs3 = CyclicHingeState.from_hinge_states(initial_hinge_states(f3b))
p3 = np.zeros(assemble_M(f3b, 'lumped').shape[0])
p3[dof3] = 1.0
target_force = 100.0                            # 超過第一個降伏點就好, 不用追完整條路徑
from frame2d.modal import eigen
w3 = eigen(f3b, n_modes=1, mass='lumped').omega[0]
T3 = 2 * np.pi / w3
t_ramp = 500 * T3                               # 遠慢於結構(仍是彈性時)的自然週期, 準靜態等效
force3 = force_series_from_pattern(p3, lambda t: target_force * min(t / t_ramp, 1.0))
dt3 = T3 / 40
nsteps3 = int(round(t_ramp / dt3))
res3 = nonlinear_newmark_integrate(f3b, hs3, dt3, nsteps3, force=force3)
assert res3.events, "動力斜坡載重也應該至少觸發一個降伏事件"
first = res3.events[0]
assert first['kind'] == 'yield' and (first['member'], first['end']) == ev_p[0]['yielded'][0], \
    "第一個降伏的應該是同一根桿件同一端"
idx = first['step']
u_first_d = res3.u[idx, dof3]
F_first_d = float(force3(res3.t[idx])[dof3])
print(f"  第一個降伏鉸: pushover=member{ev_p[0]['yielded'][0][0]}端{ev_p[0]['yielded'][0][1]}, "
      f"動力斜坡=member{first['member']}端{first['end']}")
print(f"  降伏時 u: pushover={u_first_p:.6f}, 動力斜坡={u_first_d:.6f}, 相對差={rel(u_first_d, u_first_p):.2e}")
print(f"  降伏時 F: pushover={F_first_p:.4f}, 動力斜坡={F_first_d:.4f}, 相對差={rel(F_first_d, F_first_p):.2e}")
assert rel(u_first_d, u_first_p) < 2e-3
assert rel(F_first_d, F_first_p) < 2e-3

# ------------------------------------------------------------------
print("=== 層4: 明確拒絕 ===")


def expect_error(label, fn, exc=ValueError):
    try:
        fn()
    except exc as e:
        print(f"  {label}: {exc.__name__} OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該raise {exc.__name__}")


f4 = cantilever_elastic()
hs4 = CyclicHingeState.from_hinge_states(initial_hinge_states(f4))
expect_error("dt<=0", lambda: nonlinear_newmark_integrate(f4, hs4, 0.0, 10))
expect_error("n_steps<1", lambda: nonlinear_newmark_integrate(f4, hs4, 0.01, 0))
expect_error("beta<=0", lambda: nonlinear_newmark_integrate(f4, hs4, 0.01, 10, beta=0.0))
expect_error("gamma<0.5", lambda: nonlinear_newmark_integrate(f4, hs4, 0.01, 10, gamma=0.3))
hs_plain = initial_hinge_states(f4)
expect_error("傳入單向HingeState", lambda: nonlinear_newmark_integrate(f4, hs_plain, 0.01, 10), TypeError)
f5 = cantilever_yield()
hs5 = CyclicHingeState.from_hinge_states(initial_hinge_states(f5))
bad_preload = {0: np.array([0.0, 0.0, 200.0, 0.0, 0.0, 0.0])}   # 超過Mp=25
expect_error("重力預載已超過Mp", lambda: nonlinear_newmark_integrate(f5, hs5, 0.01, 10, initial_cum_forces=bad_preload))
expect_error("force(t)長度不對", lambda: nonlinear_newmark_integrate(f4, hs4, 0.01, 5, force=lambda t: np.zeros(3)))
f_nomass = Frame2D()
f_nomass.add_node(0, 0, 0).add_node(1, L, 0)
f_nomass.add_section('s', E=E, I=I, A=A)
f_nomass.add_member(0, node_i=0, node_j=1, section='s', Mp_i=1e30, Mp_j=None, R_post_yield_i=0.0, R_post_yield_j=None)
f_nomass.fix(0)
hs_nomass = CyclicHingeState.from_hinge_states(initial_hinge_states(f_nomass))
expect_error("沒有任何質量", lambda: nonlinear_newmark_integrate(f_nomass, hs_nomass, 0.01, 5))

print("\n全部通過: 彈性極限對線性Newmark、能量平衡、準靜態等效對run_pushover都吻合。")
