"""
驗證案例: frame2d.seismic -- 非線性地震反應一站式入口 (動力分析 D8), 只依賴 numpy

這個模組不含新的求解邏輯, 純粹是把 D2/D5/D6/D7 已經各自驗證過的元件組裝起來, 所以測試重點是
「組裝有沒有接對」, 不是重新驗證底層數學(那些在各自的 test_modal.py / test_damping.py /
test_nonlinear_newmark.py / test_cyclic.py 已經做過)。

層1 能量平衡(無重力預載): 外力作功 = 動能+阻尼耗能+桿件內力作功, 精確成立。**這一層抓到一個
    真實 bug**: `energy_balance()` 第一版把每根桿件的功增量累加後忘記對時間做 `np.cumsum()`
    (Wext/Wdamp 都有 cumsum, Wint 少了), 導致回傳的是「每個時間點當下的瞬時累加值」不是「累積
    到那個時間點的總和」, 殘差高達外力功尺度的 82%。修正後殘差 1e-10 量級(相對 5e-14)。
層2 能量平衡(有重力預載): 驗證重力預載狀態不影響能量平衡式本身(重力在純相對動態運動下
    不作淨功這個假設是對的)。
層3 `seismic_analysis()` 只是接線, 不是另一套邏輯: 手動照著 D2/D5/D6/D7 的元件各自組一次,
    結果要跟 `seismic_analysis()` 逐項相同。
層4 Rayleigh 阻尼比對照解析公式(`damping.rayleigh_damping_ratio`), 包含兩個控制模態剛好
    等於目標值、其他模態依公式計算。
層5 `peak_displacement()`/`peak_drift_ratio()` 對照直接從陣列索引算出來的值。
層6 明確拒絕。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState
from frame2d.damping import rayleigh_damping_matrix, rayleigh_damping_ratio
from frame2d.dofmanager import initial_hinge_states
from frame2d.excitation import ground_motion_force
from frame2d.modal import eigen
from frame2d.nonlinear_newmark import nonlinear_newmark_integrate
from frame2d.pushover import apply_gravity
from frame2d.seismic import seismic_analysis

E, I, A = 200e6, 8e-5, 1e-2
G = 9.80665


def portal(mp=60.0, r=1500.0):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('s', E=E, I=I, A=A)
    for mid, (i, j) in enumerate([(0, 1), (1, 2), (3, 2)]):
        f.add_member(mid, node_i=i, node_j=j, section='s', Mp_i=mp, Mp_j=mp, R_post_yield_i=r, R_post_yield_j=r)
    f.fix(0).fix(3)
    f.add_mass(1, mx=500.0, my=500.0).add_mass(2, mx=500.0, my=500.0)
    return f


def ag(t):
    return 0.4 * G * np.sin(2 * np.pi * 1.5 * t) * np.exp(-0.25 * t)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: 能量平衡(無重力預載) ===")
f = portal()
res = seismic_analysis(f, ag, direction='x', zeta=0.05)
assert len(res.events) > 0, "這個案例應該至少有降伏事件才測得到重點"
eb = res.energy_balance()
resid = eb['Wext'] - (eb['KE'] + eb['Wdamp'] + eb['Wint'])
scale = float(np.max(np.abs(eb['Wext'])))
print(f"  {len(res.events)} 次事件, 能量殘差最大值 = {np.max(np.abs(resid)):.3e}(相對尺度 {scale:.3e}, "
      f"相對誤差 {np.max(np.abs(resid)) / scale:.2e})")
assert np.max(np.abs(resid)) < 1e-8 * scale, "外力作功應該精確等於動能+阻尼耗能+桿件內力作功"
assert eb['Wint'][-1] > 0.2 * eb['Wext'][-1], "有明顯降伏時, 累積內力作功應該佔最終總能量相當比例(用Wext[-1]而非時間上的最大值scale比較, 因為Wext不是單調的)"

# ------------------------------------------------------------------
print("=== 層2: 能量平衡(有重力預載) ===")
f2 = portal()
f2.distributed_load(1, w=-10.0)
res2 = seismic_analysis(f2, ag, direction='x', zeta=0.05, apply_gravity_loads=True)
assert res2.gravity_result is not None, "有載重時應該做了重力預載"
eb2 = res2.energy_balance()
resid2 = eb2['Wext'] - (eb2['KE'] + eb2['Wdamp'] + eb2['Wint'])
scale2 = float(np.max(np.abs(eb2['Wext'])))
print(f"  能量殘差最大值 = {np.max(np.abs(resid2)):.3e}(相對尺度 {scale2:.3e}, 相對誤差 {np.max(np.abs(resid2)) / scale2:.2e})")
assert np.max(np.abs(resid2)) < 1e-8 * scale2, "有重力預載時外力作功一樣應該精確平衡"

f3 = portal()
res3 = seismic_analysis(f3, ag, direction='x', zeta=0.05, apply_gravity_loads=False)
assert res3.gravity_result is None, "apply_gravity_loads=False 時不應該做重力預載"

# ------------------------------------------------------------------
print("=== 層3: seismic_analysis() 只是接線, 手動組裝應該逐項相同 ===")
f4 = portal()
hs0 = initial_hinge_states(f4)
md = eigen(f4, n_modes=None, mass='lumped')
w1, w3 = md.omega[0], md.omega[2]
C, alpha, beta = rayleigh_damping_matrix(f4, w1, w3, 0.05, mass_kind='lumped')
hs = CyclicHingeState.from_hinge_states(hs0)
dt = md.period[0] / 50.0
nsteps = int(round(10.0 / dt))
force = ground_motion_force(f4, 'x', ag, mass_kind='lumped')
nl_ref = nonlinear_newmark_integrate(f4, hs, dt, nsteps, force=force, mass_kind='lumped', damping_matrix=C)

f5 = portal()
res5 = seismic_analysis(f5, ag, direction='x', zeta=0.05)
assert res5.nl.u.shape == nl_ref.u.shape
assert float(np.max(np.abs(res5.nl.u - nl_ref.u))) < 1e-12 * float(np.max(np.abs(nl_ref.u)))
assert float(np.max(np.abs(res5.nl.a - nl_ref.a))) < 1e-9 * float(np.max(np.abs(nl_ref.a)))
assert len(res5.events) == len(nl_ref.events)
print(f"  手動組裝 vs seismic_analysis(): 位移最大相對差 "
      f"{float(np.max(np.abs(res5.nl.u - nl_ref.u))) / float(np.max(np.abs(nl_ref.u))):.2e}, "
      f"事件數 {len(res5.events)} 對 {len(nl_ref.events)}")
assert abs(res5.alpha - alpha) < 1e-12 and abs(res5.beta - beta) < 1e-12

# ------------------------------------------------------------------
print("=== 層4: Rayleigh 阻尼比對解析公式 ===")
f6 = portal()
res6 = seismic_analysis(f6, ag, direction='x', zeta=0.06, damping_modes=(0, 2))
w1_6, w3_6 = res6.modal.omega[0], res6.modal.omega[2]
zeta1 = rayleigh_damping_ratio(res6.alpha, res6.beta, w1_6)
zeta3 = rayleigh_damping_ratio(res6.alpha, res6.beta, w3_6)
print(f"  控制模態1: ζ={zeta1:.6f}(目標0.06), 控制模態3: ζ={zeta3:.6f}(目標0.06)")
assert abs(zeta1 - 0.06) < 1e-10 and abs(zeta3 - 0.06) < 1e-10
for n in range(res6.modal.n_modes):
    zeta_n = rayleigh_damping_ratio(res6.alpha, res6.beta, res6.modal.omega[n])
    zeta_proj = float(res6.modal.phi[:, n] @ res6._C @ res6.modal.phi[:, n]) / (2 * res6.modal.omega[n])
    assert rel(zeta_proj, zeta_n) < 1e-9, f"模態{n + 1}: 投影阻尼比應該對得上解析公式"
print(f"  全部 {res6.modal.n_modes} 個模態的投影阻尼比都對得上解析公式")

# ------------------------------------------------------------------
print("=== 層5: peak_displacement / peak_drift_ratio 對照直接陣列索引 ===")
f7 = portal()
res7 = seismic_analysis(f7, ag, direction='x', zeta=0.05)
dof_ux1 = f7.dofs_of(1)[0]
want_peak = float(np.max(np.abs(res7.u[:, dof_ux1])))
got_peak = res7.peak_displacement(1, 'x')
print(f"  peak_displacement(1,'x') = {got_peak:.8f}, 直接陣列索引 = {want_peak:.8f}")
assert got_peak == want_peak
dof_ux0 = f7.dofs_of(0)[0]
want_drift = float(np.max(np.abs(res7.u[:, dof_ux1] - res7.u[:, dof_ux0]))) / 4.0
got_drift = res7.peak_drift_ratio(1, 0, 4.0, 'x')
print(f"  peak_drift_ratio(1,0,4.0,'x') = {got_drift:.8f}, 直接陣列索引 = {want_drift:.8f}")
assert got_drift == want_drift

# ------------------------------------------------------------------
print("=== 層6: 明確拒絕 ===")


def expect_error(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該raise ValueError")


f_nohinge = Frame2D()
f_nohinge.add_node(0, 0, 0).add_node(1, 4, 0)
f_nohinge.add_section('s', E=E, I=I, A=A)
f_nohinge.add_member(0, node_i=0, node_j=1, section='s')
f_nohinge.fix(0)
f_nohinge.add_mass(1, my=2.0)
expect_error("沒有任何塑鉸容量", lambda: seismic_analysis(f_nohinge, ag, direction='x'))
expect_error("非法方向", lambda: seismic_analysis(portal(), ag, direction='z'))
expect_error("damping_modes超出模態數", lambda: seismic_analysis(portal(), ag, damping_modes=(0, 999)))

print("\n全部通過: 能量平衡(含重力預載)、seismic_analysis與手動組裝一致、Rayleigh阻尼比、峰值反應都吻合。")
