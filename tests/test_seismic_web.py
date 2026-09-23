"""
驗證案例: frame2d.seismic 的網頁入口 (動力分析 D6b) -- sine_pulse_ground_motion、
custom_ground_motion、nonlinear_seismic_web_analysis()、seismic_to_dict()

層1 sine_pulse_ground_motion: 公式本身(振幅、頻率、衰減)逐點手算比對; 明確拒絕。
層2 custom_ground_motion: 分段線性內插對 np.interp 獨立手算比對, 邊界外夾在端點值; 自動排序;
    明確拒絕。
層3 nonlinear_seismic_web_analysis()+seismic_to_dict(): JSON 可序列化(無 nan)、跟直接呼叫
    `seismic_analysis()` 逐項相同(不是另一套邏輯)、frame_idx 抽稀正確(首尾保留、數量不超過
    max_frames)、能量平衡在JSON化後仍然成立。
層4 明確拒絕: ground_motion_type打錯、缺參數、超過步數上限、找不到控制節點。
"""
import json

import numpy as np

from frame2d import Frame2D
from frame2d.seismic import (custom_ground_motion, nonlinear_seismic_web_analysis,
                             seismic_analysis, seismic_to_dict, sine_pulse_ground_motion)

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0
G = 9.80665


def portal(mp=60e3, r=0.1):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    for mid, (i, j, L) in enumerate([(0, 1, H), (1, 2, W), (3, 2, H)]):
        R = r * E * I / L
        f.add_member(mid, node_i=i, node_j=j, section='H', Mp_i=mp, Mp_j=mp, R_post_yield_i=R, R_post_yield_j=R)
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: sine_pulse_ground_motion ===")
amp, freq, decay = 0.4, 2.0, 0.3
ag = sine_pulse_ground_motion(amp, freq, decay)
for t in (0.0, 0.1, 0.37, 1.0):
    want = amp * G * np.sin(2 * np.pi * freq * t) * np.exp(-decay * t)
    assert abs(ag(t) - want) < 1e-12
print("  4 個測試點跟公式逐點相同")
ag0 = sine_pulse_ground_motion(amp, freq, decay=0.0)
assert abs(ag0(5.0) - amp * G * np.sin(2 * np.pi * freq * 5.0)) < 1e-9, "decay=0應該不衰減"


def expect_err(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:36]}...)")
        return
    raise AssertionError(f"{label}: 應該報錯")


expect_err("amplitude_g<=0", lambda: sine_pulse_ground_motion(0.0, 2.0))
expect_err("freq_hz<=0", lambda: sine_pulse_ground_motion(0.4, 0.0))
expect_err("decay<0", lambda: sine_pulse_ground_motion(0.4, 2.0, -0.1))

# ------------------------------------------------------------------
print("=== 層2: custom_ground_motion ===")
pts = [(0.0, 0.0), (0.5, 0.4), (1.0, -0.2), (2.0, 0.0)]
ag_c = custom_ground_motion(pts)
ts_ref = np.array([p[0] for p in pts])
ags_ref = np.array([p[1] for p in pts]) * G
for t in (0.0, 0.25, 0.5, 0.75, 1.5, 2.0):
    want = float(np.interp(t, ts_ref, ags_ref))
    assert abs(ag_c(t) - want) < 1e-9, f"t={t}"
print("  6 個測試點跟 np.interp(乘上g) 逐點相同")
assert ag_c(-1.0) == 0.0 * G and ag_c(10.0) == 0.0 * G, "範圍外應該夾在端點值(不外插)"
ag_unsorted = custom_ground_motion([(1.0, 0.2), (0.0, 0.0), (2.0, -0.1)])
assert abs(ag_unsorted(0.5) - 0.1 * G) < 1e-9, "輸入順序打亂應該自動依時間排序"
expect_err("只有1個點", lambda: custom_ground_motion([(0.0, 0.1)]))
expect_err("時間有負值", lambda: custom_ground_motion([(-1.0, 0.1), (1.0, 0.2)]))
expect_err("時間重複", lambda: custom_ground_motion([(0.0, 0.1), (0.0, 0.2)]))

# ------------------------------------------------------------------
print("=== 層3: nonlinear_seismic_web_analysis() + seismic_to_dict() ===")
f = portal()
pkg = nonlinear_seismic_web_analysis(f, control_node=1, direction='x', zeta=0.05,
                                     ground_motion_type='pulse', pulse_amplitude_g=0.5,
                                     pulse_freq_hz=1.5, pulse_decay=0.25, max_frames=150)
d = json.loads(json.dumps(seismic_to_dict(pkg), allow_nan=False))

f_ref = portal()
ag_ref = sine_pulse_ground_motion(0.5, 1.5, 0.25)
res_ref = seismic_analysis(f_ref, ag_ref, direction='x', zeta=0.05)
assert d['n_steps'] == res_ref.nl.n_steps
assert abs(d['peak_displacement'] - res_ref.peak_displacement(1, 'x')) < 1e-9 * res_ref.peak_displacement(1, 'x')
assert abs(d['period1'] - res_ref.modal.period[0]) < 1e-9
assert abs(d['alpha'] - res_ref.alpha) < 1e-9 and abs(d['beta'] - res_ref.beta) < 1e-9
eb_ref = res_ref.energy_balance()
assert abs(d['energy']['Wext'][-1] - eb_ref['Wext'][-1]) < 1e-6 * abs(eb_ref['Wext'][-1])
dof1x = f_ref.dofs_of(1)[0]
assert float(np.max(np.abs(np.array(d['control_disp']) - res_ref.u[:, dof1x]))) < 1e-9 * res_ref.peak_displacement(1, 'x')
print(f"  跟直接呼叫 seismic_analysis() 逐項相同: n_steps={d['n_steps']}, peak_disp={d['peak_displacement'] * 1000:.3f}mm, "
      f"塑鉸 {[h['label'] for h in d['hinges']]}")
assert {h['label'] for h in d['hinges']} == {'M0 i端', 'M2 i端'}, "應該只列出真的降伏過的塑鉸"
assert all(h['n_yield'] > 0 for h in d['hinges']), "hinges清單裡每一個都必須至少降伏過1次"

# 有效慣性力合力: 對照獨立算出的「總質量 × 地面加速度」
from frame2d.mass import assemble_M as _assemble_M
r_x = np.zeros(res_ref.nl.n_dof)
for nid in f_ref.nodes:
    r_x[f_ref.dofs_of(nid)[0]] = 1.0
total_mass_x = float((_assemble_M(f_ref, 'lumped') @ r_x).sum())
want_eq = [total_mass_x * ag_ref(tt) for tt in res_ref.nl.t]
got_eq = d['equivalent_force']
err_eq = max(abs(a - b) for a, b in zip(got_eq, want_eq)) / max(abs(v) for v in want_eq)
print(f"  有效慣性力合力: 對獨立算出的「總質量×地面加速度」最大相對差 = {err_eq:.2e}")
assert err_eq < 1e-9

assert len(d['frame_idx']) <= 150 and d['frame_idx'][0] == 0 and d['frame_idx'][-1] == d['n_steps']
assert len(d['frames']) == len(d['frame_idx'])
print(f"  抽稀: {d['n_steps'] + 1} 個時間步 -> {len(d['frames'])} 個動畫幀(首尾保留)")

resid = [d['energy']['Wext'][i] - (d['energy']['KE'][i] + d['energy']['Wdamp'][i] + d['energy']['Wint'][i])
        for i in range(len(d['t']))]
scale = max(abs(v) for v in d['energy']['Wext'])
print(f"  JSON化後能量平衡殘差最大值 = {max(abs(v) for v in resid):.3e}(相對尺度 {scale:.3e})")
assert max(abs(v) for v in resid) < 1e-6 * scale

# 不超過 max_frames 的小案例: frame_idx 應該就是全部時間步, 不做抽稀
pkg_small = nonlinear_seismic_web_analysis(portal(), control_node=1, dt=0.02, n_steps=20,
                                           ground_motion_type='pulse', pulse_amplitude_g=0.1, pulse_freq_hz=2.0,
                                           max_frames=400)
d_small = seismic_to_dict(pkg_small)
assert len(d_small['frames']) == 21, "步數不多時不需要抽稀, 應該是全部時間步"
print(f"  小案例(21步): frame_idx 長度 = {len(d_small['frames'])}(不抽稀)")

# ------------------------------------------------------------------
print("=== 層4: 明確拒絕 ===")
expect_err("ground_motion_type打錯", lambda: nonlinear_seismic_web_analysis(portal(), 1, ground_motion_type='bogus'))
expect_err("pulse缺參數", lambda: nonlinear_seismic_web_analysis(portal(), 1, ground_motion_type='pulse'))
expect_err("custom缺點", lambda: nonlinear_seismic_web_analysis(portal(), 1, ground_motion_type='custom'))
expect_err("超過步數上限", lambda: nonlinear_seismic_web_analysis(
    portal(), 1, ground_motion_type='pulse', pulse_amplitude_g=0.5, pulse_freq_hz=1.5, dt=1e-5, n_steps=5000))
expect_err("找不到控制節點", lambda: nonlinear_seismic_web_analysis(
    portal(), 99, ground_motion_type='pulse', pulse_amplitude_g=0.5, pulse_freq_hz=1.5))

print("\n全部通過: sine_pulse/custom地震歷程、nonlinear_seismic_web_analysis與seismic_to_dict跟核心逐項一致。")
