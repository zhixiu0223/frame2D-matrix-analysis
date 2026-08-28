"""
驗證案例13: Phase 2(局部段均佈載重)對照 SW FEA 第三方工具

使用者上傳的 phase-02.frame / phase-02_Report.pdf: 簡支梁 L=12m, 均佈載重
20kN/m 只加在[3,7]這一段, 斷面 E=200GPa, I=8.0E7mm^4(=8e-5 m^4), 跟
tests/test_partial_udl.py 案例B、examples/load_system_v2_demo.py 的
Phase2案例是同一個模型(參數完全一致)。

這次驗證比之前更完整: SW FEA報告給了11個位置點(0.0L~1.0L)的SF/BM/dY逐點
數值, 不是只比反力/彎矩峰值, 是整條曲線逐點核對。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces, member_deformed_shape

E, I, A = 200e6, 8e-5, 1e-2   # 對應SW FEA的 E=200GPa, I=8.0E7mm^4
L, w, c, d = 12.0, -20.0, 3.0, 7.0

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.distributed_load(0, w=w, x_start=c, x_end=d)
r = solve(f)

# ---- 反力 ----
print("=== 反力比對 ===")
Rx0, Ry0, _ = f.dofs_of(0)
Rx1, Ry1, _ = f.dofs_of(1)
print(f"  node0: Ry={r.reactions[Ry0]:.4f}  (SW FEA: 46.667)")
print(f"  node1: Ry={r.reactions[Ry1]:.4f}  (SW FEA: 33.333)")
assert abs(r.reactions[Ry0] - 46.667) < 1e-2
assert abs(r.reactions[Ry1] - 33.333) < 1e-2
print("PASS\n")

# ---- SW FEA報告的11個位置點資料 (Location, SF, BM, dY[mm]) ----
sw_data = [
    (0.0, 46.667, -0.000, 0.000), (0.1, 46.667, 56.000, -52.244), (0.2, 46.667, 112.000, -99.447),
    (0.3, 34.667, 164.400, -136.578), (0.4, 10.667, 191.600, -159.121), (0.5, -13.333, 190.000, -164.637),
    (0.6, -33.333, 160.000, -153.269), (0.7, -33.333, 120.000, -127.552), (0.8, -33.333, 80.000, -91.036),
    (0.9, -33.333, 40.000, -47.320), (1.0, -33.333, 0.000, -0.003),
]

x_query = np.array([frac * L for frac, _, _, _ in sw_data])
x_full, N_full, V_full, M_full = member_internal_forces(f, r, 0, n=401)
V_interp = np.interp(x_query, x_full, V_full)
M_interp = np.interp(x_query, x_full, M_full)
X_def, Y_def = member_deformed_shape(f, r, 0, scale=1.0, n=401)
dY_interp = np.interp(x_query, X_def, Y_def) * 1000   # m -> mm

print("=== 逐點比對 (11個位置, SF/BM/dY) ===")
print(f"{'x/L':>6}{'SF(SW)':>10}{'SF(f2d)':>10}{'BM(SW)':>10}{'BM(f2d)':>10}{'dY(SW)':>12}{'dY(f2d)':>12}")
max_err_sf, max_err_bm, max_err_dy = 0, 0, 0
for (frac, sf, bm, dy), v, m, dyv in zip(sw_data, V_interp, M_interp, dY_interp):
    print(f"{frac:>6.1f}{sf:>10.3f}{v:>10.3f}{bm:>10.3f}{m:>10.3f}{dy:>12.3f}{dyv:>12.3f}")
    max_err_sf = max(max_err_sf, abs(v - sf))
    max_err_bm = max(max_err_bm, abs(m - bm))
    max_err_dy = max(max_err_dy, abs(dyv - dy))

print(f"\n最大誤差: SF={max_err_sf:.3f}kN  BM={max_err_bm:.3f}kNm  dY={max_err_dy:.3f}mm")
assert max_err_sf < 0.01, "剪力誤差過大"
assert max_err_bm < 0.01, "彎矩誤差過大"
assert max_err_dy < 0.01, "撓度誤差過大"
print("PASS: 11個位置點的SF/BM/dY全數跟SW FEA吻合(誤差<0.01)\n")

print("PASS: Phase 2局部段均佈載重, 對照SW FEA逐點驗證完成")
