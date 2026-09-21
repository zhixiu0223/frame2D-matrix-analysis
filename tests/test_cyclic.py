"""
驗證案例: frame2d.cyclic -- 循環(遲滯)塑鉸與反覆載重 (動力分析路線 D7)

層1 獨立參考實作: 一維雙線性運動硬化「回歸映射」(return mapping, 標準的塑性積分演算法),
    完全不看 cyclic.py 的事件到事件邏輯。帶塑鉸底端的懸臂柱在系統層級就是一個
    (k_e, F_y, k_p) 的雙線性運動硬化彈簧:
        k_e = 3EI/L³,  F_y = Mp/L,  k_p = 1/(L³/(3EI) + L²/R)   (R = 塑鉸硬化剛度)
    對「每一個記錄點」比對力-位移, 涵蓋: 多圈等幅、幅值遞增(背應力累積)、不對稱/部分卸載、
    完全塑性(R=0)、不同步長。
層2 解析頂點: 反向降伏發生在 F = F_a - 2F_y、u = a - 2F_y/k_e(不是 -F_y —— 包辛格型運動
    硬化); 穩態迴圈面積 = 塑性耗能, 用解析頂點的鞋帶公式對照數值。
層3 能量: 完整迴圈內 外力功 ∮F du = 塑鉸累積塑性功。
層4 與既有 pushover 的一致性: 單調載重時, run_cyclic 與 run_pushover(HingeState) 給出同一條
    曲線(門型剛架、兩個塑鉸)。
層5 塑鉸內部記帳: 塑性狀態下 M - R·θp 恆在降伏面上(彎矩記帳與轉角記帳互相獨立算出來, 必須一致)。
層6 多塑鉸門型剛架反覆載重: 迴圈反對稱、穩態迴圈能量守恆、卸載一致性疊代收斂。
層7 明確拒絕。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState, make_protocol, run_cyclic
from frame2d.dofmanager import initial_hinge_states
from frame2d.hinge import HingeState
from frame2d.pushover import run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0
MP, R = 100.0, 2000.0


def column(mp=MP, r=R):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='s', Mp_i=mp, Mp_j=None,
                 R_post_yield_i=r, R_post_yield_j=None)
    f.fix(0)
    return f


KE = 3 * E * I / L**3
FY = MP / L
UY = FY / KE


def system_params(r):
    kp = 1.0 / (1.0 / KE + L**2 / r) if r > 0 else 0.0
    return KE, FY, kp


def return_mapping(ke, fy, kp, u_path):
    """一維雙線性運動硬化(獨立參考): 依序給位移點, 回傳每點的力。"""
    H = ke * kp / (ke - kp) if kp > 0 else 0.0
    F, alpha, out = 0.0, 0.0, [0.0]
    for du in np.diff(u_path):
        Ft = F + ke * du
        f = abs(Ft - alpha) - fy
        if f > 1e-12:
            s = np.sign(Ft - alpha)
            dg = f / (ke + H)
            F = Ft - ke * dg * s
            alpha += H * dg * s
        else:
            F = Ft
        out.append(F)
    return np.array(out)


def run_column(protocol, r=R, step_frac=1 / 3):
    f = column(r=r)
    hs = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
    res = run_cyclic(f, hs, [f.dofs_of(1)[0]], [1.0], protocol, UY * step_frac, [f.dofs_of(0)[0]])
    return res


# ------------------------------------------------------------------
print("=== 層1: 對獨立的回歸映射參考實作, 每個記錄點比對力-位移 ===")
cases = [
    ("等幅 ±6uy × 3圈", make_protocol([6 * UY], 3), R, 1 / 3),
    ("幅值遞增 1.5/3/6uy 各2圈(背應力累積)", make_protocol([1.5 * UY, 3 * UY, 6 * UY], 2), R, 1 / 3),
    ("不對稱與部分卸載", [4 * UY, -1 * UY, 6 * UY, 2 * UY, 5 * UY, -3 * UY, 3 * UY, 0.0], R, 1 / 3),
    ("完全塑性 R=0, ±5uy × 3圈", make_protocol([5 * UY], 3), 0.0, 1 / 3),
    ("步長 uy/17(不同的分步)", make_protocol([6 * UY], 2), R, 1 / 17),
]
for name, proto, r, sf in cases:
    res = run_column(proto, r=r, step_frac=sf)
    ke, fy, kp = system_params(r)
    F_ref = return_mapping(ke, fy, kp, res.u)
    err = np.max(np.abs(res.F - F_ref)) / fy
    print(f"  {name}: {len(res.u)} 個記錄點, 力誤差(相對F_y) 最大 = {err:.2e}")
    assert err < 1e-6, f"{name}: 與獨立回歸映射不符"

# ------------------------------------------------------------------
print("=== 層2: 解析頂點(包辛格型: 反向降伏在 F_a - 2F_y) ===")
a = 6 * UY
ke, fy, kp = system_params(R)
F_a = fy + kp * (a - fy / ke)
u_r = a - 2 * fy / ke
F_r = F_a - 2 * fy
res = run_column(make_protocol([a], 2))
yields = [e for e in res.events if e['kind'] == 'yield']
print(f"  首次降伏: u = {yields[0]['u'] / UY:.4f} uy (解析 1), F = {yields[0]['F'] / fy:.4f} F_y (解析 1)")
print(f"  反向降伏: u = {yields[1]['u'] / UY:.4f} uy (解析 {u_r / UY:.4f}), F = {yields[1]['F'] / fy:.4f} F_y (解析 {F_r / fy:.4f})")
assert abs(yields[0]['u'] - UY) < 1e-6 * UY and abs(yields[0]['F'] - fy) < 1e-6 * fy
assert abs(yields[1]['u'] - u_r) < 1e-6 * UY and abs(yields[1]['F'] - F_r) < 1e-6 * fy, \
    "反向降伏點不在 F_a - 2F_y (不是包辛格型運動硬化)"
assert abs(F_r / fy - (-0.2857)) < 1e-3 and abs(yields[1]['F'] / fy) < 1.0, \
    "反向降伏時的力遠小於 F_y: 這是運動硬化(不是等向硬化, 等向硬化會在 -F_a 才反向降伏)"

# ------------------------------------------------------------------
print("=== 層3: 穩態迴圈面積 = 塑性耗能 (解析頂點鞋帶公式 vs 數值梯形積分 vs 塑鉸累積塑性功) ===")
verts = np.array([[a, F_a], [u_r, F_r], [-a, -F_a], [-u_r, -F_r]])
x, y = verts[:, 0], verts[:, 1]
area_analytic = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
res = run_column(make_protocol([a], 4))
u_arr = res.u
# 找第3圈起點(正向轉折點 u=a 的第3次出現)與第4圈起點
turn_pos = [i for i in range(1, len(u_arr) - 1) if u_arr[i] >= a - 1e-12 and u_arr[i] >= u_arr[i - 1] and u_arr[i] >= u_arr[i + 1]]
assert len(turn_pos) >= 4, f"應該有 4 個正向轉折點, 實際 {len(turn_pos)}"
i0, i1 = turn_pos[2], turn_pos[3]                    # 第3圈: a -> -a -> a
W_ext = res.loop_work(i0, i1)
W_hinge = res.hinge_work[(0, 0)][i1] - res.hinge_work[(0, 0)][i0]
print(f"  第3圈: 解析面積 = {area_analytic:.6f}, 外力功 ∮F du = {W_ext:.6f}, 塑鉸累積塑性功 = {W_hinge:.6f}")
assert abs(W_ext - area_analytic) < 1e-6 * area_analytic, "數值迴圈面積與解析鞋帶面積不符"
assert abs(W_hinge - area_analytic) < 1e-6 * area_analytic, "塑鉸累積塑性功與迴圈面積(耗能)不符"
# 前一圈與這一圈完全相同(運動硬化第2圈起就是穩態)
W2 = res.loop_work(turn_pos[1], turn_pos[2])
assert abs(W2 - W_ext) < 1e-6 * W_ext, "等幅運動硬化: 穩態迴圈應逐圈相同"

# ------------------------------------------------------------------
print("=== 層5: 塑鉸內部記帳 —— 塑性狀態下 M - R·θp 恆在降伏面上 ===")
res = run_column(make_protocol([1.5 * UY, 4 * UY, 6 * UY], 2))
M = res.hinge_M[(0, 0)]
tp = res.hinge_theta_p[(0, 0)]
resid = M - R * tp                                   # 運動硬化: 塑性時 M - α = ±Mp, α = R·θp
sgn = np.sign(resid)
assert np.max(np.abs(resid)) <= MP * (1 + 1e-6), "彎矩超出降伏面(M - R·θp 的絕對值 > Mp)"
# 進入塑性(有降伏事件)之後的點, 塑性階段的點必須在降伏面上
plastic_pts = np.where(np.abs(np.abs(resid) - MP) < 1e-6 * MP)[0]
print(f"  最大 |M - R·θp| = {np.max(np.abs(resid)):.9f} (Mp = {MP}); 落在降伏面上的點 {len(plastic_pts)}/{len(M)}")
assert len(plastic_pts) > 0.3 * len(M), "應該有相當比例的記錄點在塑性(降伏面上)"
assert res.hinge_states[0].theta_p[0] >= abs(res.hinge_states[0].theta_p_signed[0]) - 1e-15, "∑|Δθp| 不能小於 |∑Δθp|"

# ------------------------------------------------------------------
def portal(mp=60.0, r=1500.0):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='s', Mp_i=mp, Mp_j=mp, R_post_yield_i=r, R_post_yield_j=r)
    f.add_member(1, node_i=1, node_j=2, section='s', Mp_i=mp, Mp_j=mp, R_post_yield_i=r, R_post_yield_j=r)
    f.add_member(2, node_i=3, node_j=2, section='s', Mp_i=mp, Mp_j=mp, R_post_yield_i=r, R_post_yield_j=r)
    f.fix(0).fix(3)
    return f


print("=== 層4: 單調載重與既有 run_pushover(HingeState) 一致 (門型剛架, 6 個塑鉸) ===")
f = portal()
dof = f.dofs_of(1)[0]
base = [f.dofs_of(0)[0], f.dofs_of(3)[0]]
target, dn = 0.06, 0.004
u_p, F_p, ev_p, hs_p, mech_p = run_pushover(f, initial_hinge_states(f), [dof], [1.0], target, dn, base)
f2 = portal()
hs_c = CyclicHingeState.from_hinge_states(initial_hinge_states(f2))
res_m = run_cyclic(f2, hs_c, [dof], [1.0], [target], dn, base)
print(f"  pushover: {len(u_p)} 點, 最終 F = {F_p[-1]:.6f}; cyclic: {len(res_m.u)} 點, 最終 F = {res_m.F[-1]:.6f}")
assert len(u_p) == len(res_m.u), "單調載重時兩者記錄的點數應該相同(同樣的事件切分)"
assert np.max(np.abs(np.array(u_p) - res_m.u)) < 1e-12
assert np.max(np.abs(np.array(F_p) - res_m.F)) < 1e-9 * np.max(np.abs(F_p)), "單調載重時 cyclic 必須與 pushover 逐點一致"
n_yield = sum(1 for e in res_m.events if e['kind'] == 'yield')
assert n_yield >= 4, f"這個載重量應該至少形成 4 個塑鉸, 實際 {n_yield}"
print(f"  逐點一致(最大差 {np.max(np.abs(np.array(F_p) - res_m.F)):.2e}), 降伏事件 {n_yield} 個")

# ------------------------------------------------------------------
print("=== 層6: 多塑鉸門型剛架反覆載重 ===")
f = portal()
hs = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
proto = make_protocol([0.02, 0.05], 3)
res = run_cyclic(f, hs, [dof], [1.0], proto, 0.004, base)
n_y = sum(1 for e in res.events if e['kind'] == 'yield')
n_u = sum(1 for e in res.events if e['kind'] == 'unload')
print(f"  {len(res.u)} 點, 降伏事件 {n_y}, 卸載事件 {n_u}")
assert n_y > 6 and n_u > 6, "多圈反覆載重應該有多次降伏與卸載"
# 對稱結構 + 對稱位移歷程: 正向與反向對應點的力互為反號(穩態圈)
u_arr, F_arr = res.u, res.F
i_pos = [i for i in range(1, len(u_arr) - 1) if abs(u_arr[i] - 0.05) < 1e-12 and u_arr[i] >= u_arr[i - 1] and u_arr[i] >= u_arr[i + 1]]
i_neg = [i for i in range(1, len(u_arr) - 1) if abs(u_arr[i] + 0.05) < 1e-12 and u_arr[i] <= u_arr[i - 1] and u_arr[i] <= u_arr[i + 1]]
print(f"  幅值 0.05: 正向峰 F = {F_arr[i_pos[-1]]:.6f}, 反向峰 F = {F_arr[i_neg[-1]]:.6f}")
assert abs(F_arr[i_pos[-1]] + F_arr[i_neg[-1]]) < 1e-6 * abs(F_arr[i_pos[-1]]), "穩態下正反向峰值力應該互為反號"
# 穩態迴圈: 兩個連續 (+a→-a→+a) 迴圈的能量相同, 且外力功 = 全部塑鉸塑性功
ip = i_pos[-3:]
W1, W2 = res.loop_work(ip[0], ip[1]), res.loop_work(ip[1], ip[2])
Wh = [sum(res.hinge_work[k][ip[j + 1]] - res.hinge_work[k][ip[j]] for k in res.hinge_work) for j in range(2)]
print(f"  最後兩圈: 外力功 {W1:.6f}, {W2:.6f}; 塑鉸累積塑性功 {Wh[0]:.6f}, {Wh[1]:.6f}")
assert abs(W1 - W2) < 1e-6 * W1, "等幅穩態迴圈能量應逐圈相同"
assert abs(W2 - Wh[1]) < 1e-6 * W2, "穩態迴圈: 外力功必須等於全部塑鉸的塑性功(能量守恆)"
# 每個塑鉸的記帳: 彎矩不超出降伏面
for (mid, e), Mh in res.hinge_M.items():
    tph = res.hinge_theta_p[(mid, e)]
    assert np.max(np.abs(Mh - 1500.0 * tph)) <= 60.0 * (1 + 1e-5), f"鉸 ({mid},{e}) 彎矩超出降伏面"
print("  6 個塑鉸的 M - R·θp 全部在降伏面內 OK")

# ------------------------------------------------------------------
print("=== 層7: 明確拒絕 ===")
def expect(label, exc, fn):
    try:
        fn()
    except exc as e:
        print(f"  {label}: {exc.__name__} OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該丟出 {exc.__name__}")

f = column()
hs_plain = initial_hinge_states(f)
expect("傳入單向 HingeState", TypeError,
       lambda: run_cyclic(f, hs_plain, [f.dofs_of(1)[0]], [1.0], [0.1], 0.01, [f.dofs_of(0)[0]]))
expect("direction[0]=0", ValueError,
       lambda: run_cyclic(f, CyclicHingeState.from_hinge_states(initial_hinge_states(f)), [f.dofs_of(1)[0]], [0.0], [0.1], 0.01, [f.dofs_of(0)[0]]))
expect("呼叫單向 check_yield", NotImplementedError,
       lambda: CyclicHingeState.from_hinge_states(initial_hinge_states(f))[0].check_yield(1.0, 1.0))
expect("步數安全閥", RuntimeError,
       lambda: run_cyclic(f, CyclicHingeState.from_hinge_states(initial_hinge_states(f)), [f.dofs_of(1)[0]], [1.0],
                          [10.0], 1e-6, [f.dofs_of(0)[0]], max_steps=50))
assert isinstance(HingeState(1, 1, 0, 0), HingeState) and not isinstance(HingeState(1, 1, 0, 0), CyclicHingeState)
print("  既有 HingeState 不受影響(仍是獨立的單向類別)")

print("\n全部通過: 循環塑鉸對獨立回歸映射、解析頂點、能量守恆、既有 pushover 一致性、多塑鉸記帳都吻合。")
