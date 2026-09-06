"""
驗證案例: frame2d/hinge.py -- 塑性鉸狀態機 + 含鉸樑元素勁度矩陣

案例A: 獨立用sympy從「樑ODE + 彈簧力矩-轉角關係聯立、靜力凝聚消去內部
       轉角自由度」重新推導一次凝聚公式, 通過三項極限自我檢核(R->無限大
       退化回標準剛接矩陣、R->0退化回完全機構零矩陣、單端鉸接退化回
       教科書3EI/L公式), 再跟hinge_bending_stiffness()逐項數值比對
       (20組隨機EI/L/R1/R2)。這是移植過來的公式, 不是自己憑空寫的,
       所以獨立重新推導一次比對, 不只是回歸測試自己。
案例B: 彈性狀態(未降伏)應該退化回標準4x4剛接撓曲矩陣, 這是
       hinge_bending_stiffness()跟elements.py既有的member_stiffness_local
       應該吻合的地方(用RIGID_FACTOR近似無限剛接, 允許小量數值誤差)。
案例C: HingeState.check_yield()的降伏判斷邏輯 -- 單向轉換(只能彈性->
       降伏, 不能降伏後又變回彈性), 兩端各自獨立判斷。
案例D: HingeState.performance_level()的IO/LS/CP分類 -- 沒給門檻值時
       一律回傳None; 給了門檻值時, 累積塑性轉角越過各門檻應該正確分類;
       未降伏的鉸永遠是'elastic', 不受門檻值影響。
案例E: full_6x6_with_hinge() -- 軸向自由度不受鉸/Kg影響、P=0時不含
       幾何剛度、跟elements.py的local_geometric_stiffness()疊加時
       數值正確(重用已經驗證過的Kg, 不重新推導)。
"""
import numpy as np
import sympy as sp
from frame2d.hinge import HingeState, hinge_bending_stiffness, full_6x6_with_hinge, RIGID_FACTOR
from frame2d.elements import member_stiffness_local, local_geometric_stiffness

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


# ---- 案例A: 獨立sympy重新推導 + 三項極限自我檢核 + 數值交叉比對 ----
print("=== 案例A: 獨立sympy重新推導凝聚公式 + 三項極限自我檢核 ===")
EI_s, L_s, R1_s, R2_s = sp.symbols('EI L R1 R2', positive=True)
v1, phi1, v2, phi2 = sp.symbols('v1 phi1 v2 phi2')
th1, th2 = sp.symbols('theta1 theta2')

Kb = (EI_s / L_s**3) * sp.Matrix([
    [12,      6 * L_s,      -12,      6 * L_s],
    [6 * L_s, 4 * L_s**2,   -6 * L_s, 2 * L_s**2],
    [-12,    -6 * L_s,       12,     -6 * L_s],
    [6 * L_s, 2 * L_s**2,   -6 * L_s, 4 * L_s**2],
])
disp_internal = sp.Matrix([v1, phi1, v2, phi2])
forces = Kb * disp_internal
V1, M1_beam, V2, M2_beam = forces[0], forces[1], forces[2], forces[3]
eq1 = sp.Eq(M1_beam, R1_s * (th1 - phi1))
eq2 = sp.Eq(M2_beam, R2_s * (th2 - phi2))
sol = sp.solve([eq1, eq2], [phi1, phi2])
phi1_sol, phi2_sol = sp.simplify(sol[phi1]), sp.simplify(sol[phi2])
subs_map = {phi1: phi1_sol, phi2: phi2_sol}
V1_ext, V2_ext = sp.simplify(V1.subs(subs_map)), sp.simplify(V2.subs(subs_map))
M1_ext, M2_ext = sp.simplify(R1_s * (th1 - phi1_sol)), sp.simplify(R2_s * (th2 - phi2_sol))
F_ext = sp.Matrix([V1_ext, M1_ext, V2_ext, M2_ext])
u_ext = sp.Matrix([v1, th1, v2, th2])
K_condensed = sp.simplify(F_ext.jacobian(u_ext))

K_rigid_limit = sp.simplify(K_condensed.limit(R1_s, sp.oo).limit(R2_s, sp.oo))
assert sp.simplify(K_rigid_limit - Kb) == sp.zeros(4, 4), "R->oo應退化回標準剛接矩陣"
K_pinned_limit = sp.simplify(K_condensed.limit(R1_s, 0).limit(R2_s, 0))
assert K_pinned_limit == sp.zeros(4, 4), "R->0應退化回完全機構(零矩陣)"
K_one_pinned = sp.simplify(K_condensed.limit(R2_s, 0).limit(R1_s, sp.oo))
K_textbook = (EI_s / L_s**3) * sp.Matrix([
    [3, 3 * L_s, -3, 0], [3 * L_s, 3 * L_s**2, -3 * L_s, 0],
    [-3, -3 * L_s, 3, 0], [0, 0, 0, 0],
])
assert sp.simplify(K_one_pinned - K_textbook) == sp.zeros(4, 4), "單端鉸接應退化回教科書3EI/L公式"
print("PASS: 三項極限自我檢核全部通過")

K_condensed_func = sp.lambdify((EI_s, L_s, R1_s, R2_s), K_condensed, 'numpy')
rng = np.random.default_rng(0)
for _ in range(20):
    ei = rng.uniform(1e5, 1e8)
    l_ = rng.uniform(1.0, 10.0)
    r1 = rng.uniform(1e3, 1e12)
    r2 = rng.uniform(1e3, 1e12)
    hs = HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=r1, R_post_yield_2=r2)
    hs.yielded = [True, True]   # 強制用R_post_yield(=r1,r2), 不是RIGID_FACTOR
    got = hinge_bending_stiffness(1.0, ei, l_, hs)   # E=1, I=ei => EI=ei
    want = np.array(K_condensed_func(ei, l_, r1, r2), dtype=float)
    assert np.allclose(got, want, rtol=1e-9), f"公式不符(ei={ei},l={l_},r1={r1},r2={r2})"
print("PASS: 20組隨機EI/L/R1/R2下, 封閉式公式跟獨立sympy推導逐項吻合\n")


# ---- 案例B: 彈性狀態應退化回標準剛接矩陣 ----
print("=== 案例B: 彈性(未降伏)狀態應退化回標準剛接矩陣 ===")
hs_elastic = HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1e3, R_post_yield_2=1e3)
kb_hinge = hinge_bending_stiffness(E, I, L, hs_elastic)
kb_standard = member_stiffness_local(E, I, A, L)[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])]
assert np.allclose(kb_hinge, kb_standard, rtol=1e-6), (
    f"未降伏時應該(在RIGID_FACTOR夠大的前提下)近似標準剛接矩陣\n"
    f"hinge=\n{kb_hinge}\nstandard=\n{kb_standard}"
)
print(f"PASS: 未降伏時跟標準剛接矩陣相對誤差 < 1e-6 (RIGID_FACTOR={RIGID_FACTOR:.0e})\n")


# ---- 案例C: check_yield()單向轉換 + 兩端獨立 ----
print("=== 案例C: check_yield()降伏判斷邏輯 ===")
hs = HingeState(Mp1=100.0, Mp2=150.0, R_post_yield_1=1e3, R_post_yield_2=1e3)
assert hs.yielded == [False, False]
assert hs.check_yield(M1=50.0, M2=50.0) is False, "彎矩都還沒到Mp, 不應該降伏"
assert hs.yielded == [False, False]
assert hs.check_yield(M1=100.0, M2=50.0) is True, "端1彎矩達到Mp1, 應該回傳True(新降伏)"
assert hs.yielded == [True, False], "只有端1應該降伏, 端2還沒到Mp2"
assert hs.check_yield(M1=999.0, M2=50.0) is False, "端1已經降伏過, 即使彎矩更大也不算'新'降伏"
assert hs.check_yield(M1=999.0, M2=150.0) is True, "端2現在達到Mp2, 應該回傳True"
assert hs.yielded == [True, True]
print("PASS: 單向轉換 + 兩端獨立判斷都正確\n")


# ---- 案例D: performance_level()分類 ----
print("=== 案例D: performance_level() IO/LS/CP分類 ===")
hs_no_threshold = HingeState(Mp1=100.0, Mp2=100.0, R_post_yield_1=1e3, R_post_yield_2=1e3)
hs_no_threshold.check_yield(M1=100.0, M2=0.0)
assert hs_no_threshold.performance_level(0) is None, "沒給門檻值時應該回傳None, 不是亂猜"
assert hs_no_threshold.performance_level(1) == 'elastic', "端2還沒降伏, 不管有沒有門檻值都是elastic"

hs_labeled = HingeState(Mp1=100.0, Mp2=100.0, R_post_yield_1=1e3, R_post_yield_2=1e3,
                         theta_IO=(0.005, None), theta_LS=(0.015, None), theta_CP=(0.025, None))
assert hs_labeled.performance_level(0) == 'elastic', "還沒降伏前一律elastic, 不看門檻值"
hs_labeled.check_yield(M1=100.0, M2=0.0)
hs_labeled.theta_p[0] = 0.002
assert hs_labeled.performance_level(0) == 'IO', "剛降伏、塑性轉角小於IO門檻"
hs_labeled.theta_p[0] = 0.008
assert hs_labeled.performance_level(0) == 'LS', "超過IO門檻, 還沒到LS門檻"
hs_labeled.theta_p[0] = 0.018
assert hs_labeled.performance_level(0) == 'CP', "超過LS門檻, 還沒到CP門檻"
hs_labeled.theta_p[0] = 0.030
assert hs_labeled.performance_level(0) == 'exceeds CP', "超過CP門檻"
print("PASS: IO/LS/CP分類在各個門檻區間都正確, 且純粹是事後標籤(不影響求解)\n")


# ---- 案例E: full_6x6_with_hinge() ----
print("=== 案例E: full_6x6_with_hinge() 軸向/P-Delta疊加正確性 ===")
hs = HingeState(Mp1=1e30, Mp2=1e30, R_post_yield_1=1e3, R_post_yield_2=1e3)
k_full_p0 = full_6x6_with_hinge(E, A, I, L, hs, P=0.0)
EA_L = E * A / L
assert np.isclose(k_full_p0[0, 0], EA_L) and np.isclose(k_full_p0[3, 3], EA_L)
assert np.isclose(k_full_p0[0, 3], -EA_L)
assert np.allclose(k_full_p0[np.ix_([0, 3], [1, 2, 4, 5])], 0.0), "軸向跟彎曲自由度不應該耦合"

P_test = 500.0
k_full_p = full_6x6_with_hinge(E, A, I, L, hs, P=P_test)
kg_expected = local_geometric_stiffness(P_test, L)
bend_idx = [1, 2, 4, 5]
diff = (k_full_p - k_full_p0)[np.ix_(bend_idx, bend_idx)]
assert np.allclose(diff, kg_expected[np.ix_(bend_idx, bend_idx)], rtol=1e-9), (
    "P非零時應該恰好多疊加elements.local_geometric_stiffness()算出的Kg, 不多不少"
)
print("PASS: 軸向/彎曲不耦合, P-Delta疊加量跟elements.py已驗證過的Kg完全一致\n")

print("PASS: frame2d/hinge.py 所有案例通過")
