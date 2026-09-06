"""
驗證案例: 多鉸情境下的側移機構偵測(check_mechanism())

單一懸臂柱(test_pushover.py)只測了一個鉸, 沒有真正的"機構"可以形成
(降伏後still有正的殘餘勁度)。這裡用一個對稱單跨門型鋼架(兩根柱、一根
樑), 柱兩端都設塑鉸(弱柱)、樑保持彈性(強樑), 降伏後硬化剛度設得極小
(近乎完全塑性), 同時側推兩根柱頂——這是classic的"柱側移機構"案例:
四個鉸(兩柱各兩端)全部降伏後, 結構在側推方向真的推不動了。

案例A: 彈性側向勁度跟獨立手動組裝(直接呼叫elements.py已驗證過的
       member_stiffness_local()跟transformation_matrix(), 不透過
       pushover.py的組裝路徑)算出來的12自由度縮減系統逐位吻合。
案例B: 對稱結構下, 應該恰好2個依序事件, 每個事件都同時觸發兩根柱
       "對稱的那一端"(例如都是柱底, 或都是柱頂)——因為結構跟載重都
       對稱, 不應該只有一根柱先降伏。
案例C: 第一個事件發生後(只有2/4鉸降伏), mechanism_reached還不應該是
       True——結構還沒真的變成機構, 不能太早誤判。
案例D: 全部4個鉸都降伏後, mechanism_reached應該變True, 且兩根柱的
       hinge_state.yielded都是[True, True]。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.elements import member_stiffness_local, transformation_matrix
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover

E = 200e6
I_col, A_col = 8e-5, 1e-2
I_beam, A_beam = 8e-5, 1e-2
H, Lspan = 4.0, 6.0
Mp = 100.0
R_post_yield = 1.0   # 近乎完全塑性(降伏後幾乎沒有硬化), 機構才推得出來


def build_portal_frame():
    f = Frame2D()
    f.add_node(0, 0, 0); f.add_node(1, Lspan, 0)
    f.add_node(2, 0, H); f.add_node(3, Lspan, H)
    f.add_section('col', E=E, I=I_col, A=A_col)
    f.add_section('beam', E=E, I=I_beam, A=A_beam)
    f.add_member(0, node_i=0, node_j=2, section='col',
                 Mp_i=Mp, Mp_j=Mp, R_post_yield_i=R_post_yield, R_post_yield_j=R_post_yield)
    f.add_member(1, node_i=1, node_j=3, section='col',
                 Mp_i=Mp, Mp_j=Mp, R_post_yield_i=R_post_yield, R_post_yield_j=R_post_yield)
    f.add_member(2, node_i=2, node_j=3, section='beam')   # 樑保持彈性, 不設塑鉸
    f.fix(0); f.fix(1)
    return f


# ---- 案例A: 彈性側向勁度跟獨立手動組裝(繞過pushover.py)逐位吻合 ----
print("=== 案例A: 彈性側向勁度跟獨立手動組裝12自由度系統交叉驗證 ===")
nodes = {0: (0, 0), 1: (Lspan, 0), 2: (0, H), 3: (Lspan, H)}
n = 12
K_manual = np.zeros((n, n))


def add_member_manual(ni, nj, E_, I_, A_):
    xi, yi = nodes[ni]; xj, yj = nodes[nj]
    L = np.hypot(xj - xi, yj - yi)
    angle = np.arctan2(yj - yi, xj - xi)
    k_local = member_stiffness_local(E_, I_, A_, L)
    T = transformation_matrix(angle)
    k_global = T.T @ k_local @ T
    dofs = [3 * ni, 3 * ni + 1, 3 * ni + 2, 3 * nj, 3 * nj + 1, 3 * nj + 2]
    K_manual[np.ix_(dofs, dofs)] += k_global


add_member_manual(0, 2, E, I_col, A_col)
add_member_manual(1, 3, E, I_col, A_col)
add_member_manual(2, 3, E, I_beam, A_beam)

fixed = [0, 1, 2, 3, 4, 5]
prescribed = [6, 9]
free = sorted(set(range(n)) - set(fixed) - set(prescribed))
du_p = np.array([1.0, 1.0])
K_ff = K_manual[np.ix_(free, free)]
K_fp = K_manual[np.ix_(free, prescribed)]
du_free = np.linalg.solve(K_ff, -K_fp @ du_p)
du_full = np.zeros(n)
for d, v in zip(prescribed, du_p):
    du_full[d] = v
for d, v in zip(free, du_free):
    du_full[d] = v
reaction = K_manual @ du_full
K_elastic_manual = -sum(reaction[d] for d in (0, 3))   # 底剪力 = node0,node1的ux反力加總

f = build_portal_frame()
hinge_states = initial_hinge_states(f)
control_dofs = [f.dofs_of(2)[0], f.dofs_of(3)[0]]
base_dofs = [f.dofs_of(0)[0], f.dofs_of(1)[0]]
u_hist, F_hist, event_log, hs, mechanism = run_pushover(
    f, hinge_states, prescribed_dofs=control_dofs, direction=[1.0, 1.0],
    target_total=0.2, d_nominal=0.0005, base_reaction_dofs=base_dofs,
    mechanism_ratio_limit=1e-6,
)
K_elastic_numeric = (F_hist[1] - F_hist[0]) / (u_hist[1] - u_hist[0])
assert np.isclose(K_elastic_numeric, K_elastic_manual, rtol=1e-6), (
    f"彈性側向勁度應該跟獨立手動組裝的{K_elastic_manual:.4f}吻合, 實際={K_elastic_numeric:.4f}"
)
print(f"PASS: 彈性側向勁度 {K_elastic_numeric:.4f} 跟獨立手動組裝 {K_elastic_manual:.4f} 吻合\n")


# ---- 案例B: 對稱結構應該恰好2個依序事件, 每個都同時觸發兩根柱 ----
print("=== 案例B: 對稱結構應該有2個事件, 每個事件同時觸發兩根柱同一端 ===")
assert len(event_log) == 2, f"對稱結構應該恰好2個依序事件(先兩柱底, 後兩柱頂), 實際={len(event_log)}"
for i, ev in enumerate(event_log):
    assert len(ev['yielded']) == 2, f"事件{i}應該同時觸發兩根柱(對稱), 實際={ev['yielded']}"
    end_indices = {end_idx for _, end_idx in ev['yielded']}
    assert len(end_indices) == 1, f"事件{i}兩根柱應該是同一端(對稱), 實際={ev['yielded']}"
print(f"PASS: 恰好2個事件, 每個事件都同時觸發兩根柱的同一端 "
      f"(事件1: {event_log[0]['yielded']}, 事件2: {event_log[1]['yielded']})\n")


# ---- 案例C: 只有2/4鉸降伏時, 還不應該判定為機構 ----
print("=== 案例C: 只有一半的鉸降伏時, 不應該誤判成機構 ===")
f_partial = build_portal_frame()
hinge_states_partial = initial_hinge_states(f_partial)
u_p, F_p, ev_p, hs_p, mech_p = run_pushover(
    f_partial, hinge_states_partial, prescribed_dofs=[f_partial.dofs_of(2)[0], f_partial.dofs_of(3)[0]],
    direction=[1.0, 1.0], target_total=event_log[0]['u'] * 1.15, d_nominal=0.0005,
    base_reaction_dofs=[f_partial.dofs_of(0)[0], f_partial.dofs_of(1)[0]],
    mechanism_ratio_limit=1e-6,
)
assert len(ev_p) == 1, f"這個target_total應該只推到第一個事件之後一點點, 實際事件數={len(ev_p)}"
assert mech_p is False, "只有柱底那一端降伏(2/4鉸), 柱頂還是彈性, 不應該判定為機構"
print("PASS: 只有2/4鉸降伏時, mechanism_reached正確維持False\n")


# ---- 案例D: 全部4個鉸都降伏後, 應該正確判定為機構 ----
print("=== 案例D: 全部4個鉸都降伏後, 應該正確判定為機構 ===")
assert mechanism is True, "四個鉸(兩柱各兩端)全部降伏後, 應該判定為側移機構"
assert hs[0].yielded == [True, True], f"柱0應該兩端都降伏, 實際={hs[0].yielded}"
assert hs[1].yielded == [True, True], f"柱1應該兩端都降伏, 實際={hs[1].yielded}"
print("PASS: 全部4個鉸降伏後, mechanism_reached=True, 兩柱yielded都是[True, True]\n")

print("PASS: 門型鋼架多鉸機構偵測所有案例通過")
