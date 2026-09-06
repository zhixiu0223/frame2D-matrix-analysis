"""
驗證案例: 塑鉸容量接進Frame2D模型 + dofmanager.solve_with_hinges()

案例A: Member的Mp_i/Mp_j/R_post_yield_i/R_post_yield_j欄位預設None時,
       跟完全沒有這些欄位的舊版行為逐位元一致(add_member()不傳這幾個
       參數時)。
案例B: 設定了塑鉸容量但目前(全部未降伏)的模型, solve_with_hinges()
       應該非常接近純彈性solve_dofmanager()的結果(差異只來自
       RIGID_FACTOR是很大但不是無限大的數值近似)。
案例C: initial_hinge_states()只掃兩端Mp都是None的member會被跳過, 只要
       任一端設定了就會出現在回傳字典裡。
案例D: 手動把某一端設成已降伏狀態, solve_with_hinges()算出來的側向
       位移/轉角, 跟獨立用hinge_bending_stiffness()(已在
       test_hinge_condensation.py驗證過的公式)組出縮減後2自由度系統
       手解的結果逐位吻合——這裡驗證的是dofmanager.py"有沒有正確接上
       這個公式並套用邊界條件", 不是重新驗證公式本身(那是
       test_hinge_condensation.py的責任)。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.dofmanager import solve_dofmanager, solve_with_hinges, initial_hinge_states
from frame2d.hinge import HingeState, hinge_bending_stiffness
from frame2d.pushover import run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


def cantilever(with_hinge_capacity=False):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('sec', E=E, I=I, A=A)
    if with_hinge_capacity:
        f.add_member(0, node_i=0, node_j=1, section='sec',
                     Mp_i=200.0, Mp_j=200.0, R_post_yield_i=1e4, R_post_yield_j=1e4)
    else:
        f.add_member(0, node_i=0, node_j=1, section='sec')
    f.fix(0)
    f.point_load(1, fy=-10.0)
    return f


# ---- 案例A: 不傳塑鉸參數時, Member欄位預設None不影響任何行為 ----
print("=== 案例A: 不設塑鉸容量時欄位預設None, 行為不變 ===")
f = cantilever(with_hinge_capacity=False)
m = f.members[0]
assert m.Mp_i is None and m.Mp_j is None
assert m.R_post_yield_i is None and m.R_post_yield_j is None
assert initial_hinge_states(f) == {}, "沒有member設定Mp_i時, initial_hinge_states()應該回傳空字典"
print("PASS: 預設None、initial_hinge_states()正確回傳空字典\n")


# ---- 案例B: 設定容量但全部未降伏時, 應該非常接近純彈性解 ----
print("=== 案例B: 未降伏時solve_with_hinges()應接近solve_dofmanager() ===")
f_hinge = cantilever(with_hinge_capacity=True)
f_elastic = cantilever(with_hinge_capacity=False)
r_hinge = solve_with_hinges(f_hinge)             # hinge_states=None -> 自動建立、全部未降伏
r_elastic = solve_dofmanager(f_elastic)
v_hinge = r_hinge.displacements[f_hinge.dofs_of(1)[1]]
v_elastic = r_elastic.displacements[f_elastic.dofs_of(1)[1]]
rel_diff = abs(v_hinge - v_elastic) / abs(v_elastic)
assert rel_diff < 1e-6, f"未降伏時側向位移相對差異應該極小, 實際={rel_diff:.2e}"
print(f"PASS: 未降伏時側向位移相對差異 {rel_diff:.2e} (RIGID_FACTOR近似, 應遠小於1e-6)\n")


# ---- 案例C: initial_hinge_states()只掃有設Mp_i的member ----
print("=== 案例C: initial_hinge_states()只掃有設Mp_i的member ===")
f_mixed = Frame2D()
f_mixed.add_node(0, 0, 0); f_mixed.add_node(1, L, 0); f_mixed.add_node(2, 2 * L, 0)
f_mixed.add_section('sec', E=E, I=I, A=A)
f_mixed.add_member(0, node_i=0, node_j=1, section='sec',
                    Mp_i=200.0, Mp_j=200.0, R_post_yield_i=1e4, R_post_yield_j=1e4)
f_mixed.add_member(1, node_i=1, node_j=2, section='sec')   # 沒設塑鉸容量
states = initial_hinge_states(f_mixed)
assert set(states.keys()) == {0}, f"只有member 0設了Mp_i, 應該只出現member 0, 實際={set(states.keys())}"
assert states[0].yielded == [False, False]
print("PASS: 只掃到member 0, member 1(沒設容量)正確被忽略\n")


# ---- 案例D: 手動設成已降伏, 交叉驗證位移數值(驗證wiring, 不重驗公式) ----
print("=== 案例D: 已降伏狀態下, solve_with_hinges()跟獨立手解縮減系統交叉驗證 ===")
f_yielded = cantilever(with_hinge_capacity=True)
R1_val = 5e4
hs = HingeState(Mp1=200.0, Mp2=200.0, R_post_yield_1=R1_val, R_post_yield_2=1e4)
hs.yielded = [True, False]   # 端1(在固定端node0這一側)已降伏, 端2(自由端)還是彈性
hinge_states = {0: hs}

r = solve_with_hinges(f_yielded, hinge_states=hinge_states)
v2_fem = r.displacements[f_yielded.dofs_of(1)[1]]
th2_fem = r.displacements[f_yielded.dofs_of(1)[2]]

# 獨立解: 直接用已驗證過的hinge_bending_stiffness()組4x4矩陣, 取自由度
# (v2,theta2)那個2x2子矩陣(node0的v1,theta1被fix()設成0, 對應Dirichlet
# 邊界條件消去, 只要沒有耦合項要處理的自由端就是取子矩陣), 手解
Kb_hinge = hinge_bending_stiffness(E, I, L, hs)   # 對(v1,theta1,v2,theta2)
K_free = Kb_hinge[np.ix_([2, 3], [2, 3])]
rhs = np.array([-10.0, 0.0])   # fy=-10.0, 無外加彎矩
sol = np.linalg.solve(K_free, rhs)
v2_hand, th2_hand = sol[0], sol[1]

print(f"v2: dofmanager={v2_fem:.10e}  獨立手解={v2_hand:.10e}")
print(f"th2: dofmanager={th2_fem:.10e}  獨立手解={th2_hand:.10e}")
assert abs((v2_fem - v2_hand) / v2_hand) < 1e-9, "側向位移跟獨立手解結果不符"
assert abs((th2_fem - th2_hand) / th2_hand) < 1e-9, "端點轉角跟獨立手解結果不符"

# 順便確認: 降伏後應該比完全彈性時位移更大(勁度變小)
v2_elastic_full = solve_dofmanager(cantilever(with_hinge_capacity=False)).displacements[
    f_yielded.dofs_of(1)[1]]
assert abs(v2_fem) > abs(v2_elastic_full), "降伏後側向勁度變小, 位移應該比純彈性時更大"
print(f"PASS: 降伏後位移跟dofmanager()組裝結果吻合, 且比純彈性(v2={v2_elastic_full:.6e})更大\n")


# ---- 案例E: 只填一端Mp、另一端留空None, 應該正確處理成"那端永遠彈性"
#      (2026-09修正: 這裡原本有兩個bug, 見initial_hinge_states()的
#      docstring說明; 這裡補上這兩個具體情境的回歸測試, 避免以後改壞) ----
print("=== 案例E: 只填一端Mp, 另一端留空應正確視為\"永遠彈性\" ===")

f_only_j = Frame2D()
f_only_j.add_node(0, 0, 0); f_only_j.add_node(1, 0, L)
f_only_j.add_section('sec', E=E, I=I, A=A)
f_only_j.add_member(0, node_i=0, node_j=1, section='sec',
                     Mp_j=200.0, R_post_yield_j=1e4)   # i端(Mp_i)故意留空
f_only_j.fix(0)
states_only_j = initial_hinge_states(f_only_j)
assert 0 in states_only_j, "只填j端Mp時, 這根桿件應該還是要出現在hinge_states裡(之前的bug會整根被忽略)"
assert states_only_j[0].Mp[0] == float('inf'), "留空的i端應該轉換成Mp=inf(永遠不會降伏)"
assert states_only_j[0].Mp[1] == 200.0
print("PASS: 只填j端Mp, 桿件正確出現在hinge_states裡, 留空的i端Mp=inf")

f_only_i = Frame2D()
f_only_i.add_node(0, 0, 0); f_only_i.add_node(1, 0, L)
f_only_i.add_section('sec', E=E, I=I, A=A)
f_only_i.add_member(0, node_i=0, node_j=1, section='sec',
                     Mp_i=100.0, R_post_yield_i=1e4)   # j端(Mp_j)故意留空
f_only_i.fix(0)
states_only_i = initial_hinge_states(f_only_i)
u_e, F_e, ev_e, hs_e, mech_e = run_pushover(
    f_only_i, states_only_i, prescribed_dofs=[f_only_i.dofs_of(1)[0]], direction=[1.0],
    target_total=0.05, d_nominal=0.001, base_reaction_dofs=[f_only_i.dofs_of(0)[0]],
)
assert hs_e[0].yielded == [True, False], "i端應該正常降伏, j端(留空=Mp=inf)永遠不該降伏"
print("PASS: 只填i端Mp, run_pushover()不再crash, i端正常降伏、j端(留空)永遠不降伏\n")

print("PASS: 塑鉸容量wiring所有案例通過")
