"""
驗證案例: run_pushover(..., include_snapshots=True) -- 逐步回放用的快照

案例A: 預設(include_snapshots=False)完全不影響既有5個回傳值的呼叫方式。
案例B: 快照數量跟history_u/history_F逐一對齊, 第0筆是u=0的初始狀態
       (全部未降伏、彎矩0)。
案例C: 降伏事件發生的那一步, 快照裡的M1應該恰好等於Mp(這是event-to-event
       演算法的定義: 精確走到M剛好等於Mp的那個比例, 不是隨便一個時間點)。
案例D: 快照是深拷貝, 不會被後續的HingeState物件原地修改牽動——拿走
       snapshot[0]之後, 即使hinge_states之後繼續被run_pushover原地
       修改, snapshot[0]裡的數字應該維持u=0時的原始值不變。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0
Mp_base = 100.0


def cantilever_column():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec',
                 Mp_i=Mp_base, Mp_j=1e30, R_post_yield_i=100.0, R_post_yield_j=100.0)
    f.fix(0)
    return f


# ---- 案例A: 預設不影響既有呼叫方式 ----
print("=== 案例A: 預設呼叫方式不受影響 ===")
f = cantilever_column()
hs = initial_hinge_states(f)
result = run_pushover(f, hs, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
                       target_total=0.05, d_nominal=0.01, base_reaction_dofs=[f.dofs_of(0)[0]])
assert len(result) == 5, f"預設呼叫應該回傳5個值(跟既有測試一致), 實際={len(result)}"
print("PASS: 預設呼叫回傳5個值, 沒有被新參數影響\n")


# ---- 案例B: 快照數量對齊 + 初始狀態正確 ----
print("=== 案例B: 快照數量跟history_u對齊, 初始狀態正確 ===")
f2 = cantilever_column()
hs2 = initial_hinge_states(f2)
u, F, event_log, hs_final, mech, snapshots = run_pushover(
    f2, hs2, prescribed_dofs=[f2.dofs_of(1)[0]], direction=[1.0],
    target_total=0.05, d_nominal=0.01, base_reaction_dofs=[f2.dofs_of(0)[0]],
    include_snapshots=True,
)
assert len(snapshots) == len(u), f"快照數量應該跟history_u逐一對齊, {len(snapshots)} vs {len(u)}"
assert snapshots[0]['hinge_states'][0]['yielded'] == [False, False]
assert snapshots[0]['member_forces'][0] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
print(f"PASS: {len(snapshots)}筆快照跟history_u逐一對齊, 初始狀態(u=0)正確是全0/未降伏\n")


# ---- 案例C: 降伏事件那一步的M1應該恰好等於Mp ----
print("=== 案例C: 降伏事件那一步的快照M1應該恰好等於Mp ===")
assert len(event_log) == 1
event_u = event_log[0]['u']
event_step_idx = np.where(np.isclose(u, event_u))[0]
assert len(event_step_idx) == 1, "應該剛好對應到一筆快照"
M1_at_event = snapshots[event_step_idx[0]]['member_forces'][0][2]   # 索引2 = M1
assert abs(abs(M1_at_event) - Mp_base) < 1e-6, (
    f"降伏事件那一步的M1絕對值應該恰好等於Mp={Mp_base}, 實際={M1_at_event}"
)
print(f"PASS: 降伏事件步驟(u={event_u:.6f})的快照M1={M1_at_event:.6f}, 恰好等於Mp={Mp_base}\n")


# ---- 案例D: 快照是深拷貝, 不受後續原地修改牽動 ----
print("=== 案例D: 快照深拷貝安全性 ===")
snapshot0_before = dict(snapshots[0]['hinge_states'][0])
hs_final[0].yielded[0] = False   # 手動亂改, 模擬"如果快照沒有深拷貝會被牽連"的情況
hs_final[0].theta_p[0] = 999.0
assert snapshots[0]['hinge_states'][0]['yielded'] == [False, False], (
    "快照0(u=0時)不應該被後續對hinge_states的修改影響"
)
assert snapshots[0]['hinge_states'][0]['theta_p'] == [0.0, 0.0], (
    "快照0的theta_p不應該被後續修改影響"
)
print("PASS: 快照是深拷貝, 後續修改hinge_states不會牽動已經記錄的舊快照\n")

print("PASS: 逐步回放快照所有案例通過")
