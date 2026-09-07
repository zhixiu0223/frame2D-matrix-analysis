"""
驗證案例: run_pushover(..., control_mode='force') -- 力控制

案例A: 純彈性範圍內(不設塑鉸容量, 或塑鉸容量設很大不會降伏), 力控制
       解出來的位移應該恰好等於標準懸臂梁公式 F*L^3/(3EI)——這是力控制
       存在的主要理由: 直接施加已知的力, 解出轉角/位移, 拿去跟手算比對,
       不用像位移控制那樣先猜位移。
案例B: 力控制下history_F應該恰好等於施加的力(不是像位移控制那樣要
       解出來的量, 力控制下"力"才是直接給定的量, 這裡驗證確實原樣
       反映在回傳的history_F裡)。
案例C: 塑鉸完全塑性(R_post_yield=0, 真正的機構)時, 力控制應該在
       F恰好等於Mp/L(靜力學算出的降伏底剪力)時優雅停止(mechanism_reached
       =True), 不是丟例外讓呼叫端崩潰——這是力控制方法本身的極限
       (沒辦法穿越沒有殘餘勁度的狀態), 不是bug。
案例D: 位移控制模式(control_mode預設'displacement')的行為完全不受
       這次重構影響, 跟之前的行為逐位元一致。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.dofmanager import initial_hinge_states, solve_dofmanager
from frame2d.pushover import run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


def cantilever_no_hinge():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec')
    f.fix(0)
    return f


# ---- 案例A: 彈性範圍內, 力控制結果應該恰好等於標準懸臂梁公式 ----
print("=== 案例A: 力控制彈性範圍結果對得上F*L^3/(3EI) ===")
f = cantilever_no_hinge()
hs = initial_hinge_states(f)   # 空字典(沒有member設定Mp_i), 純彈性
F_target = 500.0
u, F, event_log, hs_final, mech = run_pushover(
    f, hs, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
    target_total=F_target, d_nominal=50.0, base_reaction_dofs=[f.dofs_of(0)[0]],
    control_mode='force',
)
u_theory = F_target * L**3 / (3 * E * I)
assert abs(u[-1] - u_theory) / u_theory < 1e-6, (
    f"力控制解出的位移應該對得上標準公式F*L^3/(3EI)={u_theory}, 實際={u[-1]}"
)
print(f"PASS: 力控制數值解={u[-1]:.10f}, 手算F*L^3/(3EI)={u_theory:.10f}, 相對誤差<1e-6\n")


# ---- 案例B: history_F應該恰好等於施加的力 ----
print("=== 案例B: history_F應該恰好等於直接施加的力(不是解出來的) ===")
assert np.isclose(F[-1], F_target, rtol=1e-9), f"力控制的最終底剪力應該恰好等於target_total={F_target}, 實際={F[-1]}"
print(f"PASS: 最終底剪力={F[-1]:.6f}恰好等於施加的目標力{F_target}\n")


# ---- 案例C: 完全塑性(真正機構)時, 力控制應該優雅停止在Mp/L ----
print("=== 案例C: 完全塑性時力控制優雅停止, 不crash ===")
Mp = 100.0
f2 = Frame2D()
f2.add_node(0, 0, 0); f2.add_node(1, 0, L)
f2.add_section('sec', E=E, I=I, A=A)
f2.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=Mp, Mp_j=1e30,
              R_post_yield_i=0.0, R_post_yield_j=0.0)
f2.fix(0)
hs2 = initial_hinge_states(f2)
u2, F2, ev2, hs2_final, mech2 = run_pushover(
    f2, hs2, prescribed_dofs=[f2.dofs_of(1)[0]], direction=[1.0],
    target_total=1e9, d_nominal=5.0, base_reaction_dofs=[f2.dofs_of(0)[0]],
    control_mode='force', max_steps=1000,
)
assert mech2 is True, "完全塑性(R=0)真正變成機構後, 力控制應該回報mechanism_reached=True"
F_yield_theory = Mp / L
assert abs(F2[-1] - F_yield_theory) < 1e-6, (
    f"力控制在完全塑性下應該恰好停在Mp/L={F_yield_theory}, 實際={F2[-1]}"
)
print(f"PASS: 完全塑性時力控制優雅停在F={F2[-1]:.6f}(=Mp/L={F_yield_theory}), mechanism_reached=True, 沒有crash\n")


# ---- 案例D: 位移控制模式(預設)完全不受這次重構影響 ----
print("=== 案例D: 位移控制模式行為不受force控制重構影響 ===")
f3 = cantilever_no_hinge()
hs3 = initial_hinge_states(f3)
u3, F3, ev3, hs3f, mech3 = run_pushover(
    f3, hs3, prescribed_dofs=[f3.dofs_of(1)[0]], direction=[1.0],
    target_total=0.05, d_nominal=0.01, base_reaction_dofs=[f3.dofs_of(0)[0]],
)   # control_mode沒傳, 預設'displacement'
K_elastic_theory = 3 * E * I / L**3
K_numeric = (F3[1] - F3[0]) / (u3[1] - u3[0])
assert np.isclose(K_numeric, K_elastic_theory, rtol=1e-6), "位移控制模式(預設)結果不應該被這次重構影響"
print(f"PASS: 位移控制(預設)彈性斜率{K_numeric:.4f}跟理論值{K_elastic_theory:.4f}吻合, 沒有被重構影響\n")

print("PASS: 力控制所有案例通過")
