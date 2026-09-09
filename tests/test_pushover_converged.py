"""
驗證案例: run_pushover_converged() -- 幾何平衡疊代版本

跟run_pushover()是完全獨立的兩個函式, 這裡的驗證重點是:
1. 純彈性小位移時對得上解析解(跟run_pushover()一樣的基本要求)。
2. 小位移時跟run_pushover()(開geometry_update)給出幾乎一樣的結果
   (疊代應該幾乎立刻收斂, 兩者在這個範圍內在數學上該收斂到同一個答案)。
3. 大位移/大轉角時, 不會出現run_pushover()(一次線性外推)那種"突然
   變硬"的假象——這裡驗證的是曲線本身應該維持合理的單調趨勢, 不是
   驗證數值本身多精確(那需要跟更完整的大變形理論解比較, 這裡沒有做)。
4. 步長設太大導致疊代不收斂時, 會誠實回報mechanism_reached=True提前
   停止, 不會給一個不可信的答案。
5. 誠實記錄一個重要限制: 疊代收斂只保證"幾何/軸力自洽", 不保證轉角
   還在小角度假設的有效範圍內——這是兩個獨立的問題, 疊代收斂不會
   自動讓大轉角失真的疑慮消失, max_rotation這個既有的警告機制在這裡
   還是必要的。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover, run_pushover_converged

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


def cantilever_no_hinge():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec',
                 Mp_i=1e30, Mp_j=1e30, R_post_yield_i=1.0, R_post_yield_j=1.0)
    f.fix(0)
    return f


# ---- 案例1: 純彈性小位移對得上解析解 ----
print("=== 案例1: 力控制彈性範圍結果對得上F*L^3/(3EI) ===")
f = cantilever_no_hinge()
hs = initial_hinge_states(f)
F_target = 500.0
u, F, event_log, hs_final, mech = run_pushover_converged(
    f, hs, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
    target_total=F_target, d_nominal=50.0, base_reaction_dofs=[f.dofs_of(0)[0]],
    control_mode='force', use_pdelta=False, geometry_update=False,
)
u_theory = F_target * L**3 / (3 * E * I)
assert abs(u[-1] - u_theory) / u_theory < 1e-6, f"應該對得上解析解, 實際相對誤差={abs(u[-1]-u_theory)/u_theory}"
assert mech is False, "純彈性小位移不應該回報mechanism_reached"
print(f"PASS: 數值解={u[-1]:.10f}, 解析解={u_theory:.10f}, 相對誤差<1e-6\n")


# ---- 案例2: 小位移時跟run_pushover(geometry_update=True)幾乎一致 ----
print("=== 案例2: 小位移時跟run_pushover()(開geometry_update)結果幾乎一致 ===")
f2 = cantilever_no_hinge()
hs2 = initial_hinge_states(f2)
u_old, F_old, _, _, _ = run_pushover(
    f2, hs2, prescribed_dofs=[f2.dofs_of(1)[0]], direction=[1.0],
    target_total=0.01, d_nominal=0.001, base_reaction_dofs=[f2.dofs_of(0)[0]],
    use_pdelta=True, geometry_update=True,
)
f3 = cantilever_no_hinge()
hs3 = initial_hinge_states(f3)
u_new, F_new, _, _, mech3 = run_pushover_converged(
    f3, hs3, prescribed_dofs=[f3.dofs_of(1)[0]], direction=[1.0],
    target_total=0.01, d_nominal=0.001, base_reaction_dofs=[f3.dofs_of(0)[0]],
    use_pdelta=True, geometry_update=True,
)
rel_diff = abs(F_new[-1] - F_old[-1]) / F_old[-1]
assert rel_diff < 1e-4, f"小位移時新舊版本應該幾乎一致, 實際相對差異={rel_diff:.2e}"
assert mech3 is False
print(f"PASS: 小位移(柱長的0.25%)時相對差異={rel_diff:.2e}(<1e-4)\n")


# ---- 案例3: 大位移時曲線維持單調(不會有"突然變硬"的假象) ----
print("=== 案例3: 大位移時F-u曲線維持單調遞增, 沒有run_pushover()那種假象 ===")
f4 = cantilever_no_hinge()
hs4 = initial_hinge_states(f4)
u4, F4, ev4, hsf4, mech4 = run_pushover_converged(
    f4, hs4, prescribed_dofs=[f4.dofs_of(1)[0]], direction=[1.0],
    target_total=1.5, d_nominal=0.02, base_reaction_dofs=[f4.dofs_of(0)[0]],
    control_mode='displacement', use_pdelta=True, geometry_update=True, max_steps=1000,
)
du = np.diff(u4)
dF = np.diff(F4)
assert np.all(du >= -1e-12), "位移控制下位移歷程應該單調不減"
# 切線剛度dF/du不應該出現數量級跳升(這是之前發現的假象特徵)
tangent = dF[du > 1e-9] / du[du > 1e-9]
assert tangent.max() / max(tangent.min(), 1e-9) < 10, (
    f"切線剛度不應該有數量級跳升(之前的假象是好幾倍甚至幾十倍暴衝), "
    f"實際最大/最小比值={tangent.max()/max(tangent.min(),1e-9):.2f}"
)
print(f"PASS: 全程{len(u4)}步位移單調, 切線剛度沒有數量級跳升(最大/最小比值={tangent.max()/max(tangent.min(),1e-9):.2f})\n")


# ---- 案例4: 步長設太大時, 誠實回報不收斂, 不會硬給答案 ----
print("=== 案例4: 步長太大時應該誠實回報mechanism_reached(=不收斂) ===")
f5 = cantilever_no_hinge()
hs5 = initial_hinge_states(f5)
u5, F5, ev5, hsf5, mech5 = run_pushover_converged(
    f5, hs5, prescribed_dofs=[f5.dofs_of(1)[0]], direction=[1.0],
    target_total=3.0, d_nominal=3.0,  # 一步就想推完整個3m(柱長的75%), 步長極端大
    base_reaction_dofs=[f5.dofs_of(0)[0]],
    control_mode='displacement', use_pdelta=True, geometry_update=True,
    max_geom_iter=5,  # 疊代次數也故意壓低, 增加不收斂機率
)
assert mech5 is True, "步長極端大+疊代次數壓低時, 應該真的偵測到不收斂"
print(f"PASS: 步長太大時正確回報mechanism_reached=True, 在u={u5[-1]:.4f}處提前停止, 沒有給出不可信的答案\n")


# ---- 案例5: 誠實記錄的限制 -- 疊代收斂不代表轉角還在小角度範圍內 ----
print("=== 案例5(誠實記錄限制): 疊代收斂不保證轉角仍在小角度有效範圍 ===")
f6 = cantilever_no_hinge()
hs6 = initial_hinge_states(f6)
u6, F6, ev6, hsf6, mech6, maxrot6 = run_pushover_converged(
    f6, hs6, prescribed_dofs=[f6.dofs_of(1)[0]], direction=[1.0],
    target_total=1.5, d_nominal=0.02, base_reaction_dofs=[f6.dofs_of(0)[0]],
    control_mode='displacement', use_pdelta=True, geometry_update=True, max_steps=1000,
    include_max_rotation=True,
)
assert mech6 is False, "這個案例的疊代本身應該是收斂的(數值上自洽)"
assert np.degrees(maxrot6) > 5, (
    "刻意驗證: 就算疊代收斂(數值自洽), 轉角依然可能遠超過小角度假設的"
    "有效範圍——這是兩個獨立的問題, max_rotation這個既有警告機制在"
    "這個新solver裡仍然必要, 不能因為疊代收斂就以為結果一定可信"
)
print(f"PASS(誠實記錄): 疊代收斂(mechanism_reached=False), 但最大轉角={np.degrees(maxrot6):.1f}度, "
      f"仍然超出小角度假設有效範圍——收斂性≠物理有效性, 這是刻意驗證的限制, 不是bug\n")

print("PASS: run_pushover_converged()所有案例通過")
