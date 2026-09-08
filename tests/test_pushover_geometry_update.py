"""
驗證案例: run_pushover(..., geometry_update=True) -- 幾何更新模式

案例A: 預設(geometry_update=False)完全不影響既有行為(已經被
       tests/目錄下所有原本的pushover測試間接驗證過, 因為它們都沒有
       傳這個新參數, 全部維持通過就是最直接的證據)。
案例B: 小位移時, geometry_update開關應該給出幾乎相同的結果(幾何本來
       就沒怎麼變, 兩種模式理論上應該收斂到同一個答案)。
案例C: 大位移時, 兩種模式應該有明顯差異(這才是這個功能存在的意義:
       幾何更新模式在大變形時應該更接近真實行為, 跟凍結幾何的線性化
       P-Delta不一樣)。
案例D: 幾何更新模式下, 位移控制應該還是精確推到target_total(不會因為
       幾何更新就讓位移控制本身失準), 底剪力應該是有限的正常數值
       (不會NaN或發散)。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover

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


def run(target, geometry_update):
    f = cantilever_no_hinge()
    hs = initial_hinge_states(f)
    return run_pushover(
        f, hs, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
        target_total=target, d_nominal=target / 100, base_reaction_dofs=[f.dofs_of(0)[0]],
        geometry_update=geometry_update, max_steps=2000,
    )


# ---- 案例B: 小位移時兩種模式應該幾乎相同 ----
print("=== 案例B: 小位移(柱長的0.25%)時, 幾何更新開關應該幾乎不影響結果 ===")
u_off, F_off, _, _, _ = run(0.01, False)
u_on, F_on, _, _, _ = run(0.01, True)
rel_diff_small = abs(F_on[-1] - F_off[-1]) / F_off[-1]
assert rel_diff_small < 1e-4, f"小位移時兩種模式應該幾乎相同, 實際相對差異={rel_diff_small:.2e}"
print(f"PASS: 小位移時相對差異={rel_diff_small:.2e}(<1e-4)\n")


# ---- 案例C: 大位移時應該有明顯差異 ----
print("=== 案例C: 大位移(柱長的37.5%)時, 兩種模式應該有明顯差異 ===")
u_off2, F_off2, _, _, _ = run(1.5, False)
u_on2, F_on2, _, _, _ = run(1.5, True)
rel_diff_large = abs(F_on2[-1] - F_off2[-1]) / F_off2[-1]
assert rel_diff_large > 0.01, f"大位移時兩種模式應該有明顯差異(>1%), 實際={rel_diff_large:.2e}"
print(f"PASS: 大位移時相對差異={rel_diff_large:.2%}(關F={F_off2[-1]:.2f}, 開F={F_on2[-1]:.2f})\n")


# ---- 案例D: 幾何更新模式下位移控制依然精確、底剪力是正常有限值 ----
print("=== 案例D: 幾何更新模式下位移控制精確、結果不是NaN/發散 ===")
assert np.isclose(u_on2[-1], 1.5, rtol=1e-9), f"幾何更新模式下位移控制應該還是精確推到target, 實際={u_on2[-1]}"
assert np.isfinite(F_on2[-1]) and F_on2[-1] > 0, f"底剪力應該是正常的正有限值, 實際={F_on2[-1]}"
print(f"PASS: 位移精確推到{u_on2[-1]}, 底剪力={F_on2[-1]:.2f}(有限正值, 沒有發散)\n")

print("PASS: 幾何更新模式所有案例通過")
