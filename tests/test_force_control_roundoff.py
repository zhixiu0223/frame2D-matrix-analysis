"""
迴歸測試: 力控制的機構偵測不能依賴「恰好為零」的主元

背景(在 Termux / arm64 / Python 3.14 上實際踩到):
test_pushover_force_control.py 案例C(完全塑性、R_post_yield=0)在 x86 上會通過,
在 Termux 上卻丟出
    RuntimeError: 側推超過1000步仍未達到target_total ...

原因: run_pushover(control_mode='force') 判斷機構的唯一方式是
`solve_force_increment()` 裡 np.linalg.solve 丟出 LinAlgError, 而 LAPACK 只有在
遇到**恰好為零**的主元時才會丟。完全塑性的鉸在 x86 上剛好留下恰好0的主元;
換一個平台(FMA、不同版本的 BLAS/LAPACK), 同樣的運算會留下 1e-13 量級的殘餘,
solve 照常成功、解出巨大的位移增量, 迴圈永遠停不下來。

這支測試不依賴平台: 它把勁度矩陣加上 1e-16~1e-13 相對量級的對稱雜訊(正比於彈性勁度的
幾何平均, 模擬不同平台的捨入差異; 1 個 eps 是 2e-16, 這個範圍涵蓋到 500 個 eps), 再跑同一個案例。修正前這會丟出 1000 步的 RuntimeError, 修正後
(solve_force_increment 另外用「彈性參考對角」縮放後的最小特徵值判斷奇異)應該優雅停在 Mp/L;
同時確認合法但跨度大、降伏後仍有 0.1% 勁度的結構不會被誤判成機構。
"""
import numpy as np

import frame2d.pushover as pushover
from frame2d import Frame2D
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0
Mp = 100.0


def plastic_cantilever():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=Mp, Mp_j=1e30,
                 R_post_yield_i=0.0, R_post_yield_j=0.0)
    f.fix(0)
    return f


def run_case():
    f = plastic_cantilever()
    hs = initial_hinge_states(f)
    return run_pushover(
        f, hs, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
        target_total=1e9, d_nominal=5.0, base_reaction_dofs=[f.dofs_of(0)[0]],
        control_mode='force', max_steps=1000,
    )


print("=== 無雜訊(x86上恰好奇異的情況) ===")
u, F, ev, hs, mech = run_case()
assert mech is True and abs(F[-1] - Mp / L) < 1e-6
print(f"PASS: 停在 F={F[-1]:.6f} (=Mp/L={Mp / L}), mechanism_reached=True")

print("=== 加上對稱雜訊(每個元素 ∝ √(彈性對角_i·對角_j) × 1e-16 ~ 1e-13, 模擬捨入差異) ===")
orig = pushover._assemble_stiffness_with_hinges
rng = np.random.default_rng(20260921)
NOISE = [1e-16]


_REF = {}


def noisy(*args, **kwargs):
    K, member_dofs, member_T, member_L = orig(*args, **kwargs)
    if 'scale' not in _REF:                       # 捨入誤差正比於「相消的那些項」的量級 = 彈性勁度對角的幾何平均
        d = np.diag(orig(args[0], None)[0])
        _REF['scale'] = np.sqrt(np.outer(d, d))
    N = rng.standard_normal(K.shape)
    K = K + NOISE[0] * _REF['scale'] * (N + N.T) / 2.0
    return K, member_dofs, member_T, member_L


pushover._assemble_stiffness_with_hinges = noisy
try:
    for level in (1e-16, 1e-15, 1e-14, 1e-13):
        NOISE[0] = level
        for trial in range(3):                   # 每個量級3組不同的雜訊, 都必須優雅停止
            u, F, ev, hs, mech = run_case()
            assert mech is True, f"雜訊{level:g}第{trial + 1}組: 完全塑性後力控制應該回報mechanism_reached=True"
            assert abs(F[-1] - Mp / L) < 1e-4, f"雜訊{level:g}第{trial + 1}組: 應該停在Mp/L={Mp / L}, 實際{F[-1]}"
        print(f"PASS: 雜訊量級 {level:g}: 3組都停在 F={F[-1]:.6f}, mechanism_reached=True, 沒有跑滿1000步")
finally:
    pushover._assemble_stiffness_with_hinges = orig

print("=== 不誤判: 合法但勁度跨度大的結構(A預設1e8), 後降伏勁度很小但不是零, 力控制不能被當成機構 ===")
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, 0, L)
f.add_section('sec', E=E, I=I)                     # A預設1e8: 軸向勁度比彎曲大十幾個數量級
f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=Mp, Mp_j=1e30,
             R_post_yield_i=1e-3, R_post_yield_j=0.0)   # 降伏後還有 0.1% 的勁度
f.fix(0)
hs = initial_hinge_states(f)
u, F, ev, hs, mech = run_pushover(
    f, hs, prescribed_dofs=[f.dofs_of(1)[0]], direction=[1.0],
    target_total=40.0, d_nominal=5.0, base_reaction_dofs=[f.dofs_of(0)[0]],
    control_mode='force', max_steps=1000)
assert mech is False and abs(F[-1] - 40.0) < 1e-6, f"有後降伏勁度的結構不應被誤判成機構: mech={mech}, F={F[-1]}"
assert len(ev) == 1, "應該剛好一次降伏事件"
print(f"PASS: 後降伏勁度0.1%(A=1e8): 順利推到目標力 F={F[-1]:.1f}, 沒有被誤判成機構, 降伏事件 {len(ev)} 次")

print("\nPASS: 力控制的機構偵測不再依賴恰好為零的主元")
