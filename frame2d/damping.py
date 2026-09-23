"""
Rayleigh 阻尼 (動力分析 D5)。

Rayleigh 阻尼: C = α·M + β·K, 對模態 n 的阻尼比是
    ζ_n = α/(2ω_n) + β·ω_n/2
只在兩個控制頻率 ω_i、ω_j 上精確等於指定值, 中間頻率偏低、兩端外的頻率偏高(這是 Rayleigh
阻尼本身的物理限制, 不是這裡實作的問題)。

**這裡的 C 建立在凝縮後的系統上(C_dd = α·M_dd + β·K_eff), 不是原始的 K**
--------------------------------------------------------------------------------
`newmark.py` 模組開頭已經預告過這個設計: 如果直接用未凝縮的原始 K 組
`C_full = α·M + β·K`, 無質量DOF那幾列/行會因為 β·K 這一項而不是0(質量矩陣M在那裡才是0,
K不是), newmark_integrate() 的靜力凝縮假設(無質量DOF是純代數約束、不含阻尼)就會被打破,
變成微分-代數方程式(DAE), 遠比現在複雜。

解法: 用跟 `newmark.condense_for_dynamics()` 完全同一套凝縮公式, 先把 K 縮到 K_eff(只在
有質量的DOF上), 阻尼矩陣建立在 K_eff 上(`C_dd = α·M_dd + β·K_eff`), 再嵌回一個
n×n 的全域矩陣(只有 dyn×dyn 那個區塊非零, 無質量DOF那幾列/行**精確是0**, 因為根本没有
放任何東西進去)。這樣組出來的 `C_full` 送進 `newmark_integrate(damping_matrix=...)`
一定會通過它的「無質量DOF上必須是0」檢查, 不會有 DAE 問題。

單位: 跟輸入一致(質量=力·s²/長度, 跟D1~D4同一套); ω 單位 rad/s(跟 modal.eigen() 一致)。
"""
import numpy as np

from .newmark import condense_for_dynamics


def rayleigh_coefficients(wi, wj, zeta_i, zeta_j=None):
    """由兩個控制頻率(rad/s)與目標阻尼比反解 Rayleigh 係數 (α, β)。

    zeta_j: None(預設) = 兩個頻率用同一個阻尼比 zeta_i(最常見的用法, 化簡成經典公式
        α=2ζωᵢωⱼ/(ωᵢ+ωⱼ), β=2ζ/(ωᵢ+ωⱼ)); 給定的話, 解一般的 2×2 線性方程式
        [1/(2ωᵢ), ωᵢ/2; 1/(2ωⱼ), ωⱼ/2]·[α;β] = [ζᵢ;ζⱼ]。
    """
    if not (wi > 0 and np.isfinite(wi)):
        raise ValueError(f"wi必須是正數, 收到{wi}")
    if not (wj > 0 and np.isfinite(wj)):
        raise ValueError(f"wj必須是正數, 收到{wj}")
    if wi == wj:
        raise ValueError("wi跟wj不能相等(兩個控制頻率必須不同, 否則反解α、β的方程式是奇異的)")
    if zeta_j is None:
        zeta_j = zeta_i
    if not (zeta_i >= 0 and np.isfinite(zeta_i)) or not (zeta_j >= 0 and np.isfinite(zeta_j)):
        raise ValueError(f"zeta_i、zeta_j必須是不小於0的有限數字, 收到{zeta_i}, {zeta_j}")
    A = np.array([[1.0 / (2 * wi), wi / 2.0], [1.0 / (2 * wj), wj / 2.0]])
    alpha, beta = np.linalg.solve(A, np.array([zeta_i, zeta_j]))
    return float(alpha), float(beta)


def rayleigh_damping_ratio(alpha, beta, omega):
    """給定 (α, β), 某個頻率 ω(rad/s, 純量或陣列) 實際會得到的阻尼比 ζ(ω)=α/(2ω)+β·ω/2。
    用來檢核「控制頻率之外的其他模態」實際的阻尼比是多少(Rayleigh 阻尼在兩個控制頻率之間
    偏低、兩端之外偏高, 這是它本身的物理限制)。"""
    omega = np.asarray(omega, dtype=float)
    if np.any(omega <= 0):
        raise ValueError("omega必須都是正數")
    return alpha / (2.0 * omega) + beta * omega / 2.0


def rayleigh_damping_matrix(frame, wi, wj, zeta_i, zeta_j=None, mass_kind='lumped'):
    """組 Rayleigh 阻尼矩陣, 可以直接傳給 `newmark.newmark_integrate(damping_matrix=...)`。

    wi, wj, zeta_i, zeta_j: 見 `rayleigh_coefficients()`。
    mass_kind: 跟後面呼叫 `newmark_integrate()` 用的要一樣('lumped' 或 'consistent'), 否則
        凝縮公式用的 M/K 跟積分時實際用的不一致。

    回傳 (C_full, alpha, beta): C_full 是 n_dof×n_dof 的全域阻尼矩陣(只有有質量的自由DOF
    那個區塊非零), alpha/beta 是反解出來的 Rayleigh 係數(可以自己拿去call
    `rayleigh_damping_ratio()` 檢查其他模態的阻尼比)。
    """
    alpha, beta = rayleigh_coefficients(wi, wj, zeta_i, zeta_j)
    cond = condense_for_dynamics(frame, mass_kind)
    Cdd = alpha * cond.Mdd + beta * cond.Keff
    C_full = np.zeros((cond.n, cond.n))
    C_full[np.ix_(cond.dyn, cond.dyn)] = Cdd
    return C_full, alpha, beta
