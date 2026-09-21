"""
特徵值 / 模態分析 (動力分析 D2)。

解廣義特徵值問題  (K - ω²M) φ = 0, 只用 numpy(不需要 scipy):

1. 取 assemble_K() 與 assemble_M() (DOF編號相同)
2. 分類DOF:
     固定DOF      -- 支承拘束的DOF, φ = 0
     不活動DOF    -- K與M的對角都是0(例如純桁架節點的轉角), 直接排除
     動力DOF      -- 有質量的自由DOF
     無質量DOF    -- 自由但沒有質量(集中質量下的轉角、release端專屬DOF...)
3. 對無質量DOF做靜力凝縮 (對無質量DOF是**精確**的, 不是Guyan近似):
     無質量DOF的動力方程是  K_md φ_d + K_mm φ_m = 0  (慣性項為0)
     => φ_m = -K_mm⁻¹ K_md φ_d,   K_eff = K_dd - K_dm K_mm⁻¹ K_md
4. 對角縮放: S = diag(1/√K_dd,ii) (凝縮**前**的對角), K_s = S K_eff S, M_s = S M_d S。
   結構矩陣常常是「梯度式」的(例如 A 預設 1e8 讓軸向勁度比彎曲勁度大 1e13 倍),
   不縮放的話 eigh 的絕對誤差 ~ eps·‖K‖ 會吃掉低階模態; 縮放是同餘變換,
   特徵值 ω² 完全不變, 但低階模態的精度大幅改善。
5. 化成標準特徵值問題: M_s = L Lᵀ (Cholesky),
     (L⁻¹ K_s L⁻ᵀ) y = ω² y,   φ_d = S L⁻ᵀ y
   => φ_dᵀ M_d φ_d = yᵀy = 1: 模態形狀自動對質量矩陣正規化
   機構偵測看的是「縮放後的勁度矩陣 K_s 的最小特徵值」(單位對角, 健康結構是O(1)以上,
   機構/剛體模態趨近0), 不是 ω² 的相對大小——後者會把「很剛的軸向 + 很軟的彎曲」
   這種合法但跨度大的結構誤判成機構。
6. 回填 φ_m, 固定DOF補0, 得到全DOF的模態形狀

模態性質 (對方向 d ∈ {x, y}, r = 影響向量):
  參與係數          Γ_i  = φ_iᵀ M r          (φ已對M正規化, 所以分母=1; r是完整的
                                              地面運動向量, 含支承DOF, 所以與支承DOF
                                              耦合的一致質量項也算進去)
  有效模態質量      M*_i = Γ_i²
  總質量            M_total = rᵀ M r          (所有質量, 含直接掛在支承上的)
  可動質量          M_free  = p_dᵀ M_dd⁻¹ p_d,  p_d = (M r) 在動力DOF上的分量
                    = 取全部模態時 Σ M*_i 的精確值(動力完備性)。
                    集中質量下等於 M_total 扣掉支承DOF上的節點質量; 一致質量下還要
                    扣掉與支承DOF耦合的部分(那部分質量直接跟著地面動, 模態激發不到)。
  累積質量比        cum_ratio       = Σ_{j<=i} M*_j / M_free   (取全部模態必為1, 完備性檢核)
                    cum_ratio_total = Σ_{j<=i} M*_j / M_total  (規範檢核「≥90%總質量」用這個, 較保守)

正負號慣例: 每個模態形狀正規化成「絕對值最大的動力DOF分量為正」, 讓結果可重現。
單位: 跟輸入一致(質量 = 力·s²/長度), ω 單位 rad/s, T 單位 s, f 單位 Hz。

明確拒絕(ValueError): 結構有機構或支承不足(零頻率/負特徵值)、沒有任何質量、
質量矩陣在動力DOF上不正定、以及 mass.check_dynamic_supported() 的限制
(cable / equal_dof / 非零指定位移)。
"""
from dataclasses import dataclass, field

import numpy as np

from .assembly import assemble_K
from .mass import assemble_M, influence_vector
from .model import Frame2D


@dataclass
class Modal:
    """eigen()的結果。所有陣列的DOF維度跟 assemble_K() 相同(節點DOF在前、
    release端專屬DOF在後)。"""
    omega: np.ndarray             # (n_modes,) 圓頻率 rad/s
    period: np.ndarray            # (n_modes,) 週期 s
    frequency: np.ndarray         # (n_modes,) 頻率 Hz
    phi: np.ndarray               # (n_dof, n_modes) 對M正規化的模態形狀 (φᵀMφ = I)
    gamma: dict                   # {'x': (n_modes,), 'y': (n_modes,)} 參與係數
    eff_mass: dict                # {'x': ..., 'y': ...} 有效模態質量 Γ²
    cum_ratio: dict               # {'x': ..., 'y': ...} 累積有效質量比(相對可動質量)
    cum_ratio_total: dict         # {'x': ..., 'y': ...} 累積有效質量比(相對總質量)
    mass_free: dict               # {'x': 可動質量, 'y': ...}
    mass_total: dict              # {'x': 全部質量(含直接掛在支承上的), 'y': ...}
    kind: str                     # 'lumped' | 'consistent'
    n_dof: int
    n_node_dof: int
    frame: Frame2D = field(repr=False, default=None)

    @property
    def n_modes(self) -> int:
        return len(self.omega)

    def node_shape(self, mode: int, node_id: int):
        """第 mode 個模態(0起算)在某節點的 (ux, uy, rot)。"""
        ux, uy, rot = self.frame.dofs_of(node_id)
        return (float(self.phi[ux, mode]), float(self.phi[uy, mode]), float(self.phi[rot, mode]))

    def modes_needed(self, ratio: float = 0.9, direction: str = 'x'):
        """累積有效質量比達到 ratio 所需的最少模態數; 取到的模態數不夠回傳None。"""
        cum = self.cum_ratio[direction]
        idx = np.where(cum >= ratio - 1e-12)[0]
        return int(idx[0]) + 1 if len(idx) else None

    def table(self) -> str:
        """人讀的模態摘要表。"""
        lines = [f"{'mode':>4} {'T (s)':>11} {'f (Hz)':>11} {'ω (rad/s)':>12} "
                 f"{'Γx':>10} {'M*x/M':>8} {'ΣM*x/M':>8} {'M*y/M':>8} {'ΣM*y/M':>8}"]
        for i in range(self.n_modes):
            def frac(d):
                mf = self.mass_free[d]
                return self.eff_mass[d][i] / mf if mf > 0 else float('nan')
            lines.append(
                f"{i + 1:>4} {self.period[i]:>11.5f} {self.frequency[i]:>11.4f} {self.omega[i]:>12.4f} "
                f"{self.gamma['x'][i]:>10.4f} {frac('x'):>8.4f} {self.cum_ratio['x'][i]:>8.4f} "
                f"{frac('y'):>8.4f} {self.cum_ratio['y'][i]:>8.4f}")
        return "\n".join(lines)


def _fixed_mask(frame: Frame2D, n: int) -> np.ndarray:
    fixed = np.zeros(n, dtype=bool)
    for s in frame.supports:
        for dof, val in zip(frame.dofs_of(s.node), (s.ux, s.uy, s.rot)):
            if val is not None:       # assemble_M已檢查過: 只可能是0.0
                fixed[dof] = True
    return fixed


def eigen(frame: Frame2D, n_modes: int = None, mass: str = 'lumped',
          zero_tol: float = 1e-10) -> Modal:
    """模態分析。

    n_modes: 要幾個模態(從最低頻開始); None=全部(動力DOF數)。
    mass: 'lumped'(預設) 或 'consistent', 見 frame2d.mass。
    zero_tol: 縮放後(單位對角)的勁度矩陣最小特徵值低於這個值就視為機構/支承不足,
        直接報錯, 不會把零頻率模態混進結果。
    """
    M = assemble_M(frame, mass)            # 同時做動力分析的模型檢查
    asm = assemble_K(frame)
    K = asm.K
    n = asm.n_dof
    fixed = _fixed_mask(frame, n)

    # 「完全沒有勁度/質量」一律用整列是否恰為0判斷, 不用相對門檻: 結構勁度常常跨越
    # 十幾個數量級(A預設1e8), 相對門檻會把合法的彎曲DOF誤判成不活動。
    k_zero = np.abs(K).max(axis=1) == 0.0
    m_zero = np.abs(M).max(axis=1) == 0.0
    inactive = (~fixed) & k_zero & m_zero
    free = (~fixed) & (~inactive)
    dyn_mask = free & (~m_zero)
    mless_mask = free & m_zero
    dyn = np.where(dyn_mask)[0]
    mless = np.where(mless_mask)[0]
    if len(dyn) == 0:
        raise ValueError("沒有任何有質量的自由DOF: 請用 Section.rho 或 add_mass() 給結構質量, "
                         "且質量不能全部落在被支承拘束的DOF上。")

    Kdd = K[np.ix_(dyn, dyn)]
    if len(mless):
        Kdm = K[np.ix_(dyn, mless)]
        Kmm = K[np.ix_(mless, mless)]
        try:
            G = np.linalg.solve(Kmm, Kdm.T)          # K_mm⁻¹ K_md
        except np.linalg.LinAlgError:
            raise ValueError("無質量DOF的勁度矩陣奇異, 無法做靜力凝縮: 結構可能有機構"
                             "(例如某個節點只靠release端連接、沒有任何勁度)。請檢查模型。")
        Keff = Kdd - Kdm @ G
    else:
        G = np.zeros((0, len(dyn)))
        Keff = Kdd
    Keff = 0.5 * (Keff + Keff.T)

    Md = M[np.ix_(dyn, dyn)]
    dk = np.diag(Keff)
    if np.any(dk <= 0):
        bad = dyn[np.where(dk <= 0)[0]]
        raise ValueError(f"有質量的DOF {bad.tolist()[:5]} 的有效勁度為0或負: 支承不足(剛體運動)或結構有機構"
                         "(例如質量掛在只靠鉸接連著的節點上)。請檢查支承與桿件連接。")
    # 縮放參考用「凝縮前」的對角 K_dd,ii, 不是凝縮後的 K_eff,ii: 若用凝縮後自己的對角, 勁度已經
    # 掉到雜訊量級(1e-14)的DOF會被放大成單位對角, 看起來反而「健康」, 機構偵測就漏掉了
    # (實測: 鉸支承懸臂在 1e-14 相對雜訊下回傳 ω=4e-6 的假模態)。用凝縮前的對角, 凝縮把
    # 勁度吃光的DOF縮放後就是 ~1e-14, 一眼看得出來; 合法的跨度大結構(A預設1e8)仍然是良態。
    kdd = np.diag(Kdd)
    if np.any(kdd <= 0):
        bad = dyn[np.where(kdd <= 0)[0]]
        raise ValueError(f"有質量的DOF {bad.tolist()[:5]} 完全沒有勁度: 結構有機構。請檢查支承與桿件連接。")
    S = 1.0 / np.sqrt(kdd)
    Ks = (S[:, None] * Keff) * S[None, :]
    Ks = 0.5 * (Ks + Ks.T)
    mu = np.linalg.eigvalsh(Ks)
    if mu[0] <= zero_tol:
        n_bad = int(np.sum(mu <= zero_tol))
        raise ValueError(
            f"結構有 {n_bad} 個機構/剛體模態(縮放後勁度矩陣最小特徵值 = {mu[0]:.3e}): "
            "支承不足、有機構、或某個有質量的DOF沒有被勁度約束住。零頻率模態沒有意義, "
            "請先補足支承/桿件。")
    Ms = (S[:, None] * Md) * S[None, :]
    Ms = 0.5 * (Ms + Ms.T)
    try:
        Lc = np.linalg.cholesky(Ms)
    except np.linalg.LinAlgError:
        raise ValueError("質量矩陣在動力DOF上不是正定的(可能有DOF質量為0卻被判成有質量, "
                         "或一致質量矩陣退化)。請檢查 rho / add_mass 的設定。")
    Linv = np.linalg.inv(Lc)
    A = Linv @ Ks @ Linv.T
    A = 0.5 * (A + A.T)
    lam, Y = np.linalg.eigh(A)                       # 由小到大
    if lam[0] <= 0:
        raise ValueError(f"最小特徵值 ω² = {lam[0]:.3e} <= 0: 勁度矩陣不是正定的(機構或數值問題)。")

    n_avail = len(lam)
    nm = n_avail if n_modes is None else min(int(n_modes), n_avail)
    if nm < 1:
        raise ValueError("n_modes 至少要 1")
    lam, Y = lam[:nm], Y[:, :nm]
    Phi_d = S[:, None] * (Linv.T @ Y)                # 對Md正規化: Phi_dᵀ Md Phi_d = I

    # 正負號慣例: 動力DOF裡絕對值最大的分量為正
    for j in range(nm):
        if Phi_d[np.argmax(np.abs(Phi_d[:, j])), j] < 0:
            Phi_d[:, j] = -Phi_d[:, j]

    phi = np.zeros((n, nm))
    phi[dyn, :] = Phi_d
    if len(mless):
        phi[mless, :] = -G @ Phi_d                   # 回填無質量DOF

    omega = np.sqrt(lam)
    period = 2 * np.pi / omega
    frequency = omega / (2 * np.pi)

    gamma, eff, cum, cum_tot, m_free, m_tot = {}, {}, {}, {}, {}, {}
    for d in ('x', 'y'):
        r = influence_vector(frame, d)
        Mr = M @ r
        g = phi.T @ Mr                               # φ已對M正規化
        gamma[d] = g
        eff[d] = g * g
        p = Mr[dyn]                                  # 有效地震力 -M r a_g 在動力DOF上的分量
        m_free[d] = float(p @ np.linalg.solve(Md, p))    # 動力完備性: 取全部模態時 ΣM* 的精確值
        m_tot[d] = float(r @ Mr)
        cum[d] = (np.cumsum(eff[d]) / m_free[d]) if m_free[d] > 0 else np.full(nm, np.nan)
        cum_tot[d] = (np.cumsum(eff[d]) / m_tot[d]) if m_tot[d] > 0 else np.full(nm, np.nan)

    return Modal(omega=omega, period=period, frequency=frequency, phi=phi,
                 gamma=gamma, eff_mass=eff, cum_ratio=cum, cum_ratio_total=cum_tot,
                 mass_free=m_free, mass_total=m_tot,
                 kind=mass, n_dof=n, n_node_dof=asm.n_node_dof, frame=frame)
