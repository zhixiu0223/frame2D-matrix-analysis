"""
反應譜分析 (動力分析 D3)。

反應譜分析完全建立在模態分析(D2, `modal.eigen()`)上, 是一種**線性**方法, 給的是各模態
最大反應的估計, 不是時間歷程。流程:

  K φ_i = ω_i² M φ_i             (D2 已經做好, 這裡直接拿 Modal 物件用)
  T_i = 2π/ω_i  ->  查表/公式得到 S_a(T_i)
  第 i 個模態的等效地震靜力:  f_i = Γ_i · S_a(T_i) · M r        (r = 影響向量, D2 的 Γ、M*
                                                                  已經是這個定義下的量)
  第 i 個模態的位移:          D_i = Γ_i · S_d(T_i) · φ_i,  S_d(T_i) = S_a(T_i)/ω_i²
                              (D_i 精確滿足 K D_i = f_i, 因為 K φ_i = ω_i² M φ_i,
                              代入即得 —— 不需要另外解一次線性方程式)
  第 i 個模態的桿件內力:      k_local @ (T_member @ D_i[member_dofs])(跟靜力後處理同一套
                              局部勁度矩陣, 純線性彈性, 不含塑鉸/幾何非線性)
  第 i 個模態的反力:          R_i = K @ D_i − f_i(在自由DOF上這一項精確為0, 是內部
                              一致性檢查; 在固定DOF上就是反力)

**每一個反應量(位移分量、桿件內力分量、反力、基底剪力)分別算出各模態的訊號值(有正負號),
再對這組模態值做 SRSS 或 CQC 組合** —— 不是先把位移組合好, 再拿組合後的位移去反推內力或
反力(那樣算出來的桿件內力不會滿足平衡, 也不是文獻定義的反應譜法)。這是本模組最重要的
實作原則, 底下每個函式的回傳值都遵守它。

總基底剪力有一個不需要組裝 K 就能算出來的簡潔公式(標準結果, 見 Chopra《Dynamics of
Structures》): 第 i 模態的基底總剪力 = 有效模態質量 × 譜加速度:
    V_i = Γ_i · L_i · S_a(T_i) = Γ_i² · S_a(T_i) = M*_i · S_a(T_i)
(L_i = Γ_i, 因為 D2 的模態形狀已經對 M 正規化, φᵀMφ=1)。這是全域動量平衡的結果, 對「地面
運動方向上、結構全部支承一起構成剛性地盤」這個假設成立, 不需要知道各支承怎麼分擔。
`response_spectrum()` 同時回傳這個總基底剪力(analytic), 也回傳逐支承的反力(算出來的),
兩者在測試裡互相校驗。

模態組合
--------
SRSS(平方和開根號): combined = sqrt(Σ v_i²)。頻率分離良好(相鄰模態週期比 < 0.9)時的
  標準做法, 只需要各模態的訊號值, 不需要阻尼比。

CQC(完全二次項組合, Der Kiureghian & Rosenblueth 1981):
  combined = sqrt(Σᵢ Σⱼ ρᵢⱼ vᵢ vⱼ)
  ρᵢⱼ = 8ζ²(1+r)r^1.5 / ((1-r²)² + 4ζ²r(1+r)²),  r = ωⱼ/ωᵢ
  頻率相近(密集模態)時比 SRSS 準確; 頻率分離良好時 ρᵢⱼ→0 (i≠j), CQC 收斂到 SRSS。
  這裡假設所有模態同一個阻尼比 ζ(工程上最常見的簡化; 每個模態各自的阻尼比要等 D5
  Rayleigh 阻尼才有)。ρ 是有效的相關係數矩陣(半正定), 所以 Σρvv 恆 ≥ 0, 不需要另外取
  絕對值。

單位: 跟輸入一致(質量=力·s²/長度); `spectrum` 只吃 S_a(T) 的 callable, 這裡不對單位
做任何假設或換算。台灣規範反應譜的 S_a(T) 查表/公式做成 adapter, 放在使用端(例如
taiwan-seismic-code-calc), 不進這個核心模組。
"""
from dataclasses import dataclass, field

import numpy as np

from .assembly import assemble_K
from .elements import member_stiffness_local, member_stiffness_local_truss
from .modal import Modal

VALID_COMBINE = ('SRSS', 'CQC')


def cqc_rho(omega: np.ndarray, damping) -> np.ndarray:
    """CQC 相關係數矩陣 ρ (n×n, 對稱, 對角線=1)。

    damping: 純量(全部模態同一個阻尼比)或長度 n 的陣列(每個模態各自的阻尼比,
    用兩者的平均 —— Der Kiureghian & Rosenblueth 原始公式是單一阻尼比推導的, 不同阻尼比
    時取平均是常見的工程近似)。
    """
    n = len(omega)
    zeta = np.broadcast_to(np.asarray(damping, dtype=float), (n,))
    zi = zeta[:, None]
    zj = zeta[None, :]
    z = 0.5 * (zi + zj)
    r = omega[None, :] / omega[:, None]
    num = 8.0 * z * z * (1.0 + r) * np.power(np.abs(r), 1.5)
    den = (1.0 - r * r) ** 2 + 4.0 * z * z * r * (1.0 + r) ** 2
    rho = num / den
    di = np.arange(n)
    rho[di, di] = 1.0
    return rho


def combine_srss(values: np.ndarray) -> float:
    """SRSS 組合一組模態訊號值(1D 陣列, 沿最後一軸組合)。"""
    return float(np.sqrt(np.sum(np.asarray(values) ** 2)))


def combine_cqc(values: np.ndarray, omega: np.ndarray, damping) -> float:
    """CQC 組合一組模態訊號值。"""
    v = np.asarray(values, dtype=float)
    rho = cqc_rho(np.asarray(omega, dtype=float), damping)
    val = float(v @ rho @ v)
    return float(np.sqrt(max(val, 0.0)))                 # 理論上 rho 半正定, val>=0; 保險截斷浮點負值


def _combine_array(modal_values: np.ndarray, omega: np.ndarray, combine: str, damping) -> np.ndarray:
    """modal_values: (n_modes, ...) 每個模態的訊號值(可以是多維, 例如每個DOF一欄)。
    回傳跟 modal_values.shape[1:] 相同形狀的組合結果, 沿 axis=0 (模態) 組合。"""
    if combine == 'SRSS':
        return np.sqrt(np.sum(modal_values ** 2, axis=0))
    rho = cqc_rho(omega, damping)
    flat = modal_values.reshape(modal_values.shape[0], -1)    # (n_modes, k) k=攤平後的分量數
    # combined_k = sqrt( sum_i sum_j rho_ij * flat[i,k] * flat[j,k] )
    val = np.einsum('ik,ij,jk->k', flat, rho, flat)
    val = np.sqrt(np.clip(val, 0.0, None))
    return val.reshape(modal_values.shape[1:])


@dataclass
class RSAResult:
    """response_spectrum() 的回傳值。全部是**組合後**(SRSS/CQC)的峰值反應估計(一律非負,
    因為模態時間歷程的相對正負號已經丟失, 這是反應譜法的本質限制, 不是實作問題)。
    modal_* 開頭的欄位是組合前、各模態自己的訊號值(有正負號), 供檢查與客製化組合用。
    """
    direction: str
    combine: str
    damping: float
    periods: np.ndarray                  # (n_modes,) 用到的模態週期
    sa: np.ndarray                       # (n_modes,) 各模態的 S_a(T_i)
    displacements: np.ndarray            # (n_dof,) 組合後的位移(每個分量各自 SRSS/CQC, 非負)
    member_forces: dict                  # {member_id: (6,) 組合後的局部端點內力(非負)}
    reactions: dict                      # {dof: 組合後反力(非負)} 只包含支承DOF
    base_shear: float                    # 總基底剪力(SRSS/CQC組合, analytic公式, 見模組說明)
    modal_displacements: np.ndarray      # (n_modes, n_dof) 組合前, 有正負號
    modal_member_forces: dict            # {member_id: (n_modes, 6)} 組合前
    modal_reactions: dict                # {dof: (n_modes,)} 組合前
    modal_base_shear: np.ndarray         # (n_modes,) 組合前, = eff_mass * Sa (有正負號: Γ可正可負)
    cum_ratio_total: float               # 用到的模態相對總質量的累積有效質量比(規範≥90%檢核)
    frame: object = field(repr=False, default=None)


def response_spectrum(modal: Modal, spectrum, direction: str = 'x', damping: float = 0.05,
                      combine: str = 'SRSS', n_modes: int = None) -> RSAResult:
    """反應譜分析。

    modal: `modal.eigen()` 的回傳值(模態分析要先做好)。
    spectrum: callable, spectrum(T) -> S_a(T)(純量, 單位跟 modal 用的質量/長度/時間制一致)。
    direction: 'x' 或 'y', 要跟 modal 算的同一個方向(D2 的 Modal 同時存了 x、y 兩個方向的
        參與係數/有效質量, 這裡直接取用, 不用重算)。
    damping: 阻尼比(0~1), SRSS 不會用到, CQC 用來算相關係數矩陣。
    combine: 'SRSS' 或 'CQC'。
    n_modes: 只用前幾個模態(週期最長、頻率最低的那些); None = modal 裡全部的模態。
    """
    if direction not in ('x', 'y'):
        raise ValueError(f"direction必須是'x'或'y', 收到'{direction}'")
    if combine not in VALID_COMBINE:
        raise ValueError(f"combine必須是{VALID_COMBINE}其中之一, 收到'{combine}'")
    if not (0.0 <= damping < 1.0):
        raise ValueError(f"damping必須在[0, 1)之間, 收到{damping}")
    nm = modal.n_modes if n_modes is None else int(n_modes)
    if nm < 1 or nm > modal.n_modes:
        raise ValueError(f"n_modes必須在1~{modal.n_modes}之間, 收到{n_modes}")

    frame = modal.frame
    omega = modal.omega[:nm]
    periods = modal.period[:nm]
    gamma = modal.gamma[direction][:nm]
    eff_mass = modal.eff_mass[direction][:nm]
    phi = modal.phi[:, :nm]

    sa = np.array([float(spectrum(T)) for T in periods])
    if np.any(~np.isfinite(sa)) or np.any(sa < 0):
        raise ValueError("spectrum(T) 回傳了負值或非有限值(nan/inf), 請檢查反應譜函式")
    sd = sa / omega**2                                     # S_d(T) = S_a(T)/ω²

    # 位移: D_i = Γ_i · S_d(T_i) · φ_i  (精確滿足 K D_i = f_i, 見模組說明, 不用另外解)
    modal_disp = (gamma * sd)[:, None] * phi.T             # (nm, n_dof)

    # 桿件內力: 直接用局部勁度矩陣, 純線性彈性(不含塑鉸/幾何非線性)
    asm = assemble_K(frame)
    modal_forces = {}
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        L = asm.member_L[mid]
        T = asm.member_T[mid]
        idx = np.array(asm.member_dofs[mid])
        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            k_local = member_stiffness_local(section.E, section.I, section.A, L)
        modal_forces[mid] = (k_local @ (T @ modal_disp[:, idx].T)).T     # (nm, 6)

    # 反力: R_i = K @ D_i − f_i_full (自由DOF上應精確為0, 固定DOF上是反力)。
    # 這裡的等效靜力是 Chopra 的「模態慣性力分佈」 s_i* = Γ_i · M · φ_i(不是 M·r!
    # M·r 只用在算 Γ_i/M*_i 本身。K·D_i 的右手邊必須是 s_i*, 因為
    # K D_i = Γ_i S_d,i (K φ_i) = Γ_i S_d,i ω_i² M φ_i = Γ_i S_a,i M φ_i, 這是從
    # K φ_i = ω_i² M φ_i 直接代入得到的恆等式, 見模組開頭說明)。
    fixed_dofs = sorted({d for s in frame.supports for d, val in
                         zip(frame.dofs_of(s.node), (s.ux, s.uy, s.rot)) if val is not None})
    from .mass import assemble_M
    M = assemble_M(frame, modal.kind)
    f_full = (gamma * sa)[:, None] * (M @ phi).T                       # (nm, n_dof), s_i* = Γ_i S_a,i M φ_i
    kd = modal_disp @ asm.K.T                                          # (nm, n_dof), K@D_i
    resid = kd - f_full
    free_dofs = [d for d in range(asm.n_dof) if d not in fixed_dofs]
    if free_dofs:
        max_resid = float(np.max(np.abs(resid[:, free_dofs])))
        scale = max(1.0, float(np.max(np.abs(f_full))))
        if max_resid > 1e-6 * scale:
            raise RuntimeError(
                f"內部一致性檢查失敗: 自由DOF上 K@D_i - f_i 應該精確為0, 實際最大殘差 "
                f"{max_resid:.3e}(相對尺度 {scale:.3e})。這通常代表 assemble_K() 跟 modal.eigen() "
                "用的模型不一致(例如中途改了 frame), 請回報這個訊息。")
    modal_reactions = {d: resid[:, d] for d in fixed_dofs}

    modal_base_shear = eff_mass * sa                                    # V_i = Γ_i² · S_a(T_i)

    if combine == 'SRSS':
        disp_c = _combine_array(modal_disp, omega, 'SRSS', damping)
        forces_c = {mid: _combine_array(v, omega, 'SRSS', damping) for mid, v in modal_forces.items()}
        reactions_c = {d: combine_srss(v) for d, v in modal_reactions.items()}
        base_shear_c = combine_srss(modal_base_shear)
    else:
        disp_c = _combine_array(modal_disp, omega, 'CQC', damping)
        forces_c = {mid: _combine_array(v, omega, 'CQC', damping) for mid, v in modal_forces.items()}
        reactions_c = {d: combine_cqc(v, omega, damping) for d, v in modal_reactions.items()}
        base_shear_c = combine_cqc(modal_base_shear, omega, damping)

    return RSAResult(
        direction=direction, combine=combine, damping=float(damping),
        periods=periods, sa=sa, displacements=disp_c, member_forces=forces_c, reactions=reactions_c,
        base_shear=base_shear_c, modal_displacements=modal_disp, modal_member_forces=modal_forces,
        modal_reactions=modal_reactions, modal_base_shear=modal_base_shear,
        cum_ratio_total=float(modal.cum_ratio_total[direction][nm - 1]), frame=frame)
