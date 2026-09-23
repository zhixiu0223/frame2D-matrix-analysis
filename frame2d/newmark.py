"""
線性時程分析 (Newmark-β 直接積分) (動力分析 D4)。

解線性運動方程式 M ü + C u̇ + K u = F(t), 只用 numpy, 用 Newmark-β 法逐步積分。預設
β=1/4、γ=1/2(平均加速度法), 無條件穩定(任何 Δt 都不會數值發散, 但精度仍然跟 Δt/T 有關,
Δt 太大週期會被拉長, 見 tests/test_newmark.py 的收斂測試)。

無質量自由度: 一定要先靜力凝縮(不是可有可無的細節, 是這個模組唯一容易做錯的地方)
--------------------------------------------------------------------------------
集中質量下轉角 DOF、release 端專屬 DOF 通常沒有質量。這種 DOF 的運動方程式那一列質量恆為0,
代表它們**不是**動態自由度, 而是每個時刻都要滿足的靜力平衡約束:
    0 = F_m(t) − K_mm u_m(t) − K_md u_d(t)          (m=無質量DOF, d=有質量DOF; 這裡不支援
                                                       C_mm/C_md非零, 見下面說明)
    => u_m(t) = K_mm⁻¹ (F_m(t) − K_md u_d(t))
代入有質量 DOF 的方程式, 消去 u_m, 得到**只在有質量 DOF 上**的縮小系統:
    M_dd ü_d + C_dd u̇_d + K_eff u_d = F_d(t) − K_dm K_mm⁻¹ F_m(t)
    K_eff = K_dd − K_dm K_mm⁻¹ K_md                  (跟 modal.eigen() 的靜力凝縮是同一個公式)
Newmark 只在這個縮小系統上跑, 跑完再用同一個線性關係回填 u_m(t)、v_m(t)、a_m(t)
(這個關係式對時間微分兩次仍然成立, 因為 K_mm⁻¹K_md 是常數矩陣, 不隨時間變)。

**這是修正過的設計**: 最早的版本沒有做這個凝縮, 直接把完整系統(含無質量DOF)一起丟進
Newmark 遞迴。數學上 K_hat solve 本身沒問題(K_hat 在無質量DOF那一列還是由 K 撐住, 不會奇異),
但**初始條件**(使用者給的 u0/v0)如果在無質量DOF上不是靜力平衡一致的, 從第一步就會錯: 在
一個頂端質量懸臂柱的自由振動測試中量到最大位移誤差是初始振幅的 6 倍(不是誤差百分比, 是誤差
本身比訊號還大), 一路查到才發現是這個問題(見這裡的說明與 tests/test_newmark.py 的還原測試)。

**目前不支援**: 直接施加在無質量自由度上的力(例如對一個沒有轉動慣量的節點施加彎矩時程),
以及阻尼矩陣在無質量DOF上有非零項(C_mm、C_md ≠ 0)。後者是因為 D5 的 Rayleigh 阻尼
C=αM+βK 的 βK 項在無質量DOF上不為0, 那時無質量DOF就不再是純代數約束, 而是變成
微分-代數方程式(DAE), 遠比這裡複雜; 目前的設計刻意先把這個複雜度擋在外面, 等 D5/D6 需要
時再處理(而且 D5 的 Rayleigh 阻尼應該建立在**凝縮後**的 K_eff 上, 不是建立在原始 K 上,
這樣就不會遇到這個問題——這一點先記錄在這裡)。

單位: 跟輸入一致(質量=力·s²/長度, 跟D1~D3同一套)。
"""
from dataclasses import dataclass, field

import numpy as np

from .assembly import assemble_K
from .elements import member_stiffness_local, member_stiffness_local_truss
from .mass import assemble_M
from .pushover import _fixed_dof_set  # noqa: F401  (保留給其他模組沿用同一個 import 路徑)


@dataclass
class NewmarkResult:
    """newmark_integrate() 的回傳值。所有時間序列陣列的第0軸是時間步(長度 n_steps+1, 含 t=0
    的初始狀態), 第1軸是 DOF(長度 n_dof, 跟 assemble_K 同一套編號, 含無質量DOF, 已經回填)。"""
    t: np.ndarray                 # (n_steps+1,)
    u: np.ndarray                 # (n_steps+1, n_dof) 位移
    v: np.ndarray                 # (n_steps+1, n_dof) 速度(無質量DOF由同一組線性關係回填)
    a: np.ndarray                 # (n_steps+1, n_dof) 加速度(同上)
    beta: float
    gamma: float
    mass_kind: str
    dyn_dofs: np.ndarray = field(repr=False, default=None)     # 真正被積分的(有質量)DOF
    frame: object = field(repr=False, default=None)
    K: np.ndarray = field(repr=False, default=None)
    M: np.ndarray = field(repr=False, default=None)
    C: np.ndarray = field(repr=False, default=None)
    member_dofs: dict = field(repr=False, default=None)
    member_T: dict = field(repr=False, default=None)
    member_L: dict = field(repr=False, default=None)

    @property
    def n_dof(self) -> int:
        return self.u.shape[1]

    @property
    def n_steps(self) -> int:
        return self.u.shape[0] - 1

    def dof_history(self, node_id, direction='x'):
        """某節點某方向(x/y/rot)的位移/速度/加速度時間歷程, 回傳 (u, v, a) 三條長度
        n_steps+1 的陣列。"""
        li = {'x': 0, 'y': 1, 'rot': 2}[direction]
        dof = self.frame.dofs_of(node_id)[li]
        return self.u[:, dof], self.v[:, dof], self.a[:, dof]

    def member_force_history(self, member_id):
        """某桿件兩端的局部內力時間歷程(純線性彈性, 不含塑鉸/幾何非線性), 回傳
        (n_steps+1, 6) 陣列, 欄位順序 (Fx_i, Fy_i, M_i, Fx_j, Fy_j, M_j)。"""
        m = self.frame.members[member_id]
        section = self.frame.sections[m.section]
        L = self.member_L[member_id]
        T = self.member_T[member_id]
        idx = np.array(self.member_dofs[member_id])
        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            k_local = member_stiffness_local(section.E, section.I, section.A, L)
        return (k_local @ (T @ self.u[:, idx].T)).T


def _fixed_mask(frame, n):
    fixed = np.zeros(n, dtype=bool)
    for s in frame.supports:
        for dof, val in zip(frame.dofs_of(s.node), (s.ux, s.uy, s.rot)):
            if val is not None:
                fixed[dof] = True
    return fixed


@dataclass
class Condensation:
    """`condense_for_dynamics()` 的回傳值: 一次算好的 DOF 分類與靜力凝縮結果, 給
    `newmark_integrate()` 跟 `damping.rayleigh_damping_matrix()` 共用(單一事實來源, 避免
    兩邊各自重算、以後改公式時漏改一邊)。"""
    K: np.ndarray
    M: np.ndarray
    Keff: np.ndarray          # (n_dyn, n_dyn) 凝縮後的有效勁度
    Mdd: np.ndarray           # (n_dyn, n_dyn)
    dyn: np.ndarray           # 有質量的自由DOF, 對應 assemble_K 的DOF編號
    mless: np.ndarray         # 無質量的自由DOF
    fixed: np.ndarray         # (n,) bool, 支承(固定)DOF
    Kdm: np.ndarray           # (n_dyn, n_mless), 沒有無質量DOF時是 (n_dyn, 0)
    Kmm_inv: np.ndarray       # (n_mless, n_mless), 沒有無質量DOF時是 None
    G: np.ndarray             # Kmm_inv @ Kmd, u_m = -G u_d + Kmm_inv F_m; 沒有無質量DOF時是 (0, n_dyn)
    n: int                    # 全部DOF數(跟 assemble_K 一致)
    member_dofs: dict = field(repr=False, default=None)
    member_T: dict = field(repr=False, default=None)
    member_L: dict = field(repr=False, default=None)

    @property
    def n_dyn(self) -> int:
        return len(self.dyn)


def condense_for_dynamics(frame, mass_kind='lumped') -> Condensation:
    """組 K、M, 分類 DOF(固定/有質量的自由/無質量的自由), 並對無質量DOF做靜力凝縮
    (K_eff = K_dd − K_dm K_mm⁻¹ K_md), 見模組開頭的完整說明。`newmark_integrate()` 跟
    `damping.rayleigh_damping_matrix()` 都呼叫這個函式, 確保兩邊的凝縮公式永遠一致。"""
    asm = assemble_K(frame)
    K = asm.K
    M = assemble_M(frame, mass_kind)          # 同時做動力分析的模型檢查(cable/equal_dof/支承沉陷)
    n = asm.n_dof

    fixed = _fixed_mask(frame, n)
    k_zero = np.abs(K).max(axis=1) == 0.0
    m_zero = np.abs(M).max(axis=1) == 0.0
    inactive = (~fixed) & k_zero & m_zero
    free = (~fixed) & (~inactive)
    dyn = np.where(free & (~m_zero))[0]
    mless = np.where(free & m_zero)[0]
    if len(dyn) == 0:
        raise ValueError("沒有任何有質量的自由DOF: 請用 Section.rho 或 add_mass() 給結構質量, "
                         "且質量不能全部落在被支承拘束的DOF上。")
    if len(mless):
        Kdd, Kdm, Kmd, Kmm = K[np.ix_(dyn, dyn)], K[np.ix_(dyn, mless)], K[np.ix_(mless, dyn)], K[np.ix_(mless, mless)]
        try:
            Kmm_inv = np.linalg.inv(Kmm)
        except np.linalg.LinAlgError:
            raise ValueError("無質量DOF的勁度矩陣奇異, 無法做靜力凝縮: 結構可能有機構"
                             "(例如某個節點只靠release端連接、沒有任何勁度)。請檢查模型。")
        G = Kmm_inv @ Kmd
        Keff = Kdd - Kdm @ G
    else:
        Kdm = np.zeros((len(dyn), 0))
        Kmm_inv = None
        G = np.zeros((0, len(dyn)))
        Keff = K[np.ix_(dyn, dyn)]

    return Condensation(K=K, M=M, Keff=Keff, Mdd=M[np.ix_(dyn, dyn)], dyn=dyn, mless=mless,
                        fixed=fixed, Kdm=Kdm, Kmm_inv=Kmm_inv, G=G, n=n,
                        member_dofs=asm.member_dofs, member_T=asm.member_T, member_L=asm.member_L)


def newmark_integrate(frame, dt, n_steps, force=None, mass_kind='lumped', damping_matrix=None,
                      initial_disp=None, initial_vel=None, beta=0.25, gamma=0.5) -> NewmarkResult:
    """線性時程分析(Newmark-β 直接積分)。無質量自由度會自動靜力凝縮(見模組開頭說明), 使用端
    不用自己處理。

    dt: 時間步長(s), 必須是正數。n_steps: 要走幾步(結果會有 n_steps+1 個時間點, 含 t=0)。
    force: callable, force(t) -> 長度n_dof的力向量(跟 assemble_K 同一套DOF編號), 通常用
        `excitation.force_series_from_pattern()` 組出來; None = 沒有外力(自由振動, 通常會搭配
        非零的 initial_disp/initial_vel 才有意義)。**不支援**在無質量自由度上有非零的力。
    mass_kind: 'lumped'(預設) 或 'consistent', 見 frame2d.mass。
    damping_matrix: 選用, n_dof×n_dof 的阻尼矩陣 C; None = 零阻尼。**必須**在無質量自由度上
        整列整行都是0(見模組開頭關於 DAE 複雜度的說明)。D5 會提供 Rayleigh 阻尼的建構函式。
    initial_disp / initial_vel: 選用, 長度n_dof的初始位移/速度(只有有質量的自由度上的值會被
        採用; 無質量自由度上的值會被忽略, 一律用靜力凝縮關係式重新算, 因為無質量自由度沒有
        獨立的初始條件可言)。None = 有質量自由度全部從0開始。
    beta, gamma: Newmark 参数, 預設是平均加速度法(無條件穩定、無數值阻尼)。線性加速度法是
        β=1/6, γ=1/2(條件穩定, Δt 必須夠小); γ<0.5 時 Newmark 法有負的數值阻尼, 不穩定, 不接受。

    明確拒絕: dt<=0、n_steps<1、β<=0、γ<0.5、force(t)回傳長度不對或在無質量DOF上非零、
    damping_matrix在無質量DOF上非零、模型沒有支承/沒有質量/有機構(沿用 assemble_M 與靜力
    凝縮的訊息)。
    """
    if not (dt > 0 and np.isfinite(dt)):
        raise ValueError(f"dt必須是正數, 收到{dt}")
    if int(n_steps) < 1 or int(n_steps) != n_steps:
        raise ValueError(f"n_steps必須是不小於1的整數, 收到{n_steps}")
    n_steps = int(n_steps)
    if beta <= 0:
        raise ValueError(f"beta必須是正數, 收到{beta}")
    if gamma < 0.5:
        raise ValueError(f"gamma必須 >= 0.5(Newmark法在gamma<0.5時會有負的數值阻尼, 不穩定), 收到{gamma}")

    cond = condense_for_dynamics(frame, mass_kind)
    K, M, Keff, Mdd = cond.K, cond.M, cond.Keff, cond.Mdd
    dyn, mless, fixed = cond.dyn, cond.mless, cond.fixed
    Kdm, Kmm_inv, G = cond.Kdm, cond.Kmm_inv, cond.G
    n = cond.n
    C_full = np.zeros((n, n)) if damping_matrix is None else np.asarray(damping_matrix, dtype=float)
    if C_full.shape != (n, n):
        raise ValueError(f"damping_matrix的形狀應該是({n},{n}), 收到{C_full.shape}")
    if len(mless):
        if np.abs(C_full[np.ix_(mless, dyn)]).max() > 0 or np.abs(C_full[np.ix_(mless, mless)]).max() > 0 \
                or np.abs(C_full[np.ix_(dyn, mless)]).max() > 0:
            raise ValueError("damping_matrix在無質量自由度上有非零項: 目前的靜力凝縮假設無質量"
                             "DOF是純代數約束(不含阻尼), 見 newmark.py 模組說明的 DAE 限制。")

    Cdd = C_full[np.ix_(dyn, dyn)]

    u0_full = np.zeros(n) if initial_disp is None else np.array(initial_disp, dtype=float)
    v0_full = np.zeros(n) if initial_vel is None else np.array(initial_vel, dtype=float)
    if u0_full.shape != (n,) or v0_full.shape != (n,):
        raise ValueError(f"initial_disp/initial_vel的長度應該是{n}")
    if np.any(u0_full[fixed] != 0) or np.any(v0_full[fixed] != 0):
        raise ValueError("initial_disp/initial_vel在支承(固定)自由度上必須是0")
    u0_d, v0_d = u0_full[dyn], v0_full[dyn]

    def F_dyn(t):
        """回傳 (F_d(t), F_m(t)) 拆分, 並驗證無質量DOF上的力是0。"""
        if force is None:
            f = np.zeros(n)
        else:
            f = np.asarray(force(t), dtype=float)
            if f.shape != (n,):
                raise ValueError(f"force(t)必須回傳長度{n}的向量, 在 t={t} 收到形狀{f.shape}")
        if len(mless) and np.abs(f[mless]).max() > 1e-9 * max(1.0, np.abs(f).max()):
            raise ValueError(f"force(t)在無質量自由度 {mless[np.argmax(np.abs(f[mless]))]} 上非零"
                             f"(t={t}): 目前不支援直接施加在無質量自由度上的力(例如對沒有轉動"
                             "慣量的節點施加彎矩時程)。請改成對有質量的自由度施加力, 或給該節點"
                             "一個小的轉動慣量(add_mass的Iz參數)。")
        return f[dyn], (f[mless] if len(mless) else np.zeros(0))

    def F_eff(t):
        f_d, f_m = F_dyn(t)
        if len(mless):
            return f_d - Kdm @ (Kmm_inv @ f_m)
        return f_d

    Khat = Keff + (gamma / (beta * dt)) * Cdd + (1.0 / (beta * dt**2)) * Mdd
    try:
        Khat_inv = np.linalg.inv(Khat)
    except np.linalg.LinAlgError:
        raise ValueError("有效勁度矩陣 K_hat 是奇異的, 無法積分: 結構在有質量的自由DOF上可能有"
                         "機構。")

    nd = len(dyn)
    t = np.arange(n_steps + 1) * dt
    ud = np.zeros((n_steps + 1, nd))
    vd = np.zeros((n_steps + 1, nd))
    ad = np.zeros((n_steps + 1, nd))
    ud[0], vd[0] = u0_d, v0_d
    F0 = F_eff(0.0)
    ad[0] = np.linalg.solve(Mdd, F0 - Cdd @ v0_d - Keff @ u0_d)

    c1, c2 = 1.0 / (beta * dt), 1.0 / (2 * beta)
    c3, c4 = gamma / beta, dt * (gamma / (2 * beta) - 1.0)
    Fk = F0
    for k in range(n_steps):
        Fk1 = F_eff(t[k + 1])
        dF = Fk1 - Fk
        dF_hat = dF + Mdd @ (c1 * vd[k] + c2 * ad[k]) + Cdd @ (c3 * vd[k] + c4 * ad[k])
        du = Khat_inv @ dF_hat
        da = (du - dt * vd[k] - (dt**2 / 2.0) * ad[k]) / (beta * dt**2)
        dv = gamma * dt * da + dt * ad[k]
        ud[k + 1] = ud[k] + du
        vd[k + 1] = vd[k] + dv
        ad[k + 1] = ad[k] + da
        Fk = Fk1

    u = np.zeros((n_steps + 1, n))
    v = np.zeros((n_steps + 1, n))
    a = np.zeros((n_steps + 1, n))
    u[:, dyn], v[:, dyn], a[:, dyn] = ud, vd, ad
    if len(mless):
        Fm_series = np.array([F_dyn(tt)[1] for tt in t])           # (n_steps+1, len(mless))
        particular = Fm_series @ Kmm_inv.T                          # Kmm_inv @ F_m, 逐步
        u[:, mless] = particular - ud @ G.T
        v[:, mless] = -vd @ G.T                                     # F_m 恆為0時 v_m 的時間導數項才成立
        a[:, mless] = -ad @ G.T

    return NewmarkResult(t=t, u=u, v=v, a=a, beta=beta, gamma=gamma, mass_kind=mass_kind, dyn_dofs=dyn,
                         frame=frame, K=K, M=M, C=C_full, member_dofs=cond.member_dofs,
                         member_T=cond.member_T, member_L=cond.member_L)
