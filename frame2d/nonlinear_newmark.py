"""
非線性時程分析: Newmark-β + 循環塑鉸 (動力分析 D6)。

把 D7(`cyclic.py`)已經驗證過的雙線性運動硬化循環塑鉸(`CyclicHingeState`), 接到 D4/D5
(`newmark.py`/`damping.py`)的 Newmark 直接積分, 得到非線性(有塑鉸)的地震反應歷程。

**為什麼直接接塑鉸, 不是先做「彈性 corotational」**(這是規劃時的更正, 記錄在這裡): 原始
ROADMAP 把 D6 規劃成「先用彈性 corotational 把時間迴圈本身驗證乾淨, 再接塑鉸」。但塑鉸模型
(雙線性、event-to-event)已經在 D7 用解析解跟 OpenSeesPy 驗證到機器精度, 而且 event-to-event
對分段線性系統是**精確**的, 不需要 Newton-Raphson 疊代收斂——直接把它接上 Newmark 比另外做一個
「先驗證迴圈本身」的過渡階段更省工、風險更低。程式碼刻意最大化重用 D7 的元件
(`_assemble_stiffness_with_hinges`、`_member_force_increments`、`_yield_events`、
`_spring_rotation_increment`), 這裡新增的部分只有「怎麼把它們跟 Newmark 的時間積分接起來」。

演算法: 每個時間步都是一次「力控制」的 event-to-event 增量分析
--------------------------------------------------------------------------------
D7 的 `run_cyclic()` 是**位移控制**(已知位移增量, 解對應的力), 這裡反過來是**力控制**
(已知 Newmark 有效力增量 ΔF_hat, 解對應的位移增量) —— 這其實更簡單, 因為不需要 D7 那套
「受控 vs 自由」DOF 的區分, 全部自由DOF都是要解的:

  每個時間步 t_k → t_{k+1}(Δt固定不變):
    ΔF_hat = (F(t_{k+1})-F(t_k)) + M·(c1·v_k+c2·a_k) + C·(c3·v_k+c4·a_k)     (Newmark標準式,
                                                                                跟線性版同一套)
    remaining ← ΔF_hat
    while remaining 還沒完全消化:
        用目前的塑鉸狀態組 K(彈性/硬化, 視降伏狀態而定), K_hat = K + 質量/阻尼係數項
        (跟 D7 一樣先做卸載符號一致性疊代: 若試解導致已降伏鉸的轉角增量方向反了, 該鉸切回
        彈性、重新試解)
        用力控制解 K_hat·Δu_trial = remaining(全自由度, 不是D7的部分受控)
        找出這個 Δu_trial 內第一個會被觸發的降伏事件, 得到比例 ratio∈(0,1]
        套用 ratio·Δu_trial(累加進這個時間步的 Δu_total), remaining ← (1-ratio)·remaining
        若 ratio<1: 更新剛降伏的鉸狀態, 重新組 K、重複; 若 ratio=1: 這個時間步結束
    再用標準 Newmark 運動學關係式(跟線性版完全相同的公式)由 Δu_total 算出 Δv、Δa

這個 event-to-event 子步驟*不是*近似 —— 只要塑鉸是分段線性(雙線性運動硬化), 這個演算法在
每個線性段內是精確解, 沒有 Newton-Raphson 收斂容差的問題。

**目前的範圍(誠實)**:
- 準靜態、小位移(不含 P-Delta / corotational 幾何非線性), 跟 D7 一致
- **不對無質量DOF做靜力凝縮**(跟D4/D5不同!): 為了避免塑鉸狀態改變時必須重新凝縮的複雜度,
  這裡直接在完整系統上解, 但這要求**初始條件必須物理一致**——具體來說, 動態擾動一定要從
  u=v=a=0(相對於某個已經用 `pushover.apply_gravity()` 算到平衡的靜態起點, 透過
  `initial_cum_forces` 參數代入, 概念上跟 D7 的 `initial_cum_forces` 完全一樣)開始, 不接受
  非零的初始位移/速度(D4 已經踩過這個坑, 見 newmark.py 模組說明; 這裡乾脆直接不開放這個
  選項, 而不是重蹈覆轍)
- 沒有勁度/強度劣化、沒有捏縮, 跟 D7 一致
- 阻尼矩陣(Rayleigh)是常數, 不隨塑鉸降伏而改變(業界常見做法之一, 也有人主張降伏後應該
  用切線剛度比例阻尼——這裡沒有實作那個選項)
"""
from dataclasses import dataclass, field

import numpy as np

from .assembly import assemble_K
from .cyclic import CyclicHingeState, _spring_rotation_increment, _yield_events
from .mass import assemble_M
from .pushover import (M_LOCAL_IDX, _assemble_stiffness_with_hinges, _fixed_dof_set,
                       _member_force_increments, solve_force_increment)


@dataclass
class NonlinearNewmarkResult:
    """nonlinear_newmark_integrate() 的回傳值。t/u/v/a 是完整DOF的時間序列(n_steps+1, n_dof),
    hinge_M/hinge_theta_p 是每個曾經降伏的塑鉸端的彎矩/累積塑性轉角時間序列
    (n_steps+1,)、member_force_history 是每根桿件兩端局部內力的時間序列(n_steps+1, 6)
    (直接從逐步累積的 cum_forces 記錄下來, 對降伏後的桿件仍然正確——跟線性引擎不同,
    不能用 k_local@u 反推, 因為降伏後 k_local 不再是常數)。"""
    t: np.ndarray
    u: np.ndarray
    v: np.ndarray
    a: np.ndarray
    hinge_M: dict
    hinge_theta_p: dict
    hinge_work: dict
    member_force_history: dict
    events: list
    hinge_states: dict
    beta: float
    gamma: float
    frame: object = field(repr=False, default=None)
    n_hinge_solves: int = 0        # 診斷用: 總共做了幾次 K_hat 組裝+解(給效能檢查用)

    @property
    def n_dof(self) -> int:
        return self.u.shape[1]

    @property
    def n_steps(self) -> int:
        return self.u.shape[0] - 1

    def dof_history(self, node_id, direction='x'):
        li = {'x': 0, 'y': 1, 'rot': 2}[direction]
        dof = self.frame.dofs_of(node_id)[li]
        return self.u[:, dof], self.v[:, dof], self.a[:, dof]


def nonlinear_newmark_integrate(frame, hinge_states, dt, n_steps, force=None, mass_kind='lumped',
                                damping_matrix=None, initial_cum_forces=None, beta=0.25, gamma=0.5,
                                max_unload_iter=50, max_substeps_per_step=500) -> NonlinearNewmarkResult:
    """非線性時程分析(Newmark-β + 循環塑鉸, event-to-event)。

    hinge_states: {member_id: CyclicHingeState}(用 `CyclicHingeState.from_hinge_states()` 轉換)。
    dt, n_steps, force, mass_kind, damping_matrix, beta, gamma: 跟 `newmark.newmark_integrate()`
        意義相同(force 的無質量DOF限制在這裡**沒有**, 因為這裡不做靜力凝縮——但無質量DOF
        還是不能是完全孤立、沒有任何勁度的機構)。
    initial_cum_forces: 選用, 重力預載後的桿端內力起點(`pushover.apply_gravity()` 算出來的
        那個), 跟 D7 的 `run_cyclic()` 同一個用法。動態擾動一律從 u=v=a=0 開始(相對於這個
        預載狀態), 不支援非零初始位移/速度(見模組說明)。

    明確拒絕: dt<=0、n_steps<1、hinge_states 裡有非 CyclicHingeState 的物件、模型沒有支承/
    沒有質量、K_hat 奇異(機構)、卸載疊代不收斂、單一時間步內子步數超過 max_substeps_per_step
    (通常代表 d_nominal 等級的問題不會發生在這裡, 但保留一個安全閥避免無窮迴圈)。
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
    for mid, hs in hinge_states.items():
        if not isinstance(hs, CyclicHingeState):
            raise TypeError(f"member {mid} 的塑鉸狀態必須是 CyclicHingeState "
                            "(可用 CyclicHingeState.from_hinge_states(...) 轉換)")

    n = assemble_K(frame).n_dof
    M = assemble_M(frame, mass_kind)           # 同時做動力分析的模型檢查
    C = np.zeros((n, n)) if damping_matrix is None else np.asarray(damping_matrix, dtype=float)
    if C.shape != (n, n):
        raise ValueError(f"damping_matrix的形狀應該是({n},{n}), 收到{C.shape}")
    fixed_dofs = _fixed_dof_set(frame)
    if len(fixed_dofs) == n:
        raise ValueError("結構完全被固定, 沒有任何自由度可以動態分析")
    free_dofs = [d for d in range(n) if d not in fixed_dofs]
    if np.abs(M[np.ix_(free_dofs, free_dofs)]).max() == 0.0:
        raise ValueError("沒有任何有質量的自由DOF: 請用 Section.rho 或 add_mass() 給結構質量, "
                         "且質量不能全部落在被支承拘束的DOF上。")

    cum_forces = ({mid: np.zeros(6) for mid in frame.members} if initial_cum_forces is None
                  else {mid: np.array(v, dtype=float) for mid, v in initial_cum_forces.items()})
    for mid, hs in hinge_states.items():        # 起點必須還在彈性區間內(跟 run_cyclic 相同檢查)
        for e in (0, 1):
            if np.isfinite(hs.Mp[e]):
                M0 = cum_forces[mid][M_LOCAL_IDX[e]]
                lo, hi = hs.elastic_range(e)
                if M0 > hi * (1 + 1e-9) + 1e-12 or M0 < lo * (1 + 1e-9) - 1e-12:
                    raise ValueError(
                        f"起始狀態(例如重力預載)下 member {mid} 的{'i' if e == 0 else 'j'}端彎矩 "
                        f"{M0:.6g} 已超出塑性彎矩 Mp = {hs.Mp[e]:.6g}, 請降低重力載重或提高 Mp。")

    def F(t):
        if force is None:
            return np.zeros(n)
        f = np.asarray(force(t), dtype=float)
        if f.shape != (n,):
            raise ValueError(f"force(t)必須回傳長度{n}的向量, 在 t={t} 收到形狀{f.shape}")
        return f

    t = np.arange(n_steps + 1) * dt
    u = np.zeros((n_steps + 1, n))
    v = np.zeros((n_steps + 1, n))
    a = np.zeros((n_steps + 1, n))

    keys = [(mid, e) for mid in hinge_states for e in (0, 1)]
    hist_M = {k: [float(cum_forces[k[0]][M_LOCAL_IDX[k[1]]])] for k in keys}
    hist_tp = {k: [0.0] for k in keys}
    hist_work = {k: [0.0] for k in keys}
    mf_hist = {mid: [cum_forces[mid].copy()] for mid in frame.members}
    events = []
    n_solves = 0

    c1, c2 = 1.0 / (beta * dt), 1.0 / (2 * beta)
    c3, c4 = gamma / beta, dt * (gamma / (2 * beta) - 1.0)

    Fk = F(0.0)
    for k in range(n_steps):
        Fk1 = F(t[k + 1])
        dF_hat_total = ((Fk1 - Fk) + M @ (c1 * v[k] + c2 * a[k]) + C @ (c3 * v[k] + c4 * a[k]))
        remaining = dF_hat_total.copy()
        du_total = np.zeros(n)
        n_sub = 0
        while True:
            n_sub += 1
            if n_sub > max_substeps_per_step:
                raise RuntimeError(f"第 {k + 1} 個時間步的子步數超過 {max_substeps_per_step}, 已中止"
                                   "(可能是模型在反覆降伏/卸載, 或 max_substeps_per_step 設太小)。")
            unloaded_now = []
            for _ in range(max_unload_iter):
                K, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(frame, hinge_states)
                Khat = K + (gamma / (beta * dt)) * C + (1.0 / (beta * dt**2)) * M
                try:
                    du_trial = solve_force_increment(Khat, list(range(n)), remaining, fixed_dofs)
                except RuntimeError as exc:
                    zero_R = [f"M{mid} {'i' if e == 0 else 'j'}端" for mid, hs in hinge_states.items()
                              for e in (0, 1) if hs.yielded[e] and hs.R_post_yield[e] <= 0.0]
                    hint = ""
                    if zero_R:
                        hint = (f" 常見原因: 已降伏的塑鉸 {', '.join(zero_R)} 的硬化剛度 R_post_yield = 0"
                                "(完全塑性)。當兩個 R=0 的鉸夾住同一個沒有轉動慣量的節點時, 那個節點的"
                                "轉角沒有任何勁度或質量, 矩陣就奇異了。請把 R_post_yield 改成小的正值"
                                "(例如 0.01·EI/L), 或給該節點一個小的轉動慣量。")
                    raise RuntimeError(f"第 {k + 1} 個時間步(t={t[k]:.6g})有效勁度矩陣奇異或不正定: "
                                       "結構在目前的塑鉸狀態下可能形成機構。" + hint
                                       + " 原始訊息: " + str(exc))
                n_solves += 1
                changed = False
                for mid, hs in hinge_states.items():
                    if not any(hs.yielded):
                        continue
                    dth = _spring_rotation_increment(frame, hs, mid, member_T, member_dofs, member_L, du_trial)
                    for e in (0, 1):
                        if hs.yielded[e] and hs.direction[e] * dth[e] < -1e-12:
                            M0 = cum_forces[mid][M_LOCAL_IDX[e]]
                            hs.alpha[e] = M0 - hs.direction[e] * hs.Mp[e]
                            hs.yielded[e] = False
                            unloaded_now.append(('unload', mid, e))
                            changed = True
                if not changed:
                    break
            else:
                raise RuntimeError(f"第 {k + 1} 個時間步卸載的符號一致性疊代 {max_unload_iter} 次仍不收斂。")

            df = _member_force_increments(frame, member_T, member_dofs, member_L, hinge_states, None, du_trial)
            evs = _yield_events(hinge_states, cum_forces, df)
            ratio = min((r for r, *_ in evs), default=1.0)
            du_use = du_trial * ratio

            for mid, hs in hinge_states.items():
                if not any(hs.yielded):
                    continue
                dth = _spring_rotation_increment(frame, hs, mid, member_T, member_dofs, member_L, du_use)
                for e in (0, 1):
                    if hs.yielded[e]:
                        M0 = cum_forces[mid][M_LOCAL_IDX[e]]
                        dM = ratio * df[mid][M_LOCAL_IDX[e]]
                        hs.theta_p_signed[e] += dth[e]
                        hs.theta_p[e] += abs(dth[e])
                        hs.work[e] += (M0 + 0.5 * dM) * dth[e]
            for mid in frame.members:
                cum_forces[mid] += ratio * df[mid]
            du_total += du_use
            remaining = remaining * (1.0 - ratio)

            newly = list(unloaded_now)
            if evs:
                for r, mid, e, s in evs:
                    if abs(r - ratio) < 1e-9:
                        hinge_states[mid].yielded[e] = True
                        hinge_states[mid].direction[e] = s
                        newly.append(('yield', mid, e))
            for kind, mid, e in newly:
                events.append({'kind': kind, 'member': mid, 'end': e, 'step': k + 1, 't': float(t[k + 1])})
            if ratio >= 1.0 - 1e-12:
                break

        da = (du_total - dt * v[k] - (dt**2 / 2.0) * a[k]) / (beta * dt**2)
        dv = gamma * dt * da + dt * a[k]
        u[k + 1] = u[k] + du_total
        v[k + 1] = v[k] + dv
        a[k + 1] = a[k] + da
        for key in keys:
            mid, e = key
            hist_M[key].append(float(cum_forces[mid][M_LOCAL_IDX[e]]))
            hist_tp[key].append(hinge_states[mid].theta_p_signed[e])
            hist_work[key].append(hinge_states[mid].work[e])
        for mid in frame.members:
            mf_hist[mid].append(cum_forces[mid].copy())
        Fk = Fk1

    return NonlinearNewmarkResult(
        t=t, u=u, v=v, a=a,
        hinge_M={k: np.array(v_) for k, v_ in hist_M.items()},
        hinge_theta_p={k: np.array(v_) for k, v_ in hist_tp.items()},
        hinge_work={k: np.array(v_) for k, v_ in hist_work.items()},
        member_force_history={mid: np.array(v_) for mid, v_ in mf_hist.items()},
        events=events, hinge_states=hinge_states, beta=beta, gamma=gamma, frame=frame,
        n_hinge_solves=n_solves)
