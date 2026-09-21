"""
循環(遲滯)塑鉸與反覆載重分析 (動力分析路線 D7)。

這個檔案做兩件事
----------------
1. `CyclicHingeState`: 雙線性**運動硬化**塑鉸(彈性→降伏→卸載→反向降伏), 跟既有的
   `hinge.HingeState`(只能單向:彈性→降伏)並存, **不修改** HingeState。
2. `run_cyclic()`: 位移控制的**反覆載重**(正向推到 +a1、反向拉到 -a1、再 +a2 ...), 沿用
   pushover 已驗證過的元件(含鉸勁度矩陣組裝、位移增量求解、桿端內力增量), 只換掉「單向載重
   迴圈」跟「單向降伏判斷」, 得到 力-位移遲滯迴圈 與 各塑鉸的 彎矩-塑性轉角 迴圈。

遲滯迴圈是什麼、跟動力分析的關係
--------------------------------
遲滯迴圈來自**構件的組成律**(這裡是塑鉸的彎矩-轉角關係)在反覆變形下的行為, 不是「動力」
本身產生的。動力(時程)分析只是提供一條隨時間變化的位移歷程, 把它餵進同一個非線性構件,
構件畫出來的就是遲滯迴圈。所以先用「準靜態反覆載重」(這個檔案)把構件模型做對、對得上
解析解與 OpenSeesPy, 之後 D6/D8 的非線性時程才有可靠的構件可以用。

塑鉸模型(跟 hinge.py 同一個物理: 樑本身彈性 + 端部串接零長度旋轉彈簧)
----------------------------------------------------------------------
彈簧的彎矩-轉角關係是雙線性運動硬化:
  彈性: 剛度 K_e = RIGID_FACTOR·EI/L (數值上代表「剛接」)
  降伏: 降伏面中心(背應力) α, 彈性區間 [α - Mp, α + Mp]; 到達區間邊界就降伏, 降伏後
        彈簧剛度 = R_post_yield(硬化剛度); 反向載重時先彈性卸載, 走到**另一側**邊界
        (α ∓ Mp)才反向降伏 —— 也就是包辛格(Bauschinger)型的運動硬化, 不是等向硬化。
  這跟 OpenSees Steel01(雙線性、b = R_post/K_e、不開等向硬化)是同一個組成律。
  R_post_yield = 0 就是彈塑性(完全塑性)。

狀態機(每個鉸端)
----------------
  yielded[e]  = True: 目前在塑性狀態(彈簧剛度=R_post), False: 彈性(近剛接)
  direction[e]= +1/-1: 塑性狀態的降伏方向(彎矩正/負側)
  alpha[e]    = 背應力(降伏面中心), 卸載時凍結成 M - direction·Mp
事件到事件(event-to-event)演算法:
  一段增量內勁度矩陣固定 ⇒ 反應對增量大小是**線性**的, 所以
  - 「降伏事件」(彈性→塑性)可以用線性外插找出這一步內發生的比例(跟 pushover 相同);
  - 「卸載」不可能發生在一段增量的中間(符號在段內不變), 只會在段的**開頭**被偵測: 用目前的
    塑性勁度試解, 若某塑性端的彈簧轉角增量與降伏方向相反, 就把它切回彈性、重新試解;
    因為一個鉸卸載會讓別的鉸的轉角增量跟著變, 所以這個「符號一致性」要疊代到穩定
    (跟接觸問題的 active-set 一樣), 疊代不收斂就明確報錯, 不給假答案。

範圍(誠實)
----------
- 準靜態、小位移(不含 P-Delta / 幾何更新); 這一版先把構件模型做對。
- 沒有勁度/強度劣化、沒有捏縮(pinching)、沒有等向硬化(Steel01 預設也沒有)。
- 單調載重時結果必須與既有 run_pushover(HingeState) 完全一致(見 tests/test_cyclic.py)。
"""
from dataclasses import dataclass, field

import numpy as np

from .dofmanager import build_dof_map
from .hinge import HingeState, beam_internal_rotation
from .pushover import (M_LOCAL_IDX, _assemble_stiffness_with_hinges, _fixed_dof_set,
                       _member_force_increments, solve_displacement_increment)


class CyclicHingeState(HingeState):
    """雙線性運動硬化循環塑鉸的狀態(見模組說明)。介面跟 HingeState 相容(current_R 等), 所以
    可以直接餵給 pushover 的勁度矩陣組裝; 多出來的欄位描述遲滯歷程:

      yielded[e]        目前是否在塑性狀態
      direction[e]      塑性狀態的方向 (+1/-1), 彈性時保留上一次的方向
      alpha[e]          背應力(降伏面中心); 彈性區間 = [alpha - Mp, alpha + Mp]
      theta_p_signed[e] 帶號累積塑性轉角(rad), 正負反映淨塑性變形方向
      theta_p[e]        累積塑性轉角絕對值 ∑|Δθp| (沿用 HingeState 的欄位, 供性能等級分類用)
      work[e]           累積塑性功 ∫M dθp (含硬化儲存的可回復能量; 完整迴圈內等於耗能)
    """

    def __init__(self, Mp1, Mp2, R_post_yield_1, R_post_yield_2, **kw):
        super().__init__(Mp1, Mp2, R_post_yield_1, R_post_yield_2, **kw)
        self.direction = [0, 0]
        self.alpha = [0.0, 0.0]
        self.theta_p_signed = [0.0, 0.0]
        self.work = [0.0, 0.0]

    @classmethod
    def from_hinge_states(cls, hinge_states: dict) -> dict:
        """把 dofmanager.initial_hinge_states() 產生的 {member_id: HingeState} 轉成循環版
        (容量與硬化剛度照抄, 歷程欄位全部歸零)。"""
        return {mid: cls(hs.Mp[0], hs.Mp[1], hs.R_post_yield[0], hs.R_post_yield[1],
                         theta_IO=tuple(hs.theta_IO), theta_LS=tuple(hs.theta_LS),
                         theta_CP=tuple(hs.theta_CP))
                for mid, hs in hinge_states.items()}

    def check_yield(self, M1, M2):      # 單向版的介面在這裡沒有意義, 明確擋掉避免誤用
        raise NotImplementedError("CyclicHingeState 的降伏判斷由 run_cyclic() 的事件到事件邏輯負責, "
                                  "不要呼叫 check_yield()")

    def elastic_range(self, e):
        """彈性區間 (lo, hi) = (alpha - Mp, alpha + Mp)。"""
        return self.alpha[e] - self.Mp[e], self.alpha[e] + self.Mp[e]


@dataclass
class CyclicResult:
    """run_cyclic() 的回傳值。所有 history 陣列等長, 第 0 個點是初始狀態。"""
    u: np.ndarray                       # 控制自由度的實際位移
    F: np.ndarray                       # 底剪力 (= -Σ 底部反力, 推向為正)
    work_ext: np.ndarray                # 外力累積功 ∫F du
    hinge_M: dict                       # {(member_id, end): 該鉸端彎矩歷程}
    hinge_theta_p: dict                 # {(member_id, end): 帶號累積塑性轉角歷程}
    hinge_work: dict                    # {(member_id, end): 累積塑性功歷程}
    events: list                        # [{'kind': 'yield'|'unload', 'member', 'end', 'u', 'F', 'step'}]
    hinge_states: dict                  # 最終狀態(原地更新後的同一組物件)
    u_full: np.ndarray = field(repr=False, default=None)     # 最終全域位移向量
    cum_forces: dict = field(repr=False, default=None)       # 最終桿端內力
    n_steps: int = 0

    @property
    def total_plastic_work(self) -> float:
        return float(sum(w[-1] for w in self.hinge_work.values()))

    def loop_work(self, i0: int, i1: int) -> float:
        """歷程索引 i0..i1 之間外力做的功 ∫F du (走完整個迴圈時 = 迴圈面積 = 耗能)。"""
        return float(self.work_ext[i1] - self.work_ext[i0])


def _spring_rotation_increment(frame, hs, mid, member_T, member_dofs, member_L, du_full):
    """該增量造成的兩端彈簧相對轉角增量 (Δ(θ1-φ1), Δ(θ2-φ2))。beam_internal_rotation 對
    (v1,θ1,v2,θ2) 是線性的, 所以直接把增量代進去即可(用的是「目前狀態」的 R, 跟組 K 時一致)。"""
    m = frame.members[mid]
    section = frame.sections[m.section]
    du_local = member_T[mid] @ du_full[np.array(member_dofs[mid])]
    v1, th1, v2, th2 = du_local[1], du_local[2], du_local[4], du_local[5]
    phi1, phi2 = beam_internal_rotation(section.E, section.I, member_L[mid], hs, v1, th1, v2, th2)
    return th1 - phi1, th2 - phi2


def _yield_events(hinge_states, cum_forces, df_by_member):
    """彈性狀態的鉸端: 這一步是否會把彎矩推出彈性區間 [alpha-Mp, alpha+Mp]。回傳
    [(ratio, member_id, end, direction)], ratio 是事件發生在這一步的比例 (0~1)。"""
    events = []
    for mid, hs in hinge_states.items():
        for e in (0, 1):
            if hs.yielded[e] or not np.isfinite(hs.Mp[e]):
                continue
            idx = M_LOCAL_IDX[e]
            M0 = cum_forces[mid][idx]
            dM = df_by_member[mid][idx]
            lo, hi = hs.elastic_range(e)
            if dM > 0 and M0 + dM >= hi:
                events.append((float(np.clip((hi - M0) / dM, 0.0, 1.0)), mid, e, +1))
            elif dM < 0 and M0 + dM <= lo:
                events.append((float(np.clip((M0 - lo) / (-dM), 0.0, 1.0)), mid, e, -1))
    return events


def run_cyclic(frame, hinge_states, prescribed_dofs, direction, protocol, d_nominal,
               base_reaction_dofs, initial_cum_forces=None, max_steps=200000,
               rot_tol=1e-12, max_unload_iter=50):
    """位移控制的反覆載重。

    prescribed_dofs / direction: 跟 run_pushover 一樣(受控自由度與各自的比例); 第一個受控
        自由度是「控制點」, **protocol 是這個控制點的位移目標**(direction[0] 不能是 0)。
    protocol: 控制點位移目標序列, 例如 [0.02, -0.02, 0.04, -0.04, 0.0] 表示從 0 推到 +0.02、
        拉到 -0.02、推到 +0.04、拉到 -0.04、回到 0。相鄰目標之間依 d_nominal 分步。
    d_nominal: 每步位移增量上限(取絕對值); 降伏事件會自動在事件點切開, 不用為了抓事件調小。
    base_reaction_dofs: 底剪力 F = -Σ 這些自由度的反力。
    hinge_states: {member_id: CyclicHingeState}(可用 CyclicHingeState.from_hinge_states 轉換)。
    initial_cum_forces: 選用, 重力預載後的桿端內力起點(跟 run_pushover 相同意義)。
    max_unload_iter: 「卸載符號一致性」疊代上限, 超過就報錯(不給不一致的答案)。

    回傳 CyclicResult。
    """
    direction = np.asarray(direction, dtype=float)
    if direction[0] == 0.0:
        raise ValueError("direction[0] 是控制點的比例, 不能是 0")
    for mid, hs in hinge_states.items():
        if not isinstance(hs, CyclicHingeState):
            raise TypeError(f"member {mid} 的塑鉸狀態必須是 CyclicHingeState "
                            "(可用 CyclicHingeState.from_hinge_states(...) 轉換)")
    d_nominal = abs(float(d_nominal))
    fixed_dofs = _fixed_dof_set(frame)
    _, n_node_dof, n_extra_dof = build_dof_map(frame)
    n_total = n_node_dof + n_extra_dof
    control = prescribed_dofs[0]

    cum_forces = ({mid: np.zeros(6) for mid in frame.members} if initial_cum_forces is None
                  else {mid: np.array(f, dtype=float) for mid, f in initial_cum_forces.items()})
    cum_reaction = np.zeros(n_total)
    u_full = np.zeros(n_total)

    keys = [(mid, e) for mid in hinge_states for e in (0, 1)]
    hist_u, hist_F, hist_W = [0.0], [0.0], [0.0]
    hist_M = {k: [float(cum_forces[k[0]][M_LOCAL_IDX[k[1]]])] for k in keys}
    hist_tp = {k: [0.0] for k in keys}
    hist_work = {k: [0.0] for k in keys}
    events = []
    work_ext = 0.0
    step_count = 0
    amount = 0.0                                   # 目前「比例量」(控制點位移 = direction[0]*amount)

    def record(kind_events=()):
        nonlocal work_ext
        F = -float(sum(cum_reaction[d] for d in base_reaction_dofs))
        u = float(u_full[control])
        work_ext += 0.5 * (hist_F[-1] + F) * (u - hist_u[-1])
        hist_u.append(u); hist_F.append(F); hist_W.append(work_ext)
        for (mid, e) in keys:
            hs = hinge_states[mid]
            hist_M[(mid, e)].append(float(cum_forces[mid][M_LOCAL_IDX[e]]))
            hist_tp[(mid, e)].append(hs.theta_p_signed[e])
            hist_work[(mid, e)].append(hs.work[e])
        for kind, mid, e in kind_events:
            events.append({'kind': kind, 'member': mid, 'end': e, 'u': u, 'F': F, 'step': step_count})

    for target_u in protocol:
        target = float(target_u) / direction[0]
        while abs(target - amount) > 1e-12 * (1.0 + abs(target)):
            step_count += 1
            if step_count > max_steps:
                raise RuntimeError(f"反覆載重超過 {max_steps} 步仍未走完 protocol, 已中止(避免無窮迴圈)。"
                                   "請檢查 d_nominal 是否太小或模型是否有問題。")
            sgn = 1.0 if target > amount else -1.0
            d_step = sgn * min(d_nominal, abs(target - amount))

            # ---- (1) 用目前狀態試解, 疊代卸載的符號一致性 ----
            unloaded_now = []
            for _ in range(max_unload_iter):
                K, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(frame, hinge_states)
                try:
                    du = solve_displacement_increment(K, prescribed_dofs, direction * d_step, fixed_dofs)
                except RuntimeError as exc:
                    raise RuntimeError(f"第 {step_count} 步(控制點位移約 {hist_u[-1]:.6g})勁度矩陣奇異: "
                                       "結構在目前的塑鉸狀態下形成機構且不受位移控制約束。原始訊息: "
                                       + str(exc))
                changed = False
                for mid, hs in hinge_states.items():
                    if not any(hs.yielded):
                        continue
                    dth = _spring_rotation_increment(frame, hs, mid, member_T, member_dofs, member_L, du)
                    for e in (0, 1):
                        if hs.yielded[e] and hs.direction[e] * dth[e] < -rot_tol:
                            M0 = cum_forces[mid][M_LOCAL_IDX[e]]
                            hs.alpha[e] = M0 - hs.direction[e] * hs.Mp[e]     # 卸載點: 降伏面中心凍結
                            hs.yielded[e] = False
                            unloaded_now.append(('unload', mid, e))
                            changed = True
                if not changed:
                    break
            else:
                raise RuntimeError(f"第 {step_count} 步卸載的符號一致性疊代 {max_unload_iter} 次仍不收斂"
                                   "(多個塑鉸互相牽動), 已中止。")

            # ---- (2) 這一步內的降伏事件(彈性→塑性) ----
            df = _member_force_increments(frame, member_T, member_dofs, member_L, hinge_states, None, du)
            evs = _yield_events(hinge_states, cum_forces, df)
            ratio = min((r for r, *_ in evs), default=1.0)
            du_use = du * ratio
            d_used = d_step * ratio

            # ---- (3) 套用(線性: 增量對 ratio 成正比) ----
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
            cum_reaction += K @ du_use
            u_full += du_use
            amount += d_used

            newly = list(unloaded_now)
            if evs:
                for r, mid, e, s in evs:
                    if abs(r - ratio) < 1e-9:
                        hinge_states[mid].yielded[e] = True
                        hinge_states[mid].direction[e] = s
                        newly.append(('yield', mid, e))
            record(newly)

    return CyclicResult(
        u=np.array(hist_u), F=np.array(hist_F), work_ext=np.array(hist_W),
        hinge_M={k: np.array(v) for k, v in hist_M.items()},
        hinge_theta_p={k: np.array(v) for k, v in hist_tp.items()},
        hinge_work={k: np.array(v) for k, v in hist_work.items()},
        events=events, hinge_states=hinge_states, u_full=u_full, cum_forces=cum_forces,
        n_steps=step_count)


def make_protocol(amplitudes, n_cycles=1):
    """由「位移幅值序列」產生標準的反覆位移歷程: 每個幅值 a 走 n_cycles 個 (+a, -a) 完整循環,
    最後回到 0。例如 make_protocol([0.01, 0.02], 2) = [0.01,-0.01,0.01,-0.01, 0.02,-0.02,0.02,-0.02, 0]。"""
    proto = []
    for a in amplitudes:
        for _ in range(n_cycles):
            proto += [float(a), -float(a)]
    proto.append(0.0)
    return proto
