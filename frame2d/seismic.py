"""
非線性地震反應分析: 一站式整合入口 (動力分析 D8)。

把 D2(模態)、D5(Rayleigh阻尼+地震輸入)、D6(非線性時程+循環塑鉸)、選用的重力預載
(沿用 D7 的 `pushover.apply_gravity()`)串成一次呼叫, 加上能量平衡診斷、峰值反應摘要、
絕對加速度重建。這個檔案**不新增任何求解邏輯**, 純粹是把已經各自驗證過的元件組裝起來的
「一站式入口」, 跟 `cyclic.cyclic_analysis()`(D7 給網頁用的一站式入口)是同一種角色。

能量平衡是什麼、為什麼在地震輸入下一樣成立
--------------------------------------------
地震輸入的等效力 `P_eff(t) = -M r ag(t)`(D5 的 `ground_motion_force()`)本質上就是一個
「隨時間變化的外力」, 跟 D6 測試過的一般外力(諧和力等)沒有本質差異, 所以同一條能量平衡式

    外力作功 W_ext(t) = 動能 KE(t) + 阻尼耗能 W_damp(t) + 桿件內力作功 W_int(t)

(`W_int` 用桿件局部節點力對局部節點位移做功算, 同時涵蓋彈性儲能與塑性耗能, 見
`nonlinear_newmark.py` 模組說明)在地震輸入下一樣精確成立, 不需要另外推導。有重力預載時,
預載的靜力平衡狀態是動態擾動的起點(t=0 時 u=0, 但桿件內力已經是重力平衡值), 這個常數狀態
本身不影響**後續動態增量**的能量平衡(重力在純粹相對運動下不作淨功, 這是慣用假設, 這裡驗證
過確實不影響平衡式本身)。

單位: 跟輸入一致(質量=力·s²/長度, 跟D1~D7同一套)。
"""
from dataclasses import dataclass, field

import numpy as np

from .cyclic import CyclicHingeState
from .damping import rayleigh_damping_matrix
from .dofmanager import initial_hinge_states
from .excitation import absolute_acceleration as _absolute_acceleration
from .excitation import ground_motion_force
from .mass import assemble_M
from .modal import eigen
from .nonlinear_newmark import nonlinear_newmark_integrate
from .pushover import _assemble_stiffness_with_hinges, apply_gravity


@dataclass
class SeismicResult:
    """seismic_analysis() 的回傳值。大部分欄位直接轉發自內部的
    `nonlinear_newmark.NonlinearNewmarkResult`(`nl`), 這裡加的是地震分析特有的便利方法。"""
    nl: object                            # NonlinearNewmarkResult
    modal: object                         # eigen() 的結果(拿來算 Rayleigh 阻尼用的模態)
    direction: str
    ag: object = field(repr=False, default=None)          # 地面加速度 callable
    alpha: float = 0.0
    beta: float = 0.0
    zeta_target: float = 0.0
    mass_kind: str = 'lumped'
    gravity_result: object = field(repr=False, default=None)   # apply_gravity() 的靜力解, 沒做重力預載時是 None

    @property
    def t(self):
        return self.nl.t

    @property
    def u(self):
        return self.nl.u

    @property
    def v(self):
        return self.nl.v

    @property
    def a(self):
        return self.nl.a

    @property
    def frame(self):
        return self.nl.frame

    @property
    def events(self):
        return self.nl.events

    def peak_displacement(self, node, direction='x'):
        """某節點某方向(x/y/rot)相對位移的峰值(絕對值)。"""
        li = {'x': 0, 'y': 1, 'rot': 2}[direction]
        dof = self.frame.dofs_of(node)[li]
        return float(np.max(np.abs(self.u[:, dof])))

    def peak_drift_ratio(self, node_top, node_bottom, height, direction='x'):
        """兩節點之間相對位移峰值 ÷ 高度(層間位移角常用定義)。"""
        li = {'x': 0, 'y': 1}[direction]
        d1 = self.frame.dofs_of(node_top)[li]
        d2 = self.frame.dofs_of(node_bottom)[li]
        return float(np.max(np.abs(self.u[:, d1] - self.u[:, d2]))) / height

    def absolute_acceleration(self):
        """重建絕對加速度(樓層加速度/設備需求用), 見 `excitation.absolute_acceleration()`。"""
        return _absolute_acceleration(self.nl, self.direction, self.ag)

    def energy_balance(self):
        """算出四條能量時間序列(跟輸入單位一致的能量單位, 例如 N·m): {'KE','Wdamp','Wint',
        'Wext'}。物理上 Wext 應該精確等於 KE+Wdamp+Wint(見模組說明), 這個方法只負責算出來,
        不驗證——驗證在 tests/test_seismic.py。"""
        M = assemble_M(self.frame, self.mass_kind)
        KE = 0.5 * np.einsum('ti,ij,tj->t', self.v, M, self.v)
        Fser = np.array([self._force(tt) for tt in self.t])
        dWext = np.concatenate([[0.0], 0.5 * np.einsum('ti,ti->t', Fser[:-1] + Fser[1:], self.u[1:] - self.u[:-1])])
        Wext = np.cumsum(dWext)
        Cv = self.v @ self._C.T
        dWdamp = np.concatenate([[0.0], 0.5 * np.einsum('ti,ti->t', Cv[:-1] + Cv[1:], self.u[1:] - self.u[:-1])])
        Wdamp = np.cumsum(dWdamp)
        dWint = np.zeros(len(self.t))
        for mid in self.frame.members:
            idx = np.array(self._member_dofs[mid])
            T = self._member_T[mid]
            u_local = (T @ self.u[:, idx].T).T
            f_local = self.nl.member_force_history[mid]
            dWint += np.concatenate([[0.0], 0.5 * np.einsum('ti,ti->t', f_local[:-1] + f_local[1:],
                                                             u_local[1:] - u_local[:-1])])
        Wint = np.cumsum(dWint)
        return {'KE': KE, 'Wdamp': Wdamp, 'Wint': Wint, 'Wext': Wext}


def seismic_analysis(frame, ag, direction='x', dt=None, n_steps=None, zeta=0.05, damping_modes=(0, 2),
                     mass_kind='lumped', apply_gravity_loads=True) -> SeismicResult:
    """非線性地震反應分析(一站式入口)。

    ag: callable, ag(t) -> 地面加速度(跟模型一致的加速度單位, 例如 m/s²)。
    direction: 'x' 或 'y', 地面運動方向。
    dt, n_steps: 時間步長與步數; None 時分別預設 T1/50(T1=第一模態週期)跟涵蓋 10 秒的步數
        (只是合理的起點, 實務上應該自己依地震歷程長度跟精度需求設定)。
    zeta: Rayleigh 阻尼比目標值(兩個控制模態用同一個值)。
    damping_modes: 長度2的tuple, 指定用「模態分析排序後」第幾個模態(0起算)當 Rayleigh 阻尼的
        兩個控制頻率, 預設 (0,2) = 第1、第3模態(常見經驗法則)。
    mass_kind: 'lumped'(預設) 或 'consistent'。
    apply_gravity_loads: True(預設) 時, 如果模型上有點載重/分佈載重, 先用
        `pushover.apply_gravity()` 做一次靜力重力預載, 動態分析從這個平衡狀態的**相對**位移
        0 開始(見模組說明); False 或模型上沒有載重時跳過, 從無應力狀態開始。

    回傳 SeismicResult。
    """
    if direction not in ('x', 'y'):
        raise ValueError(f"direction必須是'x'或'y', 收到'{direction}'")
    hs0 = initial_hinge_states(frame)
    if not hs0:
        raise ValueError("模型裡沒有任何桿件設定塑鉸容量(Mp_i/Mp_j), 無法做非線性地震分析")

    md = eigen(frame, n_modes=None, mass=mass_kind)
    if max(damping_modes) >= md.n_modes:
        raise ValueError(f"damping_modes={damping_modes} 超出模態數(只有 {md.n_modes} 個模態)")
    wi, wj = md.omega[damping_modes[0]], md.omega[damping_modes[1]]
    C, alpha, beta = rayleigh_damping_matrix(frame, wi, wj, zeta, mass_kind=mass_kind)

    gravity_result = None
    initial_cum = None
    if apply_gravity_loads and (frame.point_loads or frame.distributed_loads or frame.member_point_loads):
        initial_cum, gravity_result = apply_gravity(frame, hs0)

    hs = CyclicHingeState.from_hinge_states(hs0)

    if dt is None:
        dt = md.period[0] / 50.0
    if n_steps is None:
        n_steps = int(round(10.0 / dt))

    force = ground_motion_force(frame, direction, ag, mass_kind=mass_kind)
    nl = nonlinear_newmark_integrate(frame, hs, dt, n_steps, force=force, mass_kind=mass_kind,
                                     damping_matrix=C, initial_cum_forces=initial_cum)

    res = SeismicResult(nl=nl, modal=md, direction=direction, ag=ag, alpha=alpha, beta=beta,
                        zeta_target=zeta, mass_kind=mass_kind, gravity_result=gravity_result)
    res._C = C
    _, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(frame, hs)
    res._member_dofs, res._member_T, res._member_L = member_dofs, member_T, member_L

    res._force = force
    return res


# ============================================================================
# 網頁 /nonlinear_seismic 端點用的一站式入口與地震歷程 adapter (動力分析 D6b)
# ============================================================================

MAX_WEB_STEPS = 4000      # 網頁一站式入口的步數上限, 避免使用者輸入的dt太小把後端卡住/回應過大


def sine_pulse_ground_motion(amplitude_g, freq_hz, decay=0.0):
    """簡化的衰減正弦波地震脈衝(**不是真實地震紀錄**, 示範/教學用): ag(t) = amplitude_g·g·
    sin(2π·freq_hz·t)·e^(-decay·t)。decay=0 就是不衰減的純諧和波。回傳 ag(t) callable(SI, m/s²)。
    """
    if not (amplitude_g > 0 and np.isfinite(amplitude_g)):
        raise ValueError(f"amplitude_g必須是正數, 收到{amplitude_g}")
    if not (freq_hz > 0 and np.isfinite(freq_hz)):
        raise ValueError(f"freq_hz必須是正數, 收到{freq_hz}")
    if decay < 0 or not np.isfinite(decay):
        raise ValueError(f"decay必須是不小於0的有限數字, 收到{decay}")
    amp = amplitude_g * 9.80665
    omega = 2 * np.pi * freq_hz
    return lambda t: amp * np.sin(omega * t) * np.exp(-decay * t)


def custom_ground_motion(points_g):
    """由一組 (t, ag) 資料點(ag 用 g 的倍數)做分段線性內插的地面加速度歷程。points_g:
    [(t0,ag0), (t1,ag1), ...], t 嚴格遞增、至少 2 點。範圍外(t < t0 或 t > 最後一點)夾在
    邊界值(不外插)。回傳 ag(t) callable(SI, m/s²)。"""
    pts = sorted((float(t), float(a)) for t, a in points_g)
    if len(pts) < 2:
        raise ValueError(f"自訂地震歷程至少需要 2 個資料點, 收到 {len(pts)} 個")
    ts = np.array([p[0] for p in pts])
    ags = np.array([p[1] for p in pts]) * 9.80665
    if np.any(ts < 0):
        raise ValueError("自訂地震歷程的時間 t 必須都是不小於0的數字")
    if np.any(np.diff(ts) <= 0):
        raise ValueError("自訂地震歷程的時間 t 必須嚴格遞增(不能有重複或反向的點)")
    return lambda t: float(np.interp(t, ts, ags))


@dataclass
class SeismicWebResult:
    """nonlinear_seismic_web_analysis() 的回傳值: SeismicResult + 給網頁畫圖/動畫用的取樣時間
    索引(`frame_idx`, 完整時間序列可能超過網頁想要傳輸/畫的點數, 用這組索引取子集)。"""
    seismic: object          # SeismicResult
    control_node: int
    frame_idx: np.ndarray    # 要傳給前端(動畫/圖表)的時間步索引, 遞增, 含0跟最後一步


def nonlinear_seismic_web_analysis(frame, control_node, direction='x', dt=None, n_steps=None,
                                   zeta=0.05, damping_modes=(0, 2), mass_kind='lumped',
                                   apply_gravity_loads=True, ground_motion_type='pulse',
                                   pulse_amplitude_g=None, pulse_freq_hz=None, pulse_decay=0.0,
                                   custom_points_g=None, peer_nga_text=None,
                                   peer_nga_max_points=2000, max_frames=400) -> SeismicWebResult:
    """網頁 /nonlinear_seismic 端點用的一站式入口: 建地震歷程 callable → `seismic_analysis()`
    → 取樣動畫用的時間索引(全部時間步太多時, 均勻抽稀到最多 max_frames 個, 首尾一定保留)。

    control_node: 動畫/時間歷程圖表要追蹤的節點(通常是頂層節點)。
    ground_motion_type: 'pulse'(`sine_pulse_ground_motion`)、'custom'(`custom_ground_motion`)、
        或 'peer_nga'(讀取 PEER NGA .AT2 格式的真實強震紀錄文字, 見
        `ground_motion_io.peer_nga_to_points()`; `peer_nga_max_points` 是抽稀上限, 真實紀錄
        常有上萬個取樣點, 直接全部塞進網頁傳輸不必要地肥大)。
    其餘參數見 `seismic_analysis()`。
    """
    if ground_motion_type == 'pulse':
        if pulse_amplitude_g is None or pulse_freq_hz is None:
            raise ValueError("ground_motion_type='pulse' 需要指定 pulse_amplitude_g 與 pulse_freq_hz")
        ag = sine_pulse_ground_motion(pulse_amplitude_g, pulse_freq_hz, pulse_decay)
    elif ground_motion_type == 'custom':
        if not custom_points_g:
            raise ValueError("ground_motion_type='custom' 需要指定 custom_points_g(至少2個(t, ag)點)")
        ag = custom_ground_motion(custom_points_g)
    elif ground_motion_type == 'peer_nga':
        if not peer_nga_text:
            raise ValueError("ground_motion_type='peer_nga' 需要指定 peer_nga_text(.AT2檔案的文字內容)")
        from .ground_motion_io import peer_nga_to_points
        ag = custom_ground_motion(peer_nga_to_points(peer_nga_text, max_points=peer_nga_max_points))
    else:
        raise ValueError(f"ground_motion_type必須是'pulse'、'custom'或'peer_nga', 收到'{ground_motion_type}'")

    if dt is not None and n_steps is not None and n_steps > MAX_WEB_STEPS:
        raise ValueError(f"n_steps={n_steps} 超過網頁上限 {MAX_WEB_STEPS}: 請加大 dt 或縮短總時間。")

    res = seismic_analysis(frame, ag, direction=direction, dt=dt, n_steps=n_steps, zeta=zeta,
                           damping_modes=damping_modes, mass_kind=mass_kind,
                           apply_gravity_loads=apply_gravity_loads)
    if res.nl.n_steps + 1 > MAX_WEB_STEPS:
        raise ValueError(f"這個設定總共要走 {res.nl.n_steps + 1} 個時間點, 超過網頁上限 {MAX_WEB_STEPS}: "
                         "請加大 dt 或縮短總時間(n_steps)。")
    try:
        frame.dofs_of(control_node)
    except KeyError:
        raise ValueError(f"找不到控制節點 {control_node}")

    n_total = res.nl.n_steps + 1
    if n_total <= max_frames:
        frame_idx = np.arange(n_total)
    else:
        frame_idx = np.unique(np.round(np.linspace(0, n_total - 1, max_frames)).astype(int))
    return SeismicWebResult(seismic=res, control_node=control_node, frame_idx=frame_idx)


def _finite(x):
    x = float(x)
    return x if np.isfinite(x) else None


def seismic_to_dict(pkg: SeismicWebResult) -> dict:
    """把 nonlinear_seismic_web_analysis() 的結果轉成可直接 json.dumps 的 dict。單位跟輸入一致
    (網頁後端固定 SI): 位移 m、力 N、彎矩 N·m、加速度 m/s²、時間 s、能量 N·m。

    只在 `frame_idx` 抽樣的時間點上傳輸節點位移(動畫用); 塑鉸 M/θp 曲線、時間歷程圖、能量
    平衡則傳完整的(未抽稀的)時間序列, 因為那些是靜態曲線圖, 資料量遠小於「每個時間步存全部
    節點位移」。"""
    res = pkg.seismic
    nl = res.nl
    frame = nl.frame
    fi = pkg.frame_idx
    li = {'x': 0, 'y': 1}[res.direction]
    ctrl_dof = frame.dofs_of(pkg.control_node)[li]

    frames = []
    for k in fi:
        nodes = {}
        for nid in frame.nodes:
            ux, uy, rot = frame.dofs_of(nid)
            nodes[str(nid)] = {'ux': _finite(nl.u[k, ux]), 'uy': _finite(nl.u[k, uy]), 'rot': _finite(nl.u[k, rot])}
        frames.append({'t': _finite(nl.t[k]), 'nodes': nodes})

    ag_series = [_finite(res.ag(tt)) for tt in nl.t]
    control_u = [_finite(v) for v in nl.u[:, ctrl_dof]]
    # 有效慣性力合力(等效地震力總和): 不是嚴格意義的支承反力(nonlinear_newmark 沒有算反力),
    # 是「這個方向上全部質量的慣性力合力」, 對單層結構等同基底剪力demand, 對多層結構是全部樓層
    # 合力(明確標註在前端文案裡, 不叫「基底剪力」)
    r = np.zeros(nl.n_dof)
    for nid in frame.nodes:
        r[frame.dofs_of(nid)[li]] = 1.0
    total_mass_dir = float((assemble_M(frame, res.mass_kind) @ r).sum())
    eq_force = [_finite(total_mass_dir * res.ag(tt)) for tt in nl.t]

    hinges = []
    for (mid, e), M_hist in nl.hinge_M.items():
        n_yield = sum(1 for ev in nl.events if ev['kind'] == 'yield' and ev['member'] == mid and ev['end'] == e)
        if n_yield == 0:
            continue
        hinges.append({
            'member': int(mid), 'end': int(e), 'label': f"M{mid} {'i' if e == 0 else 'j'}端",
            'M': [_finite(v) for v in M_hist], 'theta_p': [_finite(v) for v in nl.hinge_theta_p[(mid, e)]],
            'work_final': _finite(nl.hinge_work[(mid, e)][-1]), 'n_yield': n_yield,
        })
    hinges.sort(key=lambda h: -h['n_yield'])

    eb = res.energy_balance()
    events = [{'kind': e['kind'], 'member': int(e['member']), 'end': int(e['end']), 't': _finite(e['t'])}
             for e in nl.events]

    return {
        'analysis_type': 'nonlinear_seismic', 'direction': res.direction, 'control_node': pkg.control_node,
        'mass_kind': res.mass_kind, 'zeta_target': res.zeta_target, 'alpha': _finite(res.alpha),
        'beta': _finite(res.beta), 'period1': _finite(res.modal.period[0]), 'dt': _finite(nl.t[1] - nl.t[0]),
        'n_steps': int(nl.n_steps), 'apply_gravity_loads': res.gravity_result is not None,
        't': [_finite(v) for v in nl.t], 'ground_motion': ag_series, 'control_disp': control_u,
        'equivalent_force': eq_force,
        'peak_displacement': _finite(res.peak_displacement(pkg.control_node, res.direction)),
        'energy': {'t': [_finite(v) for v in nl.t], 'KE': [_finite(v) for v in eb['KE']],
                   'Wdamp': [_finite(v) for v in eb['Wdamp']], 'Wint': [_finite(v) for v in eb['Wint']],
                   'Wext': [_finite(v) for v in eb['Wext']]},
        'hinges': hinges, 'events': events,
        'frame_idx': [int(k) for k in fi], 'frames': frames,
    }
