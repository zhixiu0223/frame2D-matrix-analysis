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
