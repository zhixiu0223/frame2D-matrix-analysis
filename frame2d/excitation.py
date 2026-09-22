"""
外力時間歷程的小工具 (動力分析 D4)。

時程分析的外力慣例是「空間分佈 × 時間函式」: `F(t) = p * f(t)`, `p` 是固定的空間力向量
(哪個節點的哪個方向承受多少比例的力), `f(t)` 是純量的時間函式(諧和、脈衝、…)。這個檔案
提供兩邊各自的建構工具, 組出 `newmark.newmark_integrate()` 要的 `force(t) -> 長度n_dof向量`
callable。

不確定要抓哪個時間點的力最好算, 因為 Newmark 積分本來就是逐步走過去的, 不需要抓尖點——
只要 `f(t)` 在整個時間範圍內是良好定義的函式即可。
"""
import numpy as np


def point_force_vector(frame, entries):
    """把 {node: (fx, fy, m)} 或 [(node, direction, amount), ...] 轉成長度 n_dof 的空間力向量
    (跟 `assembly.assemble_K` 用同一套 DOF 編號; 呼叫端通常會拿 `newmark_integrate` 回傳的
    `NewmarkResult.n_dof` 或自己呼叫一次 `assembly.assemble_K` 先取得 DOF 數)。

    entries 可以是:
      - {node_id: (fx, fy, m)}: 每個節點在該處的 x向力、y向力、彎矩(單位比例, 由外層的時間
        函式決定實際大小), 缺的節點視為(0,0,0)
      - [(node_id, 'x'|'y'|'rot', amount), ...]: 逐項指定, 同一個 (node, direction) 出現
        多次會加總
    """
    from .assembly import assemble_K
    asm = assemble_K(frame)
    p = np.zeros(asm.n_dof)
    if isinstance(entries, dict):
        for node, (fx, fy, m) in entries.items():
            ux, uy, rot = frame.dofs_of(node)
            p[ux] += fx
            p[uy] += fy
            p[rot] += m
    else:
        li = {'x': 0, 'y': 1, 'rot': 2}
        for node, direction, amount in entries:
            if direction not in li:
                raise ValueError(f"方向必須是 'x'、'y' 或 'rot', 收到 '{direction}'")
            p[frame.dofs_of(node)[li[direction]]] += amount
    return p


def harmonic(amplitude, omega, phase=0.0):
    """f(t) = amplitude * sin(omega*t + phase)。"""
    return lambda t: amplitude * np.sin(omega * t + phase)


def step(amplitude, t_start=0.0):
    """f(t) = amplitude, t>=t_start; 0, t<t_start(階躍載重)。"""
    return lambda t: amplitude if t >= t_start else 0.0


def pulse(amplitude, t_start, duration):
    """f(t) = amplitude, t_start<=t<t_start+duration; 0 其他時候(矩形脈衝)。"""
    if duration <= 0:
        raise ValueError(f"duration必須是正數, 收到{duration}")
    return lambda t: amplitude if t_start <= t < t_start + duration else 0.0


def ramp(amplitude, t_start, t_rise):
    """f(t): 從 t_start 開始用 t_rise 秒線性升到 amplitude, 之後維持 amplitude; t<t_start 時為0。"""
    if t_rise <= 0:
        raise ValueError(f"t_rise必須是正數, 收到{t_rise}")

    def f(t):
        if t < t_start:
            return 0.0
        if t < t_start + t_rise:
            return amplitude * (t - t_start) / t_rise
        return amplitude
    return f


def combine(*time_functions):
    """把多個純量時間函式加總成一個(例如同時有階躍載重+諧和微擾)。"""
    return lambda t: sum(f(t) for f in time_functions)


def force_series_from_pattern(spatial_vector, time_function):
    """組出 `newmark_integrate` 要的 force callable: F(t) = spatial_vector * time_function(t)。"""
    spatial_vector = np.asarray(spatial_vector, dtype=float)
    return lambda t: spatial_vector * time_function(t)
