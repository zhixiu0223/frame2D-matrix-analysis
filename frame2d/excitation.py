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


def ground_motion_force(frame, direction, ag, mass_kind='lumped'):
    """地面加速度輸入(動力分析 D5), 回傳等效地震力 callable P_eff(t) = -M r ag(t), 可以直接
    傳給 `newmark.newmark_integrate(force=...)`。

    direction: 'x' 或 'y', 地面運動的方向。ag: callable, ag(t) -> 純量地面加速度(單位跟模型
        的加速度一致, 例如 m/s²; 如果你的地震歷程是用 g 的倍數記錄的, 記得先乘上重力加速度
        常數再傳進來)。mass_kind: 跟後面呼叫 `newmark_integrate()` 用的要一樣('lumped' 或
        'consistent')。

    這個力天生就滿足 `newmark_integrate()` 「力不能直接施加在無質量DOF上」的要求, 不需要
    另外檢查: 質量矩陣是半正定的, 對角線是0的DOF整列/整行必為0(半正定矩陣的性質), 所以
    M @ r 在無質量DOF上自動是0。

    只回傳「相對」運動方程式的力(標準做法); 積分完的 `NewmarkResult.u/v/a` 都是相對於地面的
    相對值。要重建絕對加速度(例如算樓層反應譜、或設備的加速度需求), 用
    `absolute_acceleration()`。
    """
    if direction not in ('x', 'y'):
        raise ValueError(f"direction必須是'x'或'y', 收到'{direction}'")
    from .mass import assemble_M, influence_vector
    M = assemble_M(frame, mass_kind)          # 同時做動力分析的模型檢查
    r = influence_vector(frame, direction)
    Mr = M @ r
    return lambda t: -Mr * ag(t)


def absolute_acceleration(result, direction, ag):
    """把 `newmark_integrate()` 算出來的相對加速度轉成絕對加速度: a_abs(t) = a_rel(t) + r·ag(t)。
    只有在做地震輸入(`ground_motion_force()`)時才有意義; 一般外力(諧和力、脈衝力等)不要用
    這個函式, 因為那種情況下 a_rel 本來就是絕對加速度, 沒有「地面也在動」這件事。

    result: `newmark.newmark_integrate()` 的回傳值。direction: 'x'或'y', 要跟算 result 時用
    的 `ground_motion_force()` 方向一致。ag: 跟算 result 時用的同一個地面加速度函式。

    回傳 (n_steps+1, n_dof) 陣列, 每個節點在該方向的分量是 a_rel + ag(t); 其他方向的分量
    (跟輸入方向垂直的平動、轉角)不受地面運動影響, 直接等於相對加速度(因為影響向量 r 在那些
    分量上是0)。
    """
    if direction not in ('x', 'y'):
        raise ValueError(f"direction必須是'x'或'y', 收到'{direction}'")
    from .mass import influence_vector
    r = influence_vector(result.frame, direction)
    ag_series = np.array([ag(tt) for tt in result.t])
    return result.a + np.outer(ag_series, r)
