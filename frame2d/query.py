"""
Result API -- 便利查詢介面

在已經驗證過的 postprocess.py 底層函式(member_internal_forces,
member_deformed_shape)之上包一層方便查詢的介面, 純粹是介面整理, 不改變
任何計算邏輯——所有數字都還是同一套已驗證公式算出來的, 只是不用再自己
手動呼叫 member_internal_forces + np.interp 去查特定位置的內力。

用法:
    result = solve(frame)
    result.max_moment()              # 全結構彎矩最大值+發生位置
    result.member(5).moment_at(2.3)  # 5號桿件在x=2.3m處的彎矩
    result.member(5).max_moment()    # 5號桿件自己的彎矩最大值+位置

「最大變形」(max_displacement)的定義: 沿桿件全長掃描, 找出偏離原始
直線最遠的點(=plotting.py畫Δmax標籤用的同一套邏輯, 兩者數字保證一致),
不是只看節點位移——因為桿件跨中在有載重時的撓度往往比兩端節點位移還大,
只看節點會低估真正的最大變形。
"""
from dataclasses import dataclass
import numpy as np

from .elements import member_geometry
from .postprocess import member_internal_forces, member_deformed_shape


@dataclass
class ExtremeValue:
    """一次極值查詢的結果: 數值 + 發生在哪根桿件的哪個局部座標位置。"""
    value: float
    member_id: int
    x: float   # 沿桿件局部座標(0到該桿件長度L), 不是全域座標

    def __repr__(self):
        return f"ExtremeValue(value={self.value:.4g}, member_id={self.member_id}, x={self.x:.4g})"


class MemberQuery:
    """result.member(member_id) 回傳的物件, 提供對單一桿件內力/變形的
    便利查詢。內部直接呼叫已驗證過的 member_internal_forces /
    member_deformed_shape, 用np.interp在查詢點取值(桿件內部有集中力/
    局部段均佈載重時, 底層取樣點會在那些位置加密, 直接查任意x都準確,
    不會像早期版本一樣因為取樣點疏密不均對錯位置)。"""

    def __init__(self, frame, result, member_id, n_sample=201):
        self._frame = frame
        self._result = result
        self._member_id = member_id
        self._n_sample = n_sample
        self._cache_forces = None
        self._cache_defl = None

    def _forces(self):
        if self._cache_forces is None:
            self._cache_forces = member_internal_forces(
                self._frame, self._result, self._member_id, n=self._n_sample)
        return self._cache_forces

    def _deflection_offset(self):
        """回傳(x, offset): x是沿桿件局部座標, offset是該點偏離原始直線
        的距離(跟plotting.py畫Δmax標籤用的同一套算法, 數字保證一致)。"""
        if self._cache_defl is None:
            m = self._frame.members[self._member_id]
            ni = self._frame.nodes[m.node_i]
            nj = self._frame.nodes[m.node_j]
            L, angle = member_geometry(ni, nj)
            n = self._n_sample
            t = np.linspace(0, L, n)
            base_x = ni.x + t * np.cos(angle)
            base_y = ni.y + t * np.sin(angle)
            X, Y = member_deformed_shape(self._frame, self._result, self._member_id,
                                          scale=1.0, n=n)
            offset = np.hypot(X - base_x, Y - base_y)
            self._cache_defl = (t, offset)
        return self._cache_defl

    @property
    def L(self):
        return self._result.member_results[self._member_id].L

    def moment_at(self, x):
        """回傳這根桿件在局部座標x(公尺, 從node_i算起)處的彎矩。"""
        x_full, N, V, M = self._forces()
        return float(np.interp(x, x_full, M))

    def shear_at(self, x):
        """回傳這根桿件在局部座標x處的剪力。"""
        x_full, N, V, M = self._forces()
        return float(np.interp(x, x_full, V))

    def axial_at(self, x):
        """回傳這根桿件在局部座標x處的軸力(拉力為正)。"""
        x_full, N, V, M = self._forces()
        return float(np.interp(x, x_full, N))

    def deflection_at(self, x):
        """回傳這根桿件在局部座標x處偏離原始直線的距離(=撓度大小)。"""
        t, offset = self._deflection_offset()
        return float(np.interp(x, t, offset))

    def max_moment(self):
        x_full, N, V, M = self._forces()
        i = int(np.argmax(np.abs(M)))
        return ExtremeValue(value=float(M[i]), member_id=self._member_id, x=float(x_full[i]))

    def max_shear(self):
        x_full, N, V, M = self._forces()
        i = int(np.argmax(np.abs(V)))
        return ExtremeValue(value=float(V[i]), member_id=self._member_id, x=float(x_full[i]))

    def max_axial(self):
        x_full, N, V, M = self._forces()
        i = int(np.argmax(np.abs(N)))
        return ExtremeValue(value=float(N[i]), member_id=self._member_id, x=float(x_full[i]))

    def max_deflection(self):
        t, offset = self._deflection_offset()
        i = int(np.argmax(offset))
        return ExtremeValue(value=float(offset[i]), member_id=self._member_id, x=float(t[i]))


def _scan_all_members(frame, result, extractor):
    """對全部桿件呼叫extractor(取得一個ExtremeValue), 回傳|value|最大的
    那一個。跳過鬆弛的cable(slack=True), 因為那些桿件已知不受力,
    掃它們的內力/變形沒有意義。"""
    best = None
    for mid, m in frame.members.items():
        mr = result.member_results.get(mid)
        if mr is not None and mr.slack:
            continue
        ev = extractor(MemberQuery(frame, result, mid))
        if best is None or abs(ev.value) > abs(best.value):
            best = ev
    return best


def max_moment(frame, result):
    """全結構掃描, 回傳|M|最大值(數值+發生在哪根桿件的哪個位置)。"""
    return _scan_all_members(frame, result, lambda q: q.max_moment())


def max_shear(frame, result):
    """全結構掃描, 回傳|V|最大值。"""
    return _scan_all_members(frame, result, lambda q: q.max_shear())


def max_axial(frame, result):
    """全結構掃描, 回傳|N|最大值。"""
    return _scan_all_members(frame, result, lambda q: q.max_axial())


def max_displacement(frame, result):
    """全結構掃描, 回傳偏離原始直線最遠的點(=plotting.py畫Δmax用的同一套
    定義, 沿桿件全長找, 不是只看兩端節點——桿件跨中在有載重時的撓度
    往往比節點位移還大)。"""
    return _scan_all_members(frame, result, lambda q: q.max_deflection())
