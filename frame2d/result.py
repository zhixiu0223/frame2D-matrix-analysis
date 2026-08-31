"""
求解結果的共用資料型別。solve.py(靜力凝縮)跟dofmanager.py(DOFManager)
兩套獨立的求解路徑都回傳這個格式, 這樣postprocess.py/plotting.py不用
關心底層是哪一種算法算出來的。
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class MemberResult:
    member_id: int
    L: float
    angle: float
    # 局部座標系桿端內力: [Fx1, Fy1, M1, Fx2, Fy2, M2]
    # Fx = 軸力(壓力為正的節點力慣例, postprocess.py取負號校正成拉力為正),
    # Fy = 剪力, M = 彎矩
    end_forces_local: np.ndarray
    slack: bool = False   # True = 這是一條cable, 且這次求解判定為鬆弛(不受力)


@dataclass
class SolveResult:
    displacements: np.ndarray          # 全域自由度位移向量, 長度 3*n_nodes
    reactions: np.ndarray              # 全域自由度反力向量 (僅支承自由度有意義)
    member_results: dict               # member_id -> MemberResult
    slack_cables: set = None           # 最終解裡判定為鬆弛的cable member_id集合
    frame: object = None               # 對應的Frame2D模型(給Result API查詢用,
                                        # 見query.py; 舊版solve()的呼叫端如果只用
                                        # displacements/reactions/member_results,
                                        # 不受這個新欄位影響)

    def member(self, member_id):
        """回傳這根桿件的便利查詢物件(見query.MemberQuery), 例如
        result.member(5).moment_at(2.3)。需要frame欄位(呼叫solve()時
        會自動帶入, 不需要使用者自己設定)。"""
        from .query import MemberQuery
        if self.frame is None:
            raise ValueError("這個SolveResult沒有附帶frame, 無法使用member()查詢"
                              "(直接呼叫solve()/solve_condensation()取得的結果"
                              "應該都會自動附帶, 如果是手動建構的SolveResult"
                              "請自行設定frame欄位)")
        return MemberQuery(self.frame, self, member_id)

    def max_moment(self):
        """全結構|M|最大值, 回傳ExtremeValue(value, member_id, x)。"""
        from .query import max_moment as _max_moment
        return _max_moment(self.frame, self)

    def max_shear(self):
        """全結構|V|最大值。"""
        from .query import max_shear as _max_shear
        return _max_shear(self.frame, self)

    def max_axial(self):
        """全結構|N|最大值。"""
        from .query import max_axial as _max_axial
        return _max_axial(self.frame, self)

    def max_displacement(self):
        """全結構最大變形量(沿桿件全長找偏離原始直線最遠的點, 不是只看
        節點, 定義跟plot_all()畫的Δmax標籤一致)。"""
        from .query import max_displacement as _max_displacement
        return _max_displacement(self.frame, self)
