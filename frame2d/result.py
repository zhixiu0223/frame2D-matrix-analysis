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
