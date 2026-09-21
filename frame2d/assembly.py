"""
assemble_K() -- 公開的全域勁度矩陣組裝入口(動力分析 D0)。

為什麼需要這個檔案
------------------
動力分析(質量矩陣、特徵值、時間積分)都需要「只有 K, 沒有載重跟邊界
條件」的乾淨全域勁度矩陣。這個東西其實早就存在了:
`pushover._assemble_stiffness_with_hinges()` 就是純粹的組裝函式, 不含
載重、不含邊界條件。缺的只是一個**公開的、名字不綁定pushover**的入口,
所以這裡只做一層很薄的包裝, 不重寫任何組裝邏輯。

刻意「只加不改」(見ROADMAP.md 動力分析路線的設計原則1)
------------------------------------------------------
- 沒有修改 `dofmanager._solve_once_dofmanager()`(它內部仍然有自己那份
  組裝迴圈, 跟這裡呼叫的是「兩份獨立寫的組裝」, 不是同一份程式碼)
- 沒有修改 `pushover._assemble_stiffness_with_hinges()`
- 兩份組裝會不會悄悄分歧, 由 tests/test_assembly.py 守住: 用 assemble_K()
  自己組K、自己解, 對照 solve() / solve_pdelta()核心 / solve_with_hinges()
  的結果, 加上獨立的剛體運動檢核(不依賴任何求解器)

已知範圍(動力分析初期刻意的限制, 見ROADMAP.md「既有限制」表格)
------------------------------------------------------------
- 不處理「鬆弛的cable」: cable一律當成taut(跟truss同一個軸向勁度)貢獻
  進K。要做cable鬆弛判斷的呼叫端得自己處理(動力分析D1~D5會直接拒絕
  含cable的模型, 之後再處理)
- equal_dofs 仍然是懲罰法(會把高勁度彈簧疊進K): 動力分析D1~D5會拒絕
  含equal_dofs的模型, 因為懲罰彈簧會製造虛假的高頻模態

用法:
    from frame2d.assembly import assemble_K
    asm = assemble_K(frame)
    asm.K              # (n_dof, n_dof) 全域勁度矩陣
    asm.member_dofs    # {member_id: (ux_i,uy_i,rot_i, ux_j,uy_j,rot_j)}
    asm.n_node_dof     # 節點自己的DOF數(3 * 節點數), 排在最前面
    asm.n_extra_dof    # release端專屬轉角DOF數, 排在節點DOF後面
"""
from typing import NamedTuple

import numpy as np

from .dofmanager import build_dof_map
from .model import Frame2D
from .pushover import _assemble_stiffness_with_hinges


class Assembly(NamedTuple):
    """assemble_K()的回傳值。NamedTuple, 所以既可以用欄位名稱取值
    (asm.K), 也可以直接拆包(K, member_dofs, member_T, member_L,
    n_node_dof, n_extra_dof = assemble_K(frame))。"""
    K: np.ndarray
    member_dofs: dict
    member_T: dict
    member_L: dict
    n_node_dof: int
    n_extra_dof: int

    @property
    def n_dof(self) -> int:
        return self.n_node_dof + self.n_extra_dof


def assemble_K(frame: Frame2D, hinge_states: dict = None, axial_forces: dict = None) -> Assembly:
    """組裝全域勁度矩陣K(不含載重、不含邊界條件)。

    hinge_states: 可選, {member_id: HingeState}, 該桿件改用含鉸勁度矩陣
        (跟solve_with_hinges()同一個意思)。
    axial_forces: 可選, {member_id: N(拉力為正)}, 疊加P-Delta幾何勁度
        (跟solve_pdelta()同一個意思)。

    兩個參數都是None(預設)就是最基本的線彈性K。
    """
    K, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(
        frame, hinge_states, axial_forces)
    _, n_node_dof, n_extra_dof = build_dof_map(frame)
    return Assembly(K, member_dofs, member_T, member_L, n_node_dof, n_extra_dof)
