"""
2D Frame 通用資料結構 (Node / Member / DOF)

設計原則:
- 使用者描述「結構」(節點+桿件+支承+載重),不用手寫勁度矩陣或自由度編號
- DOF 查詢包成 dof_index()/dofs_of()函式,不直接讓外部程式碼寫死 3*node_id,
  未來要換成真正的 DOFManager(支援 truss/不連續id/release)時,呼叫端不用改
"""
from dataclasses import dataclass, field


@dataclass
class Node:
    id: int
    x: float
    y: float


@dataclass
class Section:
    """材料/斷面性質,可被多根桿件共用"""
    name: str
    E: float   # 楊氏模數
    I: float   # 慣性矩
    A: float   # 斷面積 (若只做彎矩分析可給大數字近似軸向剛體)


@dataclass
class Member:
    """2D 樑柱元素:每端3個自由度(ux, uy, rot),共6個自由度"""
    id: int
    node_i: int      # 近端節點id
    node_j: int      # 遠端節點id
    section: str      # 對應 Section.name


@dataclass
class Support:
    node: int
    ux: bool = False   # True = 這個方向被拘束
    uy: bool = False
    rot: bool = False


@dataclass
class PointLoad:
    node: int
    fx: float = 0.0
    fy: float = 0.0
    m: float = 0.0


@dataclass
class DistributedLoad:
    """垂直於桿件局部y方向的均佈載重 (正值 = 沿局部+y方向)"""
    member: int
    w_start: float   # kN/m 或對應單位
    w_end: float = None  # None = 均佈 (w_end = w_start)

    def __post_init__(self):
        if self.w_end is None:
            self.w_end = self.w_start


class Frame2D:
    def __init__(self):
        self.nodes: dict[int, Node] = {}
        self.sections: dict[str, Section] = {}
        self.members: dict[int, Member] = {}
        self.supports: list[Support] = []
        self.point_loads: list[PointLoad] = []
        self.distributed_loads: list[DistributedLoad] = []

    # ---- 建模 API ----
    def add_node(self, id: int, x: float, y: float):
        self.nodes[id] = Node(id, x, y)
        return self

    def add_section(self, name: str, E: float, I: float, A: float = 1e8):
        self.sections[name] = Section(name, E, I, A)
        return self

    def add_member(self, id: int, node_i: int, node_j: int, section: str):
        self.members[id] = Member(id, node_i, node_j, section)
        return self

    def fix(self, node: int):
        self.supports.append(Support(node, ux=True, uy=True, rot=True))
        return self

    def pin(self, node: int):
        self.supports.append(Support(node, ux=True, uy=True, rot=False))
        return self

    def roller_y(self, node: int):
        """只拘束 uy (最常見的滾支承方向)"""
        self.supports.append(Support(node, ux=False, uy=True, rot=False))
        return self

    def point_load(self, node: int, fx: float = 0.0, fy: float = 0.0, m: float = 0.0):
        self.point_loads.append(PointLoad(node, fx, fy, m))
        return self

    def distributed_load(self, member: int, w: float, w_end: float = None):
        self.distributed_loads.append(DistributedLoad(member, w, w_end))
        return self

    # ---- DOF 查詢 (唯一允許碰自由度編號的地方) ----
    def n_dof(self) -> int:
        return 3 * len(self.nodes)

    def dofs_of(self, node_id: int) -> tuple[int, int, int]:
        """節點node_id的三個全域自由度編號: (ux, uy, rot)
        MVP: 假設節點 id 從 0 開始連續編號。若非連續,先在這裡建 id->index 對照表。
        """
        return (3 * node_id, 3 * node_id + 1, 3 * node_id + 2)
