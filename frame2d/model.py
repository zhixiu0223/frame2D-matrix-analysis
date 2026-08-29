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
    """2D 桿件元素。member_type='frame'(預設): 樑柱元素, 每端3個自由度
    (ux,uy,rot), 可傳軸力+剪力+彎矩。member_type='truss': 桁架元素, 兩端鉸接,
    只能傳軸力(彎曲/剪力勁度為0), 用同一組6自由度格式儲存(v,theta永遠是0),
    方便共用組裝程式碼, 不用另外做一套DOF系統。
    release_i/release_j: frame元素專用, 該端是否有內部鉸接(彎矩釋放M=0)。
    用靜力凝縮處理, 不改變DOF系統, 不用切節點就能表示"桿件端點是鉸接,
    不是剛接"這件事。"""
    id: int
    node_i: int      # 近端節點id
    node_j: int      # 遠端節點id
    section: str      # 對應 Section.name
    member_type: str = 'frame'   # 'frame' 或 'truss' 或 'cable'
    release_i: bool = False
    release_j: bool = False


@dataclass
class Support:
    """支承: None=該方向自由, 數值=該方向指定位移(0.0=固定在原位,
    非0=強制位移/沉陷分析)。統一了fixed/pin/roller/settlement,
    不用另外做一個SupportDisplacement類別。"""
    node: int
    ux: float = None
    uy: float = None
    rot: float = None


@dataclass
class PointLoad:
    node: int
    fx: float = 0.0
    fy: float = 0.0
    m: float = 0.0


@dataclass
class DistributedLoad:
    """垂直於桿件局部y方向的均佈載重 (正值 = 沿局部+y方向)。
    預設整根桿件都有(x_start=None -> 0, x_end=None -> 桿件全長L);
    可以指定x_start/x_end只加在桿件的局部一段(0<=x_start<=x_end<=L)。"""
    member: int
    w_start: float   # kN/m 或對應單位
    w_end: float = None  # None = 均佈 (w_end = w_start)
    x_start: float = None   # None = 0 (從node_i開始)
    x_end: float = None     # None = 桿件全長 (到node_j為止)

    def __post_init__(self):
        if self.w_end is None:
            self.w_end = self.w_start


@dataclass
class MemberPointLoad:
    """桿件內部任意位置(不一定在節點上)的集中力/集中力矩。
    a: 距node_i沿桿軸的距離(局部座標, 0<=a<=L)。
    fx: 沿局部+x方向(軸向)的力; fy: 沿局部+y方向(橫向)的力,
    跟distributed_load的w同一套正負號慣例; m: 逆時針為正的集中力矩。
    不會新增節點、不切割桿件——用等效節點載重(work-equivalent)處理,
    跟distributed_load同一套機制(見solve.py)。
    """
    member: int
    a: float
    fx: float = 0.0
    fy: float = 0.0
    m: float = 0.0


class Frame2D:
    def __init__(self):
        self.nodes: dict[int, Node] = {}
        self.sections: dict[str, Section] = {}
        self.members: dict[int, Member] = {}
        self.supports: list[Support] = []
        self.point_loads: list[PointLoad] = []
        self.distributed_loads: list[DistributedLoad] = []
        self.member_point_loads: list[MemberPointLoad] = []
        self._node_index_cache: dict[int, int] = None   # node_id -> 緊湊的0-based索引, 延遲建立

    # ---- 建模 API ----
    def add_node(self, id: int, x: float, y: float):
        self.nodes[id] = Node(id, x, y)
        self._node_index_cache = None   # 節點集合變了, 快取失效
        return self

    def add_section(self, name: str, E: float, I: float, A: float = 1e8):
        self.sections[name] = Section(name, E, I, A)
        return self

    def add_member(self, id: int, node_i: int, node_j: int, section: str, member_type: str = 'frame',
                   release_i: bool = False, release_j: bool = False):
        self.members[id] = Member(id, node_i, node_j, section, member_type, release_i, release_j)
        return self

    def add_truss(self, id: int, node_i: int, node_j: int, section: str):
        """桁架元素的簡寫: 等同 add_member(..., member_type='truss')。
        兩端視為鉸接, 只傳軸力。可承受拉力或壓力(例如撐架的撐桿)。"""
        self.members[id] = Member(id, node_i, node_j, section, member_type='truss')
        return self

    def add_cable(self, id: int, node_i: int, node_j: int, section: str):
        """纜線元素的簡寫: 等同 add_member(..., member_type='cable')。
        跟truss一樣兩端鉸接、只傳軸力, 但只能受拉(壓力=0, 物理上代表纜線鬆弛
        退出作用)。solve()會自動偵測受壓的纜線、移除其勁度貢獻、重新求解,
        反覆直到沒有纜線受壓為止(見 solve.py 的說明)。"""
        self.members[id] = Member(id, node_i, node_j, section, member_type='cable')
        return self

    def fix(self, node: int):
        self.supports.append(Support(node, ux=0.0, uy=0.0, rot=0.0))
        return self

    def pin(self, node: int):
        self.supports.append(Support(node, ux=0.0, uy=0.0, rot=None))
        return self

    def roller_y(self, node: int):
        """只拘束 uy (最常見的滾支承方向)"""
        self.supports.append(Support(node, ux=None, uy=0.0, rot=None))
        return self

    def support(self, node: int, ux: float = None, uy: float = None, rot: float = None):
        """通用支承: None=該方向自由, 數值=該方向指定位移
        (0.0=固定在原位, 非0=強制位移/沉陷)。fix/pin/roller_y是這個的
        簡寫, 要做斜支承、沉陷分析等特殊情況直接用這個。"""
        self.supports.append(Support(node, ux=ux, uy=uy, rot=rot))
        return self

    def point_load(self, node: int, fx: float = 0.0, fy: float = 0.0, m: float = 0.0):
        self.point_loads.append(PointLoad(node, fx, fy, m))
        return self

    def distributed_load(self, member: int, w: float, w_end: float = None,
                          x_start: float = None, x_end: float = None):
        """均佈/線性變化載重。預設(x_start=x_end=None)整根桿件都有;
        指定x_start/x_end可以只加在桿件的局部一段(局部座標, 0<=x_start<=x_end<=L)。"""
        self.distributed_loads.append(DistributedLoad(member, w, w_end, x_start, x_end))
        return self

    def member_point_load(self, member: int, a: float, fx: float = 0.0, fy: float = 0.0, m: float = 0.0):
        """桿件內部任意位置(距node_i為a)加集中力/集中力矩, 不新增節點。
        跟point_load(node,...)的差別: 這裡a是桿件"局部"座標(沿桿軸距離node_i
        多遠), 不是節點id。"""
        self.member_point_loads.append(MemberPointLoad(member, a, fx, fy, m))
        return self

    # ---- DOF 查詢 (唯一允許碰自由度編號的地方) ----
    def n_dof(self) -> int:
        return 3 * len(self.nodes)

    def node_index(self, node_id: int) -> int:
        """node_id -> 緊湊的0-based索引(照節點加入順序編號, 不要求node_id本身
        連續或從0開始)。這是dofs_of()底層真正的id->index對照表, 延遲建立+
        快取(加新節點時失效重建)。"""
        if self._node_index_cache is None:
            self._node_index_cache = {nid: i for i, nid in enumerate(self.nodes.keys())}
        return self._node_index_cache[node_id]

    def dofs_of(self, node_id: int) -> tuple[int, int, int]:
        """節點node_id的三個全域自由度編號: (ux, uy, rot)。
        node_id不需要連續或從0開始(內部透過node_index()對照到緊湊索引),
        例如節點id用10,25,99也只會佔用9個DOF, 不會浪費空間到300個。"""
        i = self.node_index(node_id)
        return (3 * i, 3 * i + 1, 3 * i + 2)
