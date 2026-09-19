"""
2D Frame 通用資料結構 (Node / Member / DOF)

設計原則:
- 使用者描述「結構」(節點+桿件+支承+載重),不用手寫勁度矩陣或自由度編號
- DOF 查詢包成 dof_index()/dofs_of()函式,不直接讓外部程式碼寫死 3*node_id,
  未來要換成真正的 DOFManager(支援 truss/不連續id/release)時,呼叫端不用改
"""
from dataclasses import dataclass, field
import math


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
    不是剛接"這件事。
    Mp_i/Mp_j/R_post_yield_i/R_post_yield_j: frame元素專用, 選用的塑性鉸
    容量(見hinge.py的HingeState)。這是"這根桿件材料能承受多少"的結構性質,
    跟release不同的地方是release是"永遠沒有勁度", 塑鉸是"降伏前有勁度、
    降伏後換一個較軟的勁度"——所以這裡只存"容量"這個靜態資料, 不存"目前
    有沒有降伏"這個分析過程中才會變化的狀態(那個狀態活在HingeState實例
    裡, 由呼叫端在求解時自己建立/傳入, 不跟著Member這個結構描述走)。
    預設None時完全不影響任何現有行為, 跟標準彈性frame元素一樣。"""
    id: int
    node_i: int      # 近端節點id
    node_j: int      # 遠端節點id
    section: str      # 對應 Section.name
    member_type: str = 'frame'   # 'frame' 或 'truss' 或 'cable'
    release_i: bool = False
    release_j: bool = False
    Mp_i: float = None
    Mp_j: float = None
    R_post_yield_i: float = None
    R_post_yield_j: float = None


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
class EqualDOF:
    """跟OpenSeesPy的equalDOF同一個概念: 讓slave_node指定的自由度
    強制跟master_node對應的自由度相等(不是"連在一起變成同一個節點",
    是"這幾個自由度的值被綁定相等", 其餘沒指定的自由度還是各自獨立)。
    ux/uy/rot各自是bool, True=這個方向要綁定, False=這個方向不受
    影響、各自獨立。用高勁度彈簧懲罰法實現, 見dofmanager.py的
    _apply_equal_dof()說明——不是精確的自由度消去法(那個要動到
    _solve_once_dofmanager()的組裝核心, 風險較高), 懲罰法用一個
    遠大於結構本身勁度的虛擬彈簧硬把兩個自由度"拉在一起", 精度
    足夠工程使用(見測試案例的驗證), 但嚴格來說不是完全精確為0
    的束制, 是"非常接近"。"""
    master_node: int
    slave_node: int
    ux: bool = False
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
    """垂直於桿件局部y方向的均佈載重 (正值 = 沿局部+y方向)。
    預設整根桿件都有(x_start=None -> 0, x_end=None -> 桿件全長L);
    可以指定x_start/x_end只加在桿件的局部一段(0<=x_start<=x_end<=L)。

    direction='local'(預設, 既有行為): w就是局部+y方向的分量, 跟桿件
    本身的角度無關。
    direction='global_y': w代表"沿全域垂直方向, 大小以沿桿件長度量測"
    的均佈載重(例如屋頂重力/雪載重的標準表示方式: 不管桿件本身斜不斜,
    每公尺桿長多重, 方向永遠垂直向下)。w_start!=w_end(線性變化)已支援,
    但目前仍只支援整根桿件(不支援x_start/x_end局部段)。
    direction='global': 全域"任意角度"均佈載重(global_y的推廣版, 用
    angle_deg指定角度, global_y等同angle_deg=-90或270)。w的大小以沿
    桿件長度量測(跟global_y同一套慣例), 方向固定是全域座標下角度
    angle_deg(標準數學慣例, 0度=+x方向, 逆時針為正, 跟桿件本身的角度
    無關, 不管桿件是斜的或水平的)。跟global_y不同, 這個方向**支援
    局部段(x_start/x_end)跟線性變化(w_start!=w_end)的任意組合**——
    推導: 桿件是直的, 角度沿桿長不變, 所以全域載重向量投影到局部x/y
    座標後仍然各自是線性函數, 用跟fixed_end_forces_partial_udl()同一套
    高斯積分(對point_load公式積分, 數值精確不是近似)分開處理局部x
    (fixed_end_forces_axial_partial_udl)、局部y(既有的
    fixed_end_forces_partial_udl)兩個分量即可。用途: 例如SW FEA app
    儲存的載重角度跟桿件本身角度不完全對齊(local+90或local-90取決於
    使用者畫圖方向)時, 直接用app紀錄的絕對角度重現, 不用去猜app的
    local慣例。"""
    member: int
    w_start: float   # kN/m 或對應單位
    w_end: float = None  # None = 均佈 (w_end = w_start)
    x_start: float = None   # None = 0 (從node_i開始)
    x_end: float = None     # None = 桿件全長 (到node_j為止)
    direction: str = 'local'   # 'local'、'global_y' 或 'global'
    angle_deg: float = None   # 只有direction='global'時使用: 全域角度(度)

    def __post_init__(self):
        if self.w_end is None:
            self.w_end = self.w_start
        if self.direction == 'global_y' and (self.x_start is not None or self.x_end is not None):
            raise ValueError(
                "direction='global_y' 目前只支援整根桿件, 不支援局部段"
                "(x_start/x_end); 若需要局部段+任意角度, 改用direction='global'")
        if self.direction == 'global' and self.angle_deg is None:
            raise ValueError("direction='global' 必須指定angle_deg(全域角度, 度)")


@dataclass
class DistributedMoment:
    """桿件全長均佈的分布彎矩(kN·m/m, 逆時針為正, 跟MemberPointLoad.m
    /fixed_end_forces_point_moment()同一個符號慣例——業界對「均佈彎矩」
    的正負號沒有統一標準, 這裡刻意沿用本模組既有的集中力矩慣例, 求
    內部一致)。

    目前只支援整根桿件的常數m(不支援線性變化或局部段, 比
    DistributedLoad簡化很多)——這是刻意先做最常見/最基本的情況(見
    Frame2D.distributed_moment()的docstring說明實務上什麼情況會用到
    這個載重類型), 不是忘記做完整版, 之後有需要可以再擴充。"""
    member: int
    m: float   # kN·m/m或對應單位, 逆時針為正


@dataclass
class MemberPointLoad:
    """桿件內部任意位置(不一定在節點上)的集中力/集中力矩。
    a: 距node_i沿桿軸的距離(局部座標, 0<=a<=L)。
    fx: 沿局部+x方向(軸向)的力; fy: 沿局部+y方向(橫向)的力,
    跟distributed_load的w同一套正負號慣例; m: 逆時針為正的集中力矩。
    不會新增節點、不切割桿件——用等效節點載重(work-equivalent)處理,
    跟distributed_load同一套機制(見solve.py)。

    direction='local'(預設, 既有行為): fx/fy就是局部座標分量。
    direction='global': 用F(大小)+angle_deg(全域角度, 標準數學慣例,
    0度=+x, 逆時針為正)指定一個全域方向的集中力, 依桿件角度自動分解
    成局部fx/fy(跟distributed_load(direction='global')同一套邏輯,
    只是這裡是點載重不是分佈載重)。這個模式下fx/fy欄位不使用(必須
    是預設值0, 由F/angle_deg決定), m仍然照常填(力矩沒有方向性問題)。
    只有主要求解器solve()支援, solve_condensation()遇到會直接報錯。
    """
    member: int
    a: float
    fx: float = 0.0
    fy: float = 0.0
    m: float = 0.0
    direction: str = 'local'   # 'local' 或 'global'
    F: float = None            # 只有direction='global'時使用: 力的大小
    angle_deg: float = None    # 只有direction='global'時使用: 全域角度(度)

    def __post_init__(self):
        if self.direction == 'global':
            if self.F is None or self.angle_deg is None:
                raise ValueError("direction='global' 必須指定F(力的大小)跟angle_deg(全域角度, 度)")
            if self.fx != 0.0 or self.fy != 0.0:
                raise ValueError("direction='global' 時fx/fy不使用, 請改用F+angle_deg指定")


class Frame2D:
    def __init__(self):
        self.nodes: dict[int, Node] = {}
        self.sections: dict[str, Section] = {}
        self.members: dict[int, Member] = {}
        self.supports: list[Support] = []
        self.equal_dofs: list[EqualDOF] = []
        self.point_loads: list[PointLoad] = []
        self.distributed_loads: list[DistributedLoad] = []
        self.distributed_moments: list[DistributedMoment] = []
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
                   release_i: bool = False, release_j: bool = False,
                   Mp_i: float = None, Mp_j: float = None,
                   R_post_yield_i: float = None, R_post_yield_j: float = None):
        self.members[id] = Member(id, node_i, node_j, section, member_type, release_i, release_j,
                                   Mp_i, Mp_j, R_post_yield_i, R_post_yield_j)
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

    def equal_dof(self, master_node: int, slave_node: int,
                  ux: bool = False, uy: bool = False, rot: bool = False):
        """跟OpenSeesPy的equalDOF同一個概念: 讓slave_node指定的自由度
        (ux/uy/rot, 每個True/False獨立設定)強制等於master_node的對應
        自由度——常見用途: 剛性樓板(多個節點的水平位移綁在一起)、
        剛性連桿(兩個重疊節點的部分自由度綁死)、鉸接處只放開轉角。
        用高勁度彈簧懲罰法實現(見frame2d.model.EqualDOF、dofmanager.py
        的_apply_equal_dof()), 不是完全精確為0的束制, 但精度足夠工程
        使用。至少要指定一個方向為True, 不然這個約束沒有意義。"""
        if not (ux or uy or rot):
            raise ValueError("equal_dof()至少要指定ux/uy/rot其中一個為True, 不然這個約束沒有意義。")
        self.equal_dofs.append(EqualDOF(master_node, slave_node, ux=ux, uy=uy, rot=rot))
        return self

    def point_load(self, node: int, fx: float = 0.0, fy: float = 0.0, m: float = 0.0,
                   F: float = None, angle_deg: float = None):
        """節點集中力/力矩。fx/fy本來就是全域座標分量, 沒有local/global
        的區別問題; 如果比較習慣用「大小+角度」描述, 可以改用F+angle_deg
        (標準數學慣例, 0度=+x, 逆時針為正), 兩者算出的fx/fy會直接相加
        (通常只用其中一種, 不用同時混用)。"""
        if angle_deg is not None:
            if F is None:
                raise ValueError("指定angle_deg時必須同時指定F(力的大小)")
            ang = math.radians(angle_deg)
            fx = fx + F * math.cos(ang)
            fy = fy + F * math.sin(ang)
        self.point_loads.append(PointLoad(node, fx, fy, m))
        return self

    def distributed_load(self, member: int, w: float, w_end: float = None,
                          x_start: float = None, x_end: float = None,
                          direction: str = 'local', angle_deg: float = None):
        """均佈/線性變化載重。預設(x_start=x_end=None)整根桿件都有;
        指定x_start/x_end可以只加在桿件的局部一段(局部座標, 0<=x_start<=x_end<=L)。
        direction='local'(預設): w是局部+y方向分量。
        direction='global_y': w代表沿全域垂直方向、大小以沿桿件長度量測的
        載重(屋頂重力/雪載重的標準表示方式), 支援線性變化(w!=w_end)但
        只支援整根桿件。
        direction='global': 全域任意角度均佈載重(用angle_deg指定角度,
        度, 0=+x方向逆時針為正), 支援局部段+線性變化的任意組合。"""
        self.distributed_loads.append(
            DistributedLoad(member, w, w_end, x_start, x_end, direction, angle_deg))
        return self

    def distributed_moment(self, member: int, m: float):
        """桿件全長均佈的分布彎矩(kN·m/m, 逆時針為正, 跟
        member_point_load()的m參數同一個符號慣例)。目前只支援整根
        桿件的常數m, 不支援線性變化或局部段。

        跟distributed_load()(分布力, kN/m)不一樣: 這是分布"力矩",
        不是分布"力", 直接對桿件的內部彎矩梯度貢獻, 不像分布力那樣
        會先產生剪力再累積成彎矩——見frame2d.elements.
        fixed_end_forces_distributed_moment()的推導說明。

        實務上比較少見, 常見的動機情境: (1) 桿件截面深度方向有溫度
        梯度(例如日曬面跟背陰面溫差), 等效成沿桿長的分布彎矩載重
        來模擬溫度應力效應; (2) 一長串緊鄰、間距很密的偏心軸力/剪力
        釘, 各自產生一個小力矩, 密到可以近似成連續分布; (3) 風/水壓
        沿桿長分布但壓力中心偏離桿軸(例如壓力中心隨深度變化), 也可能
        用分布力+分布彎矩的組合來表示。不是結構分析裡最常見的載重
        類型, 大部分實務案例還是用分布力、集中力矩處理就夠。"""
        self.distributed_moments.append(DistributedMoment(member, m))
        return self

    def member_point_load(self, member: int, a: float, fx: float = 0.0, fy: float = 0.0,
                           m: float = 0.0, direction: str = 'local',
                           F: float = None, angle_deg: float = None):
        """桿件內部任意位置(距node_i為a)加集中力/集中力矩, 不新增節點。
        跟point_load(node,...)的差別: 這裡a是桿件"局部"座標(沿桿軸距離node_i
        多遠), 不是節點id; 而且direction='local'(預設)時fx/fy是該桿件
        自己的局部座標分量, 不是全域座標。
        direction='global': 改用F(大小)+angle_deg(全域角度, 標準數學慣例)
        指定, 依桿件角度自動分解成局部分量(跟distributed_load的
        direction='global'同一套邏輯), 這個模式下fx/fy不能用。"""
        self.member_point_loads.append(
            MemberPointLoad(member, a, fx, fy, m, direction, F, angle_deg))
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
