"""
Pydantic 輸入/輸出模型 -- 把 frame2d 的 dataclass 包一層,做 JSON API 用。
不改動 frame2d 核心程式碼,這層只負責「JSON <-> Frame2D / SolveResult」的轉換,
欄位名稱/預設值刻意跟 frame2d/model.py 的 dataclass 保持一致,方便對照。
"""
from typing import Optional, List, Literal
from pydantic import BaseModel, Field


class NodeIn(BaseModel):
    id: int
    x: float
    y: float


class SectionIn(BaseModel):
    name: str
    E: float
    I: float
    A: float = 1e8


class MemberIn(BaseModel):
    id: int
    node_i: int
    node_j: int
    section: str
    member_type: Literal['frame', 'truss', 'cable'] = 'frame'
    release_i: bool = False
    release_j: bool = False
    Mp_i: Optional[float] = None
    Mp_j: Optional[float] = None
    R_post_yield_i: Optional[float] = None
    R_post_yield_j: Optional[float] = None
    """塑鉸容量(選用, 見frame2d/hinge.py的HingeState), 只有'frame'元素
    有意義。預設None時完全不影響任何行為, 跟原本一樣是純彈性桿件——
    只有做pushover(analysis_type='pushover')才會用到, 線性/P-Delta分析
    完全忽略這幾個欄位(不管有沒有設定)。"""


class SupportIn(BaseModel):
    node: int
    ux: Optional[float] = None
    uy: Optional[float] = None
    rot: Optional[float] = None


class PointLoadIn(BaseModel):
    node: int
    fx: float = 0.0
    fy: float = 0.0
    m: float = 0.0


class DistributedLoadIn(BaseModel):
    member: int
    w_start: float
    w_end: Optional[float] = None
    x_start: Optional[float] = None
    x_end: Optional[float] = None
    direction: Literal['local', 'global_y', 'global'] = 'local'
    angle_deg: Optional[float] = None


class MemberPointLoadIn(BaseModel):
    member: int
    a: float
    fx: float = 0.0
    fy: float = 0.0
    m: float = 0.0
    direction: Literal['local', 'global'] = 'local'
    F: Optional[float] = None
    angle_deg: Optional[float] = None


class FrameIn(BaseModel):
    """完整結構模型的 JSON 描述,對應 frame2d.Frame2D 的所有輸入資料。"""
    nodes: List[NodeIn]
    sections: List[SectionIn]
    members: List[MemberIn]
    supports: List[SupportIn] = Field(default_factory=list)
    point_loads: List[PointLoadIn] = Field(default_factory=list)
    distributed_loads: List[DistributedLoadIn] = Field(default_factory=list)
    member_point_loads: List[MemberPointLoadIn] = Field(default_factory=list)
    units: Optional[dict] = None
    """匯出PDF/報告時, 前端目前「顯示設定」裡選的單位(E/I/A/disp/
    force/moment, 例如{"E":"GPa","force":"kN",...})——/export/pdf、
    /preview/fbd、/export/fbd_images這三個端點會用到, 讓報告顯示
    的數字/單位標籤跟使用者畫面上看到的一致, 不用另外去猜或換算;
    其他端點忽略這個欄位, 因為/solve本身進出都是SI, 跟顯示單位
    無關。"""
    member_ids: Optional[List[int]] = None
    """指定桿件(查詢畫面選一個/多個/全選)——/export/pdf、
    /preview/fbd、/export/fbd_images這三個端點會用到: 每根桿件
    附一張自由體圖(含旁邊的結構縮圖, 驗證Fx/Fy/M平衡)。其他端點
    忽略。"""
    fbd_only: Optional[bool] = False
    """搭配member_ids用: True時/export/pdf只附自由體圖(含縮圖),
    跳過每根桿件自己的N/V/M/變形圖那一頁, 讓報告更精簡。"""
    analysis_type: Literal['linear', 'pdelta', 'pushover'] = 'linear'
    """'linear'(預設, 完全等同舊行為)、'pdelta'(疊代更新軸力的線性化
    P-Delta, 見frame2d.dofmanager.solve_pdelta)、或'pushover'(遞增側推,
    見frame2d.pushover.run_pushover)。舊的呼叫端不帶這個欄位時預設
    'linear', /solve回傳格式完全不變, 不受影響。"""
    pushover_control_node: Optional[int] = None
    """analysis_type='pushover'時必填: 側推控制點的節點id(位移控制,
    這個節點會被強制推到pushover_target)。目前只支援單一控制點,
    多點同向側推(例如兩根柱頂一起推)是未來擴充。"""
    pushover_direction: Literal['x', 'y'] = 'x'
    """側推方向對應的自由度, 絕大多數情況是'x'(水平側推)。"""
    pushover_target: Optional[float] = None
    """analysis_type='pushover'時必填: 目標側推位移(m)。"""
    pushover_step: Optional[float] = None
    """analysis_type='pushover'時必填: 名目步長(m), 沒有降伏事件發生時
    每步走多少; 太大會讓降伏點定位不夠精確(event-to-event只在單一步
    內找降伏比例, 一步內如果有兩個以上獨立的降伏會漏掉後面那個),
    太小則計算變慢——一般抓target的1/50到1/200之間。"""
    pushover_use_pdelta: bool = False
    """側推過程是否同時考慮P-Delta(用目前累積軸力組幾何剛度)。"""
    pushover_mechanism_ratio_limit: float = 1e-8
    """機構偵測的縮聚剛度特徵值比例門檻, 見
    frame2d.pushover.check_mechanism()。一般不需要調整。"""


class NodeResultOut(BaseModel):
    node: int
    ux: float
    uy: float
    rot: float
    Rx: float
    Ry: float
    M: float


class MemberResultOut(BaseModel):
    member_id: int
    L: float
    angle_deg: float
    N1: float
    V1: float
    M1: float
    N2: float
    V2: float
    M2: float
    slack: bool = False


class SolveOut(BaseModel):
    nodes: List[NodeResultOut]
    members: List[MemberResultOut]
