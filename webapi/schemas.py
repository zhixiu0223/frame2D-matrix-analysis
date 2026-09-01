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
