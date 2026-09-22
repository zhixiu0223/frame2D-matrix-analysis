"""
Pydantic 輸入/輸出模型 -- 把 frame2d 的 dataclass 包一層,做 JSON API 用。
不改動 frame2d 核心程式碼,這層只負責「JSON <-> Frame2D / SolveResult」的轉換,
欄位名稱/預設值刻意跟 frame2d/model.py 的 dataclass 保持一致,方便對照。
"""
from typing import Optional, List, Literal, Dict
from pydantic import BaseModel, Field


class NodeIn(BaseModel):
    id: int
    x: float
    y: float
    mx: Optional[float] = None    # 節點質量(動力分析用, 後端固定SI: kg), x/y方向平動質量
    my: Optional[float] = None
    Iz: Optional[float] = None    # 繞z軸質量轉動慣量(kg·m²)


class SectionIn(BaseModel):
    name: str
    E: float
    I: float
    A: float = 1e8
    alpha: Optional[float] = None   # 熱膨脹係數(1/°C或1/K), 只有thermal_loads
                                     # 要用到才需要填
    depth: Optional[float] = None   # 截面深度(m), 只有溫度梯度(彎曲熱效應)
                                     # 才需要填
    rho: Optional[float] = None     # 質量密度(kg/m³, 後端固定SI), 只有動力分析
                                     # (質量矩陣)才需要填, 靜力分析完全忽略


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


class EqualDofIn(BaseModel):
    """跟frame2d.model.EqualDOF對應——見Frame2D.equal_dof()的說明,
    讓slave_node指定的自由度強制等於master_node的對應自由度(用高
    勁度彈簧懲罰法實現, 不是完全精確為0的束制, 但精度足夠工程使用)。

    支援範圍: analysis_type='linear'/'pdelta', 以及全部四種pushover
    求解器(event_to_event/converged/newton/corotational_oneshot)
    都支援。event_to_event/converged底層走frame2d.pushover.
    _assemble_stiffness_with_hinges(), 跟linear/pdelta共用同一個
    dofmanager._apply_equal_dof()函式; newton/corotational_oneshot
    底層走frame2d.newton.py完全獨立的co-rotational組裝邏輯, 用另一個
    專屬的_apply_equal_dof_newton()(除了修正切線剛度矩陣, 還要同步
    修正殘餘力向量, 因為Newton疊代檢查的是殘餘力有沒有收斂到0, 光改
    剛度矩陣殘餘力不會反映懲罰彈簧的存在)——懲罰彈簧的勁度倍率也
    刻意比linear/pdelta用的小一個數量級(1e5而不是1e6), 因為Newton
    疊代對矩陣病態比一次性線性解敏感得多, 倍率太大會讓疊代完全無法
    收斂(不是變慢, 是卡住), 這是實際測試發現、修正過的問題。"""
    master_node: int
    slave_node: int
    ux: bool = False
    uy: bool = False
    rot: bool = False


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


class DistributedMomentIn(BaseModel):
    """跟frame2d.model.DistributedMoment對應——見Frame2D.
    distributed_moment()的說明, 桿件的分布彎矩(kN·m/m, 逆時針為正,
    跟member_point_load()的m參數同一個符號慣例)。預設整根桿件都有
    固定值m; 可以指定x_start/x_end只加在局部一段(跟
    DistributedLoadIn同一套慣例), 也可以指定m_end讓強度線性變化
    (梯形分布彎矩)。彎矩沒有"方向"這個概念, 沒有對應
    DistributedLoadIn的direction/angle_deg欄位。"""
    member: int
    m: float
    m_end: Optional[float] = None
    x_start: Optional[float] = None
    x_end: Optional[float] = None


class ThermalLoadIn(BaseModel):
    """跟frame2d.model.ThermalLoad對應——見Frame2D.thermal_load()的
    說明, 桿件的溫度效應載重: delta_T(均勻溫度變化, 軸向熱效應, 需要
    斷面設定alpha)、delta_T_top/delta_T_bottom(截面頂/底溫度變化,
    彎曲熱效應, 兩者不同時才有效果, 需要斷面同時設定alpha跟depth)。
    這兩個分量各自獨立, 可以只給其中一個。"""
    member: int
    delta_T: float = 0.0
    delta_T_top: Optional[float] = None
    delta_T_bottom: Optional[float] = None


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
    equal_dofs: List[EqualDofIn] = Field(default_factory=list)
    point_loads: List[PointLoadIn] = Field(default_factory=list)
    distributed_loads: List[DistributedLoadIn] = Field(default_factory=list)
    distributed_moments: List[DistributedMomentIn] = Field(default_factory=list)
    thermal_loads: List[ThermalLoadIn] = Field(default_factory=list)
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
    pushover_step_member_forces: Optional[Dict[str, List[float]]] = None
    """只有/pushover_step_diagrams這個端點會用到: 從Pushover回放某一步
    的快照(history_snapshots[i].member_forces)直接原封不動傳回來,
    key是桿件id(字串), value是[Fx1,Fy1,M1,Fx2,Fy2,M2]。"""
    pushover_step_u_full: Optional[List[float]] = None
    """只有/pushover_step_diagrams這個端點會用到: 從Pushover回放某一步
    的快照(history_snapshots[i].u_full)直接原封不動傳回來, 長度是
    3*節點數(+release造成的額外自由度), 跟frame.dofs_of()的編號對應。"""
    fbd_only: Optional[bool] = False
    """搭配member_ids用: True時/export/pdf只附自由體圖(含縮圖),
    跳過每根桿件自己的N/V/M/變形圖那一頁, 讓報告更精簡。"""
    cyclic_control_nodes: Optional[List[int]] = None
    """只有/cyclic端點用: 反覆載重的控制節點(第一個是參考點, 遲滯迴圈的位移軸就是它的位移)。"""
    cyclic_weights: Optional[List[float]] = None
    """各控制節點的比例, None=全部1.0(跟pushover_weights同義)。"""
    cyclic_direction: Literal['x', 'y'] = 'x'
    cyclic_amplitudes: Optional[List[float]] = None
    """位移幅值序列(SI, m), 例如[0.015, 0.03, 0.045, 0.06]: 每個幅值走 cyclic_n_cycles 個(+a,-a)循環。"""
    cyclic_n_cycles: int = 2
    cyclic_step: Optional[float] = None
    """位移步長上限(SI, m); 降伏事件會自動在事件點切開, 不用為了抓事件調小。"""

    rsa_direction: Literal['x', 'y'] = 'x'
    rsa_damping: float = 0.05
    rsa_combine: Literal['SRSS', 'CQC'] = 'SRSS'
    rsa_n_modes: Optional[int] = None
    rsa_mass_kind: Literal['lumped', 'consistent'] = 'lumped'
    rsa_spectrum_type: Literal['code', 'custom'] = 'code'
    rsa_code_sds: Optional[float] = None
    """規範反應譜用: 無單位係數(g的倍數)。"""
    rsa_code_sd1: Optional[float] = None
    rsa_code_tl: Optional[float] = 6.0
    """規範反應譜用: 長週期轉角週期(s)。"""
    rsa_custom_points: Optional[List[List[float]]] = None
    """自訂反應譜用: [[T(s), Sa(SI加速度單位, m/s²)], ...], 至少2點。"""

    modal_n_modes: Optional[int] = None
    """只有/modal端點用: 要幾個模態(從最低頻算起), None=全部。"""
    modal_mass_kind: Literal['lumped', 'consistent'] = 'lumped'
    """只有/modal端點用: 集中質量(預設)或一致質量, 見frame2d/mass.py。"""
    analysis_type: Literal['linear', 'pdelta', 'pushover'] = 'linear'
    """'linear'(預設, 完全等同舊行為)、'pdelta'(疊代更新軸力的線性化
    P-Delta, 見frame2d.dofmanager.solve_pdelta)、或'pushover'(遞增側推,
    見frame2d.pushover.run_pushover)。舊的呼叫端不帶這個欄位時預設
    'linear', /solve回傳格式完全不變, 不受影響。"""
    pushover_control_node: Optional[int] = None
    """analysis_type='pushover'時必填(除非用下面的pushover_control_nodes
    指定多點): 側推控制點的節點id(這個節點會被強制推到pushover_target,
    或是力控制時在這個點施加力)。"""
    pushover_control_nodes: Optional[List[int]] = None
    """多點側推用(例如同時推兩三個樓層做倒三角形/均勻型態): 節點id清單。
    給了這個欄位時優先於pushover_control_node(單點), 兩者不會同時生效。
    清單第一個節點是"參考點"——容量曲線的位移軸(history_u)固定回報
    這個節點的位移, 不管它自己的權重是多少(業界慣例: 不管用哪種側推
    型態, 容量曲線永遠看屋頂/最高樓層位移, 所以第一個節點通常應該填
    屋頂節點)。"""
    pushover_weights: Optional[List[float]] = None
    """搭配pushover_control_nodes用, 跟它等長, 各節點的相對權重(位移
    控制時是相對位移比例, 力控制時是相對力比例)。不給的話全部預設1.0
    (等權重, 對應FEMA356/ASCE41兩種標準側推型態之一的"均勻型態");
    要做"倒三角形型態"就照各樓層高度比例填權重, 例如三層樓
    [0.33, 0.67, 1.0]。"""
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
    pushover_control_mode: Literal['displacement', 'force'] = 'displacement'
    """'displacement'(預設, 位移控制, 可以穿越極限承載力之後的軟化段,
    但多點時是"強制固定相對位移比例", 不是規範標準做法)或'force'
    (力控制, 直接施加已知的力, 多點時"固定力型態、位移形狀隨降伏自然
    演化"才是FEMA356/ASCE41的標準做法——但沒辦法穿越極限承載力之後的
    軟化段, 到達那個點會優雅停止(mechanism_reached=True), 見
    frame2d.pushover.run_pushover()。此時pushover_target/pushover_step
    代表的是力(N), 不是位移(m)。"""
    pushover_geometry_update: bool = False
    """False(預設, 完全等同原本行為, 全程用最初始的未變形幾何組裝勁度
    矩陣, 只靠線性化的P-Delta修正項近似大變形效應)或True(每一步用目前
    累積變形後的節點位置重新算桿件長度/角度, Updated Lagrangian的
    基本精神, 比線性化P-Delta更接近大變形時的真實行為, 但不是完整的
    co-rotational大轉角分析)。見frame2d.pushover.run_pushover()的
    geometry_update參數說明。"""
    pushover_solver: Literal['event_to_event', 'converged', 'newton', 'corotational_oneshot'] = 'event_to_event'
    """'event_to_event'(預設, 對應frame2d.pushover.run_pushover(): 塑鉸
    狀態凍結的每一段區間內, 只用該段開始時的幾何/軸力組一次勁度矩陣
    解一次, 不會檢查這個答案在真正的變形終點上是否還跟一開始用的幾何/
    軸力一致)、'converged'(對應frame2d.pushover.run_pushover_converged():
    同一段區間內反覆疊代到幾何/軸力自洽為止, 疊代不收斂時會誠實回報
    mechanism_reached=True提前停止, 不會給不可信的答案——但收斂只保證
    數值上自洽, 不保證轉角還在小角度假設的有效範圍內, 這是兩個獨立的
    問題, 見run_pushover_converged()的說明)、'newton'(對應
    frame2d.newton.run_pushover_newton(): 真正的大轉角co-rotational
    幾何+Newton-Raphson平衡疊代, 材料非線性跟幾何非線性都用嚴謹疊代
    解, 轉角再大也不會像前兩者那樣失真, 是四者裡最嚴謹也最慢的——
    目前不支援release端(有的話直接回400錯誤); 重力/桿件內部載重
    (均佈載重、桿件內部集中力)有支援, 用固定端反力公式轉成等效節點力
    (跟其餘求解器的apply_gravity()同一種"算一次、全程凍結"簡化), 已
    驗證跟apply_gravity()精確一致, 見run_pushover_newton()
    docstring的已知限制)、或'corotational_oneshot'(對應
    frame2d.newton.run_pushover_corotational_oneshot(): 用跟'newton'
    完全同一套co-rotational元素公式(精確大轉角), 但求解策略換成跟
    event_to_event同等級的"每一步只解一次, 不疊代到殘餘力收斂"——
    速度比'newton'快, 精確度介於'converged'(仍是小角度近似)跟
    'newton'(疊代到收斂)之間; 限制跟'newton'相同(不支援release端);
    這是驗證性質較重的選項, 見ANALYSIS_ARCHITECTURE.md的完整說明)。
    前兩者材料非線性(塑鉸降伏)邏輯完全相同, 差異只在幾何/軸力這一塊
    怎麼解;'newton'跟'corotational_oneshot'共用同一套co-rotational
    元素公式, 差異只在要不要疊代到殘餘力收斂——**重要**: 只有同一個
    幾何模型內部的比較(event_to_event↔converged, 或newton↔
    corotational_oneshot)才是「同一組方程式、只是解法嚴謹度不同」,
    理論上步長夠細時誤差會趨近於0;跨越'event_to_event'/'converged'
    (線性化幾何)跟'newton'/'corotational_oneshot'(精確大轉角)這兩組
    之間的比較, 是兩個不同的幾何非線性模型, 差距不會隨步長縮小而
    消失, 那個差距本身就是「小角度近似在這個變形量下還可不可信」的
    有意義資訊, 不是求解精確度的問題(2025-09對話裡驗證過: 線性化
    幾何跟co-rotational在中等位移下的差距穩定在0.07%附近, 步長切10倍
    幾乎沒有變化)。"""
    pushover_geom_tol: float = 1e-6
    """只有pushover_solver='converged'時有意義: 幾何/軸力疊代的相對
    收斂容忍度(無因次)。"""
    pushover_max_geom_iter: int = 30
    """只有pushover_solver='converged'時有意義: 每一段最多疊代幾次,
    超過視為這一段解不出來(不收斂, 提前停止)。"""
    pushover_newton_tol: float = 1e-6
    """只有pushover_solver='newton'時有意義: 殘餘力收斂容忍度(相對
    於目前內力量級的比值)。"""
    pushover_newton_max_iter: int = 30
    """只有pushover_solver='newton'時有意義: 每一步最多疊代幾次,
    超過視為這一步不收斂(提前停止側推)。"""


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
