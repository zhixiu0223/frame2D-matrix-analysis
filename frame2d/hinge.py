"""
塑性鉸狀態機 + 含鉸樑元素勁度矩陣。

移植自portal-frame-pushover-scratch(該repo已用OpenSeesPy的zeroLength轉角
彈簧+Steel01雙線性材料逐點驗證過, 單層框架全曲線誤差0.55%、雙層框架
0.02%)。這裡只搬"塑鉸狀態機+勁度公式"這一塊, 遞增側推的載重控制迴圈、
事件到事件子步進、容量曲線輸出還沒開始做, 屬於下一階段。

物理模型: 樑本身彈性(EI, L), 兩端各串接一個零長度旋轉彈簧(剛度R1, R2)。
全域自由度看到的是彈簧外側的轉角, 樑內部真正的端點轉角因彈簧變形而跟
外部轉角不同, 用靜力凝聚消去內部轉角自由度, 得到對外部自由度
(v1, theta1, v2, theta2)的4x4矩陣:

  R -> 很大: 彈簧不變形, 退化回標準剛接樑(彈性狀態)
  R -> 很小: 彈簧接近不能傳彎矩, 退化回該端鉸接(降伏硬化剛度很小時的極限)

封閉式解已用sympy從樑ODE跟彈簧力矩-轉角關係聯立推導, 通過三項極限自我
檢核(R->無限大退化回標準剛接矩陣、R->0退化回完全機構零矩陣、單端鉸接
退化回教科書3EI/L公式)——見tests/test_hinge_condensation.py。這裡直接用
推導出的封閉式公式, 不在runtime呼叫sympy(避免增加執行期相依與啟動成本,
跟elements.py的local_geometric_stiffness()同一套做法: 推導過程留在測試,
正式程式碼用推導完的封閉式公式)。
"""
import numpy as np

# 數值上代表「剛接」的彈簧剛度: 要遠大於典型EI/L量級, 用相對值而非絕對
# 常數, 避免不同單位制/尺度下數值病態。跟portal-frame-pushover-scratch
# 用同一個值, 不是重新調參數。
RIGID_FACTOR = 1.0e8


class HingeState:
    """單一元素兩端的鉸狀態。

    Mp1, Mp2: 兩端塑性彎矩容量(需從斷面性質算出, 不是調參數湊的)
    R_post_yield_1/2: 兩端降伏後的硬化旋轉剛度(跟EI同一套單位制,
        對應M-θ曲線降伏後的斜率, 不是alpha折減係數)

    theta_IO/theta_LS/theta_CP(可選, 各是長度2的tuple/list, 對應兩端):
        ASCE41/FEMA356表10-7一類驗收基準表格查出來的塑性轉角限值(rad)。
        這是「驗收準則」, 不是「材料力學性質」, 刻意跟Mp/R_post_yield分開:
        不給的話(預設None)完全不影響任何求解或勁度計算, 只有呼叫
        performance_level()時會一律回傳None, 表示「沒有分類依據」。
    """

    def __init__(self, Mp1, Mp2, R_post_yield_1, R_post_yield_2,
                 theta_IO=(None, None), theta_LS=(None, None), theta_CP=(None, None)):
        self.Mp = [Mp1, Mp2]
        self.R_post_yield = [R_post_yield_1, R_post_yield_2]
        self.yielded = [False, False]
        self.theta_p = [0.0, 0.0]   # 累積塑性轉角(rad); 下一階段(遞增側推
                                    # 迴圈)才會實際更新這個值, 這裡先開欄位
        self.theta_IO = list(theta_IO)
        self.theta_LS = list(theta_LS)
        self.theta_CP = list(theta_CP)

    def current_R(self, EI, L, end_idx):
        """回傳該端目前的彈簧剛度: 未降伏用剛接數值, 降伏用硬化剛度。"""
        if self.yielded[end_idx]:
            return self.R_post_yield[end_idx]
        return RIGID_FACTOR * EI / L

    def check_yield(self, M1, M2):
        """用回算出的端點彎矩檢查是否降伏。回傳這一步是否有新的鉸形成。
        只允許彈性->降伏的單向轉換(目前範疇不含卸載/循環載重)。"""
        newly_yielded = False
        if not self.yielded[0] and abs(M1) >= self.Mp[0]:
            self.yielded[0] = True
            newly_yielded = True
        if not self.yielded[1] and abs(M2) >= self.Mp[1]:
            self.yielded[1] = True
            newly_yielded = True
        return newly_yielded

    def performance_level(self, end_idx):
        """回傳這一端目前的性能等級分類: 'elastic' / 'IO' / 'LS' / 'CP' /
        'exceeds CP' / None(沒設定門檻值時)。純粹用累積塑性轉角theta_p
        跟建構時給的theta_IO/LS/CP門檻比較做事後分類貼標籤, 不影響任何
        求解邏輯, 也不會反過來改變current_R()的行為。"""
        if not self.yielded[end_idx]:
            return 'elastic'
        tp = self.theta_p[end_idx]
        io, ls, cp = self.theta_IO[end_idx], self.theta_LS[end_idx], self.theta_CP[end_idx]
        if io is None and ls is None and cp is None:
            return None
        if cp is not None and tp >= cp:
            return 'exceeds CP'
        if ls is not None and tp >= ls:
            return 'CP'
        if io is not None and tp >= io:
            return 'LS'
        return 'IO'


def hinge_bending_stiffness(E, I, L, hinge_state):
    """回傳含鉸效應的4x4撓曲勁度矩陣(對v1, theta1, v2, theta2), 用靜力
    凝聚封閉式公式(見本檔案開頭說明), 依目前鉸狀態代入對應的R1, R2。"""
    EI = E * I
    R1 = hinge_state.current_R(EI, L, 0)
    R2 = hinge_state.current_R(EI, L, 1)
    denom = 12 * EI**2 + 4 * EI * L * R1 + 4 * EI * L * R2 + L**2 * R1 * R2
    k11 = 12 * EI * (EI * R1 + EI * R2 + L * R1 * R2) / (L**2 * denom)
    k12 = 6 * EI * R1 * (2 * EI + L * R2) / (L * denom)
    k14 = 6 * EI * R2 * (2 * EI + L * R1) / (L * denom)
    k22 = 4 * EI * R1 * (3 * EI + L * R2) / denom
    k24 = 2 * EI * L * R1 * R2 / denom
    k44 = 4 * EI * R2 * (3 * EI + L * R1) / denom
    return np.array([
        [k11,  k12, -k11,  k14],
        [k12,  k22, -k12,  k24],
        [-k11, -k12,  k11, -k14],
        [k14,  k24, -k14,  k44],
    ])


def beam_internal_rotation(E, I, L, hinge_state, v1, theta1, v2, theta2):
    """回傳樑本身內部端點的真實轉角(phi1, phi2), 不是外部節點自由度
    (theta1, theta2)本身——這兩者只有在彈簧沒有變形(彈性、剛接近似)時
    才會幾乎相等; 鉸降伏後, 彈簧會有顯著的相對轉動, 樑內部端點轉角跟
    外部節點轉角會明顯不同(尤其是鉸剛好在完全固定支承上這種情況:
    外部節點轉角theta1被邊界條件釘死在0, 但樑內部端點phi1可以因為
    彈簧軟化而顯著轉動, 這正是"塑性鉸轉動"這件事在數學上的體現)。

    封閉式解跟hinge_bending_stiffness()用的是同一組聯立方程式(樑ODE
    的M1,M2 = 彈簧力矩-轉角關係R1*(theta1-phi1), R2*(theta2-phi2)),
    只是這裡解出的是phi1,phi2本身而不是消去它們之後的縮聚矩陣。用途:
    pushover.py拿(theta1-phi1), (theta2-phi2)當作該端"目前的塑性轉角"
    (theta_p), 給後續視覺化(圈圈大小)跟IO/LS/CP分類用。"""
    EI = E * I
    R1 = hinge_state.current_R(EI, L, 0)
    R2 = hinge_state.current_R(EI, L, 1)
    denom = L * (12 * EI**2 + 4 * EI * L * R1 + 4 * EI * L * R2 + L**2 * R1 * R2)
    phi1 = (-12 * EI**2 * v1 + 12 * EI**2 * v2 + 4 * EI * L**2 * R1 * theta1
            - 2 * EI * L**2 * R2 * theta2 - 6 * EI * L * R2 * v1 + 6 * EI * L * R2 * v2
            + L**3 * R1 * R2 * theta1) / denom
    phi2 = (-12 * EI**2 * v1 + 12 * EI**2 * v2 - 2 * EI * L**2 * R1 * theta1
            + 4 * EI * L**2 * R2 * theta2 - 6 * EI * L * R1 * v1 + 6 * EI * L * R1 * v2
            + L**3 * R1 * R2 * theta2) / denom
    return phi1, phi2


def full_6x6_with_hinge(E, A, I, L, hinge_state, P=0.0):
    """組出完整6x6局部剛度矩陣(軸向不受鉸影響 + 撓曲用含鉸凝聚公式 +
    選用的P-Delta幾何剛度), 自由度順序跟elements.py一致:
    [u1, v1, theta1, u2, v2, theta2]。

    P: 局部座標下的元素軸力, 跟elements.local_geometric_stiffness()同一
    慣例(拉力為正)。預設0.0時完全不含幾何剛度, 純塑鉸效應。這裡沿用
    dofmanager.py已經確立的簡化(見_solve_once_dofmanager docstring):
    Kg直接疊加, 不管端點鉸狀態, 跟已驗證過的portal-frame-pushover-scratch
    採同一種簡化。"""
    from .elements import local_geometric_stiffness

    k = np.zeros((6, 6))
    EA_L = E * A / L
    k[0, 0] = EA_L
    k[0, 3] = -EA_L
    k[3, 0] = -EA_L
    k[3, 3] = EA_L

    kb = hinge_bending_stiffness(E, I, L, hinge_state)
    idx = [1, 2, 4, 5]
    for i, ii in enumerate(idx):
        for j, jj in enumerate(idx):
            k[ii, jj] = kb[i, j]

    if P != 0.0:
        kg = local_geometric_stiffness(P, L)
        k = k + kg

    return k
