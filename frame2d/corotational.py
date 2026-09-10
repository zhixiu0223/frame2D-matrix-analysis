"""
真正的大轉角co-rotational(共旋)公式, 供frame2d.newton這個真正的
Newton-Raphson平衡疊代求解器使用。

核心概念(標準co-rotational樑元素理論, 不是這個repo自己發明的):
把節點的"總變形"拆成兩部分——(1)桿件整根一起平移+轉動的"剛體運動"
(不產生任何內力, 不管轉了多少度), (2)桿件相對於自己(相對弦線)的
"純變形"(拉伸量e、兩端相對弦線的轉角theta1_def/theta2_def)。內力
只由(2)決定, 跟(1)完全無關——這保證"桿件被剛體轉動再多度, 只要它
自己沒有真的拉伸/彎曲, 內力就一定是0", 這正是run_pushover()裡
geometry_update選項欠缺、之前在對話紀錄裡被實際案例找出來的那個限制
(拿一根桿件的原始轉角自由度直接套用小角度公式, 沒有先扣掉剛體轉動
了多少)。

跟frame2d.hinge的關係: 塑鉸(雙折線M-theta)的物理量本身沒有變(還是
同一個Mp/R_post_yield模型), 只是這裡餵給它的"theta"改成theta_def
(相對弦線的變形轉角), 不是節點的總轉角——這才是正確的做法, 因為
塑鉸真正感受到的是"這根桿件自己被凹了多少", 不是"這根桿件在全域
座標裡轉了多少度"(後者含有大量跟塑性完全無關的剛體轉動)。

已知限制(誠實記錄, 不誇大):
- 這裡的"樑本身(chord-relative)彎矩-轉角關係"沿用hinge_bending_
  stiffness()既有的封閉式公式, 只是用current L(目前變形後的桿長)
  代入, 不是完整的有限應變(finite strain)樑理論——對絕大多數結構
  分析場景(轉角到幾十度量級)這個近似足夠, 但不是連續介質力學意義上
  完全嚴謹的大應變公式。
- 切線剛度矩陣用有限差分(數值微分)算, 不是解析導數——這是刻意的
  設計選擇, 不是偷懶: 解析導數手推容易出正負號/係數錯誤(這正是這個
  co-rotational模組存在的理由: 上一版的geometry_update就是因為
  沒有處理好這類幾何細節才出問題), 有限差分保證跟"內力函式"本身
  完全一致(不會有"內力對, 但切線推錯"這種難以察覺的問題), 代價是
  疊代次數可能比解析切線多一點、每次疊代稍微慢一點——用Newton法時
  重要的是"殘餘力(內力函式)本身要對", 切線只是用來找方向, 不完美
  的切線不影響收斂到的答案是否正確, 只影響收斂速度。
"""
import numpy as np

from .hinge import hinge_bending_stiffness


def corotational_kinematics(xi, yi, ui, vi, thetai, xj, yj, uj, vj, thetaj):
    """給定桿件原始座標(xi,yi)-(xj,yj)跟目前的節點絕對位移/轉角(ui,vi,
    thetai,uj,vj,thetaj), 回傳:
        L0: 原始桿長
        L: 目前(變形後)桿長
        beta: 目前弦線(i->j)方位角(rad)
        e: 軸向"自然應變"量測, e = L - L0
        theta1_def, theta2_def: 兩端相對弦線的"變形轉角"(剛體轉動已經
            扣掉), 這才是塑鉸/彎曲勁度該吃的"轉角", 不是thetai/thetaj
            本身
        B: 3x6矩陣, 把(δe, δtheta1_def, δtheta2_def)這組"自然變形速率"
            跟6個全域自由度的虛位移(δui,δvi,δthetai,δuj,δvj,δthetaj)
            關聯起來(虛功一致的B矩陣, 標準co-rotational推導方式,
            推導過程見tests/test_corotational.py的驗證, 不是憑空寫的)
    """
    dx0, dy0 = xj - xi, yj - yi
    L0 = float(np.hypot(dx0, dy0))
    beta0 = float(np.arctan2(dy0, dx0))

    xci, yci = xi + ui, yi + vi
    xcj, ycj = xj + uj, yj + vj
    dx, dy = xcj - xci, ycj - yci
    L = float(np.hypot(dx, dy))
    beta = float(np.arctan2(dy, dx))

    # 剛體轉動量: 目前弦線方位角 - 原始弦線方位角。用atan2取出來的角度
    # 本來就落在(-pi,pi], 兩個角度相減理論上可能跨到(-2pi,2pi]這個
    # 範圍——但這裡刻意"不做"環繞修正(wrap成(-pi,pi]), 因為如果桿件
    # 真的被推轉了超過180度, 這代表模型已經深入某種物理上不合理的
    # 崩塌狀態, 讓psi維持"跨圈"的原始值, 後面的max_rotation檢查才能
    # 誠實反映出這個異常, 不要在這裡悄悄把它折回合理範圍、製造"看起來
    # 沒事"的假象。
    psi = beta - beta0

    e = L - L0
    theta1_def = thetai - psi
    theta2_def = thetaj - psi

    c, s = np.cos(beta), np.sin(beta)
    B = np.array([
        [-c, -s, 0.0,  c,  s, 0.0],
        [-s / L, c / L, 1.0, s / L, -c / L, 0.0],
        [-s / L, c / L, 0.0, s / L, -c / L, 1.0],
    ])
    return L0, L, beta, e, theta1_def, theta2_def, B


def corotational_local_forces(E, A, I, L0, L, e, theta1_def, theta2_def, hinge_state=None,
                               M_ref=None, theta_def_ref=None):
    """給定"自然變形量"(e, theta1_def, theta2_def), 回傳對應的局部內力
    (N, M1, M2)——這一步完全不含任何幾何/剛體轉動資訊(那些都已經在
    corotational_kinematics()裡處理掉了), 這裡純粹是"這根桿件被拉伸
    e、兩端相對弦線轉了theta1_def/theta2_def, 材料要出多少力/彎矩"。

    hinge_state: None時用標準彈性樑公式(4EI/L, 2EI/L)——彈性關係全程
    不變, 直接用"總量"公式即可, M_ref/theta_def_ref不會被用到。

    有給hinge_state時, 彎矩改成從(M_ref, theta_def_ref)這個參考狀態
    "增量"算出來, 不是直接把theta1_def,theta2_def的總量套進當下的
    R1,R2公式——這是為了修正一個真實發生過的bug(見對話紀錄): 塑鉸
    降伏瞬間R會從"近似剛接"(~1e8*EI/L量級)驟降到R_post_yield(通常
    只有EI/L的百分之幾), 如果直接拿同一個theta_def總量套用新的(小
    很多的)R重新算一次M, 算出來的M會遠低於Mp、瞬間"消失"大半, 這在
    物理上不合理(雙折線硬化模型, 降伏後M只會從Mp緩慢往上加, 不會
    倒退)。用增量(dM = K2(目前R) @ (theta_def - theta_def_ref), 再
    加回M_ref)才能保持彎矩在降伏瞬間連續, 這正是run_pushover()裡
    event-to-event那套"cum_forces累加df"在做的事, 這裡把同樣的精神
    帶進co-rotational的Newton版本。

    M_ref/theta_def_ref: 分別是這根桿件上一次"成功收斂、已經被接受"
    的那一步的(M1,M2)/(theta1_def,theta2_def), 不給(None)時視為
    (0,0)(也就是從沒有任何內力/變形的狀態開始算, 對第一步或者一直
    保持彈性的情況, 這樣算出來的結果會等同用總量公式直接算, 完全
    一致)。

    回傳(N, M1, M2)。
    """
    N = E * A / L0 * e
    if hinge_state is None:
        kb = E * I / L
        M1 = 4 * kb * theta1_def + 2 * kb * theta2_def
        M2 = 2 * kb * theta1_def + 4 * kb * theta2_def
    else:
        K4 = hinge_bending_stiffness(E, I, L, hinge_state)
        K2 = K4[np.ix_([1, 3], [1, 3])]
        M1_ref, M2_ref = (0.0, 0.0) if M_ref is None else M_ref
        th1_ref, th2_ref = (0.0, 0.0) if theta_def_ref is None else theta_def_ref
        dtheta = np.array([theta1_def - th1_ref, theta2_def - th2_ref])
        dM1, dM2 = K2 @ dtheta
        M1, M2 = M1_ref + dM1, M2_ref + dM2
    return float(N), float(M1), float(M2)


def corotational_global_force(xi, yi, ui, vi, thetai, xj, yj, uj, vj, thetaj,
                               E, A, I, hinge_state=None, M_ref=None, theta_def_ref=None):
    """組合kinematics+local_forces, 直接回傳這根桿件目前狀態下的6維
    全域節點內力向量(f = B^T @ [N,M1,M2]), 順序跟其餘frame2d程式碼
    一致: [Fxi,Fyi,Mi,Fxj,Fyj,Mj]。這是Newton-Raphson疊代的"內力"
    函式本身, 必須是精確的(殘餘力的正確性完全靠這個函式)。

    M_ref/theta_def_ref: 見corotational_local_forces()說明, 原封不動
    往下傳, 塑鉸的路徑相依性(降伏後彎矩要從Mp連續往上加, 不能每次
    重算總量)完全靠這兩個參數維持。"""
    L0, L, beta, e, theta1_def, theta2_def, B = corotational_kinematics(
        xi, yi, ui, vi, thetai, xj, yj, uj, vj, thetaj)
    N, M1, M2 = corotational_local_forces(
        E, A, I, L0, L, e, theta1_def, theta2_def, hinge_state, M_ref, theta_def_ref)
    f_local = np.array([N, M1, M2])
    return B.T @ f_local


def corotational_tangent_fd(xi, yi, ui, vi, thetai, xj, yj, uj, vj, thetaj,
                             E, A, I, hinge_state=None, M_ref=None, theta_def_ref=None, eps=1e-7):
    """有限差分算切線剛度矩陣(6x6): 對6個全域自由度各自微小擾動, 看
    corotational_global_force()怎麼變化, 组成Jacobian。見本檔案開頭
    docstring說明為什麼刻意選有限差分而不是手推解析式。

    M_ref/theta_def_ref在整個有限差分擾動過程中保持固定(呼叫端"這一次
    Newton疊代"的參考狀態不會因為擾動了u而跟著變, 這是正確的做法:
    參考狀態只在"一整步"被接受收斂時才更新一次, 不是每次疊代都變)。

    eps用"絕對"擾動量(不是相對), 呼叫端要注意如果模型單位制差異很大
    (例如公尺 vs 公厘), 可能要調整這個值——但對這個repo一貫使用的SI
    單位(公尺, 弧度), 1e-7是合理的量級(比典型位移小7個數量級, 比
    浮點數精度極限(~1e-16乘上典型量級)大很多, 不會被捨入誤差淹沒)。
    """
    u0 = np.array([ui, vi, thetai, uj, vj, thetaj], dtype=float)
    f0 = corotational_global_force(xi, yi, *u0[0:3], xj, yj, *u0[3:6], E, A, I,
                                    hinge_state, M_ref, theta_def_ref)
    K = np.zeros((6, 6))
    for k in range(6):
        u_pert = u0.copy()
        u_pert[k] += eps
        f_pert = corotational_global_force(xi, yi, *u_pert[0:3], xj, yj, *u_pert[3:6], E, A, I,
                                            hinge_state, M_ref, theta_def_ref)
        K[:, k] = (f_pert - f0) / eps
    return K
