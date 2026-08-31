"""
後處理模組: 結構圖 / 受力圖 / 軸力圖 / 剪力圖 / 彎矩圖 / 變形圖

設計原則(呼應規劃階段的討論): 不要讓繪圖程式自己重新計算力學,
所有數值都是從 SolveResult 這份「唯一的正確答案」推導出來的不同視角。

桿端內力 N(x)/V(x)/M(x) 公式推導與驗證:
  N(x) = -Fx1_local            (Fx1本身是"壓力為正"的節點力慣例, 取負號校正成
                                 標準工程慣例"拉力為正"。這個bug是做truss驗證
                                 時才抓到的: test_truss.py用一根單純受拉的桿件
                                 驗證軸力應該是正的, 結果早期版本(N=+Fx1)對一個
                                 A字撐架的壓力桿算出正值, 顯示成好像在受拉,
                                 已修正並在test_truss.py鎖住這個正負號)
  V(x) = Fy1_local + W(x)              W(x) = 0到x的累積分佈載重
  M(x) = -M1_local + Fy1_local*x + ∫[0,x] W(s) ds

  注意 M(x) 前面是 -M1_local (不是+M1_local): 這是 FEM 節點力慣例
  (M1,M2都是CCW正, 兩端同號) 跟古典彎矩圖慣例(全長連續, 上下拉力側正負一致)
  之間的標準轉換關係 M(0)=-M1, M(L)=+M2。之前只用簡支梁(M1=M2≈0)測過
  這條公式, 兩端剛好都是0, 正負號錯了也測不出來; 後來拿懸臂梁(M1=40非零)
  反推才抓到 M(L)應該等於M2, 而不是+M1+Fy1*L那樣算出來的80。
  已重新驗證: 懸臂梁 M(L)=0=M2 (原本錯誤版本算出80), 簡支梁跨中M=22.5=wL²/8 依然吻合。

  桿件內部集中力/力矩(member_point_load)的貢獻(見tests/test_member_point_load.py
  的簡支梁決定性反力反推驗證): 在a點處,
    N(x) 於 x>=a 處跳躍 -fx (fx沿局部+x, N拉力為正, 符號跟軸力慣例對齊)
    V(x) 於 x>=a 處跳躍 +fy (剛好就是fy本身, 因為V(x)本來就是"Fy1+累積量"這套邏輯)
    M(x) 於 x>=a 處額外貢獻 fy*(x-a) (V跳躍的積分, 跟連續分佈載重的累積量同一套機制)
    M(x) 於 x>=a 處另外再跳躍 -m (集中力矩造成的彎矩不連續, 不影響V;
      這個負號曾經寫反過, 詳見fixed_end_forces_point_moment()的說明)

變形內插用 Hermite cubic shape function (homogeneous解), 節點值精確,
但若桿件內部有分佈載重或集中力(member_point_load)、且只用單一元素代表
整根桿件, 內插的跨間撓度會低估實際下垂量(因為真實解在這些情況下是四次
甚至分段三次多項式, 不是單一三次)——均佈載重時誤差通常較小, 但**桿件內部
集中力的誤差可能相當明顯(實測case: L=10,P=30,a=3.5, 單一桿件內插的
最大撓度比精確節點解低估了約23%, 位置也有偏差)**。若要精確的跨間撓度數值
(例如要拿去跟SW FEA這類逐點計算的工具比對), 建議直接在該位置插入一個
真實節點、把桿件拆成兩段, 這樣該點的位移就是FEM節點解, 沒有內插誤差
(見examples/load_system_v2_demo.py的說明或直接在該點加node+point_load)。
"""
import numpy as np

from .elements import member_geometry


def _w_at(w_start, w_end, x, L):
    return w_start + (w_end - w_start) * x / L


def _cumulative_W(w_start, w_end, x, L):
    """W(x) = integral of w(s) ds from 0 to x"""
    dw = w_end - w_start
    return w_start * x + dw * x**2 / (2 * L)


def _cumulative_M_from_W(w_start, w_end, x, L):
    """integral of W(s) ds from 0 to x (雙重積分, 用於M(x))"""
    dw = w_end - w_start
    return w_start * x**2 / 2 + dw * x**3 / (6 * L)


def _partial_load_W(w_start, w_end, c, d, x):
    """局部段[c,d]均佈/線性變化載重, 對每個x位置的累積載重貢獻W(x)
    (x<=c時為0, c<x<=d時是部分積分, x>d時是整段積分後維持定值)。
    c=0,d=L時應退化成 _cumulative_W() 的結果。"""
    x = np.asarray(x, dtype=float)
    W = np.zeros_like(x)
    span = d - c
    if span <= 1e-12:
        return W
    dw = w_end - w_start
    mask2 = (x > c) & (x <= d)
    u = x[mask2] - c
    W[mask2] = w_start * u + dw * u**2 / (2 * span)
    mask3 = x > d
    W[mask3] = w_start * span + dw * span / 2.0
    return W


def _partial_load_M(w_start, w_end, c, d, x):
    """局部段[c,d]均佈/線性變化載重, 對每個x位置的M(x)貢獻
    (_partial_load_W的積分)。c=0,d=L時應退化成 _cumulative_M_from_W() 的結果。"""
    x = np.asarray(x, dtype=float)
    Mc = np.zeros_like(x)
    span = d - c
    if span <= 1e-12:
        return Mc
    dw = w_end - w_start
    mask2 = (x > c) & (x <= d)
    u = x[mask2] - c
    Mc[mask2] = w_start * u**2 / 2 + dw * u**3 / (6 * span)
    mask3 = x > d
    M_at_d = w_start * span**2 / 2 + dw * span**2 / 6
    W_at_d = w_start * span + dw * span / 2.0
    Mc[mask3] = M_at_d + W_at_d * (x[mask3] - d)
    return Mc


def member_internal_forces(frame, result, member_id, n=21):
    """回傳局部座標系下沿桿長分佈的 (x, N, V, M) 陣列。
    取樣點會自動包含每個member_point_load的確切位置、以及每個局部段
    distributed_load的起訖位置(左右各插一個極近的點), 這樣繪圖時跳躍
    不連續處才會畫成真正的垂直線, 不會被線性內插抹平成斜線。
    """
    mr = result.member_results[member_id]
    L = mr.L
    Fx1, Fy1, M1 = mr.end_forces_local[0], mr.end_forces_local[1], mr.end_forces_local[2]

    dls = [dl for dl in frame.distributed_loads if dl.member == member_id]
    ranges = []       # 局部y方向(橫向)的均佈載重範圍, 影響V(x)/M(x)
    axial_ranges = []  # 局部x方向(軸向)的均佈載重範圍, 影響N(x)
    m_obj = frame.members[member_id]
    for dl in dls:
        c = 0.0 if dl.x_start is None else dl.x_start
        d = L if dl.x_end is None else dl.x_end
        if dl.direction == 'global_y':
            # 全域垂直方向載重, 依桿件角度分解成局部x(軸向)+局部y(橫向)
            # 兩個分量(跟dofmanager.py組裝時用的同一套分解邏輯, 用T矩陣
            # 對(0,-w)向量做旋轉)。direction='global_y'目前只支援均佈+
            # 整根桿件(model.py的DistributedLoad.__post_init__已經擋掉
            # 其他組合), 所以這裡c,d必然是0,L, w_start=w_end。
            ni, nj = frame.nodes[m_obj.node_i], frame.nodes[m_obj.node_j]
            _, angle = member_geometry(ni, nj)
            c_ang, s_ang = np.cos(angle), np.sin(angle)
            # 局部 = R(-angle) @ 全域, R(-angle)= [[cosθ, sinθ],[-sinθ, cosθ]];
            # 全域載重向量固定是(0, -w)(垂直向下, 大小w)
            w_local_x = s_ang * (-dl.w_start)
            w_local_y = c_ang * (-dl.w_start)
            ranges.append((w_local_y, w_local_y, c, d))
            axial_ranges.append((w_local_x, w_local_x, c, d))
        elif dl.direction == 'global':
            # 全域任意角度載重(global_y的推廣版, 支援局部段+線性變化,
            # 見model.py的DistributedLoad說明), 跟dofmanager.py組裝時
            # 同一套分解邏輯: 兩端分開投影到局部x/y座標。
            ni, nj = frame.nodes[m_obj.node_i], frame.nodes[m_obj.node_j]
            _, angle = member_geometry(ni, nj)
            c_ang, s_ang = np.cos(angle), np.sin(angle)
            ang = np.radians(dl.angle_deg)
            u_global = np.array([np.cos(ang), np.sin(ang)])
            R = np.array([[c_ang, s_ang], [-s_ang, c_ang]])
            local_start = R @ (u_global * dl.w_start)
            local_end = R @ (u_global * dl.w_end)
            ranges.append((local_start[1], local_end[1], c, d))
            axial_ranges.append((local_start[0], local_end[0], c, d))
        else:
            ranges.append((dl.w_start, dl.w_end, c, d))

    x = np.linspace(0, L, n)
    eps = L * 1e-6
    extra = []
    for pl in frame.member_point_loads:
        if pl.member == member_id:
            extra += [max(0.0, pl.a - eps), min(L, pl.a + eps)]
    for _, _, c, d in ranges:
        extra += [max(0.0, c - eps), min(L, c + eps), max(0.0, d - eps), min(L, d + eps)]
    if extra:
        x = np.unique(np.concatenate([x, extra]))

    N = np.full(x.shape, -Fx1)   # 拉力為正的工程慣例 (Fx1本身是"壓力為正", 取負號校正)
    for w_start, w_end, c, d in axial_ranges:
        # 軸向均佈載重造成N(x)線性變化: N(x) = -Fx1 - (累加軸向載重)
        # (減號: Fx1本身已經是壓力為正的節點力慣例, 這裡的累加項要跟它同一套
        # 慣例, 才能先加總再一起取負號校正成拉力為正; 已用slop-roof案例
        # 對照SW FEA的N(x)逐點資料驗證過, 見tests/test_sloped_roof_global_udl.py)
        N -= _partial_load_W(w_start, w_end, c, d, x)
    W = np.zeros_like(x)
    Mcum = np.zeros_like(x)
    for w_start, w_end, c, d in ranges:
        W += _partial_load_W(w_start, w_end, c, d, x)
        Mcum += _partial_load_M(w_start, w_end, c, d, x)
    V = Fy1 + W
    M = -M1 + Fy1 * x + Mcum

    for pl in frame.member_point_loads:
        if pl.member != member_id:
            continue
        step = (x >= pl.a - 1e-12).astype(float)
        N += -pl.fx * step
        V += pl.fy * step
        M += pl.fy * np.clip(x - pl.a, 0, None)
        M += -pl.m * step

    return x, N, V, M


def member_deformed_shape(frame, result, member_id, scale=1.0, n=51):
    """回傳桿件變形後的全域座標 (X, Y) 陣列 (含放大係數scale)。

    精確作法(不是近似, 也不是另外的"精確版"): 直接對已經驗證過的 M(x)
    (member_internal_forces) 除以EI做兩次數值積分, 用兩端FEM算出來的真實
    節點位移釘住積分常數, 得到唯一一條v(x)撓度曲線。這樣桿件內部有集中力
    (member_point_load)或局部段均佈載重時撓度曲線依然精確——不會像單純
    兩端內插一條三次多項式(舊做法)那樣, 在桿件內部有集中力時可能低估撓度
    20%以上、位置也偏(見開發過程中實測的case: L=10,P=30,a=3.5, 舊內插法
    低估了約23%)。求解器本身(節點位移/反力/N,V,M)一直都是精確的, 這個
    修正只是讓「畫變形曲線」這一步的方法, 跟求解器同一個精確度等級,
    不需要另外維護一個「近似畫圖版」+「精確驗算版」。

    對truss/cable(沒有彎曲勁度, M(x)恆為0), 這個方法自動退化成兩端直線
    內插, 跟truss/cable「桿件維持直線不彎曲」的物理行為一致, 不需要
    另外特殊處理。
    """
    from .elements import member_geometry, transformation_matrix
    m = frame.members[member_id]
    node_i = frame.nodes[m.node_i]
    node_j = frame.nodes[m.node_j]
    L, angle = member_geometry(node_i, node_j)
    T = transformation_matrix(angle)

    dofs = frame.dofs_of(m.node_i) + frame.dofs_of(m.node_j)
    u_global_member = result.displacements[np.array(dofs)]
    u_local = T @ u_global_member
    u1, v1, th1, u2, v2, th2 = u_local

    section = frame.sections[m.section]
    EI = section.E * section.I

    if m.member_type in ('truss', 'cable') or abs(EI) < 1e-30:
        # 沒有彎曲勁度: 桿件維持直線, 兩端內插即為精確解
        x = np.linspace(0, L, n)
        Y_local = v1 + (v2 - v1) * (x / L)
    else:
        n_fine = max(n, 101)   # 積分用細網格(含載重/段落邊界附近的加密點), 保證精度
        x_fine, N_f, V_f, M_f = member_internal_forces(frame, result, member_id, n=n_fine)
        kappa = M_f / EI
        # 累積梯形積分: 曲率 -> 斜率(不定積分, 從0開始, 積分常數還沒定)
        slope_indef = np.concatenate([[0.0], np.cumsum((kappa[1:] + kappa[:-1]) / 2 * np.diff(x_fine))])
        # 斜率 -> 撓度(不定積分, 從0開始)
        defl_indef = np.concatenate([[0.0], np.cumsum((slope_indef[1:] + slope_indef[:-1]) / 2 * np.diff(x_fine))])
        # 邊界條件 v(0)=v1, v(L)=v2 (FEM算出的真實節點位移) 釘住兩個積分常數
        C2 = v1
        C1 = (v2 - v1 - defl_indef[-1]) / L
        Y_fine = defl_indef + C1 * x_fine + C2
        # 注意: x_fine不是均勻網格(member_internal_forces會在集中力/段落邊界
        # 附近插入額外取樣點, 保證那裡的跳躍畫得出來), 所以這裡改用內插對齊到
        # 均勻網格 x=linspace(0,L,n), 不能直接按"索引"下取樣(按索引下取樣會
        # 因為x_fine本身疏密不均, 導致取出來的x座標也跟著疏密不均、不等於
        # 呼叫端以為的均勻網格, 曾經因此在畫圖時對錯位置導致標籤數字離譜)。
        x = np.linspace(0, L, n)
        Y_local = np.interp(x, x_fine, Y_fine)

    X_local = x + scale * (u1 + (u2 - u1) * (x / L))   # 軸向位移線性內插(桿件軸向本來就是線性的)
    Y_local = scale * Y_local

    # 轉回全域座標 (局部->全域是 T的轉置, 再加上桿件原點座標)
    c, s = np.cos(angle), np.sin(angle)
    X_global = node_i.x + X_local * c - Y_local * s
    Y_global = node_i.y + X_local * s + Y_local * c
    return X_global, Y_global
