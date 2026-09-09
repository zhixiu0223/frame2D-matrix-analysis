"""
遞增側推(pushover)求解器: 位移控制 + event-to-event精確定位降伏事件。

移植自portal-frame-pushover-scratch/src/portal_frame_pushover.py(已用
OpenSeesPy驗證過, 單層框架全曲線誤差0.55%、雙層框架0.02%), 改寫成用
frame2d.dofmanager.build_dof_map()的通用DOF系統(支援release端專屬DOF、
任意拓樸), 取代原本寫死3*node_idx的簡化版Frame。

架構(跟原版一致):
  外層: 位移控制增量迴圈(每步名目步長d_nominal, 逼近target_total)
  中層: 在目前這段"分段線性"(鉸狀態固定不變)範疇下, 一個增量只需要
        解一次線性方程式就精確收斂——目前只有雙折線鉸, 沒有連續硬化
        曲線, 不需要真正的Newton-Raphson多次迭代(跟原版說明一致)
  內層: event-to-event子步長 -- 每步先用當前(常數)切線剛度試算, 檢查
        是否有尚未降伏的鉸會在這步內被觸發, 若會, 精確反算子步長, 優先
        走到事件點、更新鉸狀態, 再繼續剩下的步長

跟原始版本的關鍵差異(配合frame2d架構的必要調整, 不是簡化):
  - K組裝重用hinge.full_6x6_with_hinge() + dofmanager.build_dof_map(),
    支援release端(額外DOF), 不是寫死3*node_idx
  - 重力預載直接重用已驗證過的dofmanager.solve_with_hinges()做一次
    力控制求解, 不用另外寫一套力控制求解邏輯
  - P-Delta軸力代表值取該元素node_j端(index3), 跟frame2d自己的
    「拉力為正」慣例一致(不是原版的index0再取負號)

跟webapi/跟webapi_stdlib/兩份diagrams.py刻意重複同一個道理: 這裡的
_assemble_stiffness_with_hinges()獨立複製一份組裝邏輯, 不是重構
dofmanager.py共用, 避免一般線性/P-Delta求解路徑的維護意外牽動側推邏輯,
反之亦然。
"""
import numpy as np

from .hinge import full_6x6_with_hinge
from .elements import member_geometry, member_stiffness_local, member_stiffness_local_truss, transformation_matrix
from .dofmanager import build_dof_map, solve_with_hinges

M_LOCAL_IDX = [2, 5]   # 局部力向量[Fx1,Fy1,M1,Fx2,Fy2,M2]裡, 兩端彎矩的索引


def _fixed_dof_set(frame):
    """回傳側推增量分析用的「固定」DOF集合: frame.supports裡有設定值
    (不管是True的0.0固定還是settlement的非零值)的節點DOF——側推是位移
    「增量」分析, 支承的絕對位置在重力階段就已經定了, 增量分析裡支承的
    位移增量永遠是0, 不管它原本的絕對值是不是0。"""
    fixed = set()
    for s in frame.supports:
        ux, uy, rot = frame.dofs_of(s.node)
        for dof_idx, val in zip((ux, uy, rot), (s.ux, s.uy, s.rot)):
            if val is not None:
                fixed.add(dof_idx)
    return fixed


def _member_geometry_updated(frame, m, u_full_cum):
    """跟elements.member_geometry()一樣算L跟angle, 但用"目前累積變形後"
    的節點座標(原始座標+u_full_cum裡對應的ux,uy), 不是frame.nodes裡
    固定不變的原始座標。這是"幾何更新"(geometry-updating)側推模式的
    核心: 每一步都用上一步結束時真正的變形位置重新算桿件幾何, 不是
    全程都凍結在最初始的未變形幾何上——後者(凍結幾何, 只疊加跟軸力
    成正比的Kg修正項)就是預設的線性化P-Delta做法, 這裡是另一種模式。

    已知限制(誠實記錄, 不誇大): 這裡只更新了"用哪個幾何組裝勁度矩陣"
    這一件事(Updated Lagrangian的基本精神), 沒有做完整的co-rotational
    大轉角處理(例如沒有把桿件自己的剛體轉動從局部變形裡分離出來),
    也沒有在每個增量步內做真正的Newton-Raphson平衡疊代去消除殘餘力
    (跟不開這個選項時同一套"每步凍結切線剛度、只解一次線性方程式"
    的疊代格式一樣, 只是這次切線剛度用的幾何會更新)。對變形沒有大到
    需要考慮桿件本身大幅轉動的情況, 這個近似已經比完全凍結幾何更接近
    真實行為; 但不是完整意義上的大變形非線性分析。"""
    ni = frame.nodes[m.node_i]
    nj = frame.nodes[m.node_j]
    dofs_i = frame.dofs_of(m.node_i)
    dofs_j = frame.dofs_of(m.node_j)
    xi = ni.x + u_full_cum[dofs_i[0]]
    yi = ni.y + u_full_cum[dofs_i[1]]
    xj = nj.x + u_full_cum[dofs_j[0]]
    yj = nj.y + u_full_cum[dofs_j[1]]
    L = np.hypot(xj - xi, yj - yi)
    angle = np.arctan2(yj - yi, xj - xi)
    return L, angle


def _assemble_stiffness_with_hinges(frame, hinge_states, axial_forces=None, u_full_cum=None):
    """組裝含鉸(+選用P-Delta)的全域勁度矩陣, 純粹供側推增量分析用,
    不含任何載重/邊界條件。

    axial_forces: 可選, {member_id: N(拉力為正)}, 用於P-Delta。
    u_full_cum: 可選, 目前累積的絕對位移向量——給了這個, 桿件的L跟
        angle會用"目前變形後"的節點位置重新算(見_member_geometry_updated()),
        不是frame.nodes裡固定的原始座標。這是"幾何更新"模式用的參數,
        預設None時完全等同於原本行為(用原始未變形幾何), 不影響任何
        既有呼叫端。

    回傳: K, member_dofs, member_T, member_L(供後續增量內力回算用)。
    """
    member_dofs, n_node_dof, n_extra_dof = build_dof_map(frame)
    n = n_node_dof + n_extra_dof
    K = np.zeros((n, n))
    member_T = {}
    member_L = {}
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        if u_full_cum is not None:
            L, angle = _member_geometry_updated(frame, m, u_full_cum)
        else:
            node_i = frame.nodes[m.node_i]
            node_j = frame.nodes[m.node_j]
            L, angle = member_geometry(node_i, node_j)
        T = transformation_matrix(angle)
        member_T[mid] = T
        member_L[mid] = L
        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            P = 0.0 if axial_forces is None else axial_forces.get(mid, 0.0)
            if hinge_states is not None and mid in hinge_states:
                k_local = full_6x6_with_hinge(section.E, section.A, section.I, L, hinge_states[mid], P=P)
            else:
                k_local = member_stiffness_local(section.E, section.I, section.A, L, P=P)
        k_global = T.T @ k_local @ T
        idx = np.array(member_dofs[mid])
        K[np.ix_(idx, idx)] += k_global
    return K, member_dofs, member_T, member_L


def solve_displacement_increment(K, prescribed_dofs, du_prescribed, fixed_dofs):
    """給定全域K、受控自由度清單(全域dof編號)跟其位移增量, 解出所有
    自由度的位移增量(固定=0, 受控=給定, 其餘=求解)。跟
    portal-frame-pushover-scratch的Frame.solve_increment()同一個演算法。"""
    n = K.shape[0]
    all_dofs = set(range(n))
    prescribed_set = set(prescribed_dofs)
    free_dofs = sorted(all_dofs - fixed_dofs - prescribed_set)
    prescribed_dofs = list(prescribed_dofs)

    du_full = np.zeros(n)
    for d, val in zip(prescribed_dofs, du_prescribed):
        du_full[d] = val

    if len(free_dofs) > 0:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        K_fp = K[np.ix_(free_dofs, prescribed_dofs)]
        rhs = -K_fp @ np.array(du_prescribed)
        try:
            du_free = np.linalg.solve(K_ff, rhs)
        except np.linalg.LinAlgError:
            raise RuntimeError("勁度矩陣奇異, 無法求解增量 -- 結構可能已經達到機構狀態。")
        for d, val in zip(free_dofs, du_free):
            du_full[d] = val
    return du_full


def solve_force_increment(K, load_dofs, dF_prescribed, fixed_dofs):
    """力控制版的增量求解: 給定力增量dF(施加在load_dofs上), 解出所有
    非固定自由度的位移增量。跟位移控制不同的地方: 力控制沒有"受控但
    不求解"這個中間類別, 除了固定的自由度以外全部都是要解的自由度,
    load_dofs只是"外力施加在哪裡", 不是"這個自由度的位移已知"。

    這個函式沒有辦法解過極限承載力(peak)之後的軟化段——一旦切線
    剛度矩陣不再正定, K_ff會奇異或病態, 直接丟RuntimeError, 這是力
    控制方法本身的極限, 不是bug(這也是run_pushover()預設用位移控制
    的原因: 位移控制在軟化段依然穩定可解, 力控制不行)。這個函式主要
    用途是驗證彈性範圍內的結果對不對得上手算, 不是用來找極限承載力。
    """
    n = K.shape[0]
    all_dofs = set(range(n))
    free_dofs = sorted(all_dofs - fixed_dofs)

    dF_full = np.zeros(n)
    for d, val in zip(load_dofs, dF_prescribed):
        dF_full[d] += val

    du_full = np.zeros(n)
    if len(free_dofs) > 0:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        rhs = dF_full[free_dofs]
        try:
            du_free = np.linalg.solve(K_ff, rhs)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                "勁度矩陣奇異或接近奇異, 無法求解增量 -- 結構可能已經達到或"
                "超過極限承載力(力控制無法穿越peak之後的軟化段, 建議改用"
                "位移控制才能繼續往下推)。")
        for d, val in zip(free_dofs, du_free):
            du_full[d] = val
    return du_full


def _member_force_increments(frame, member_T, member_dofs, member_L, hinge_states, axial_forces, du_full):
    """給定全域位移增量, 回傳每個元素局部座標下的內力增量
    {member_id: 6-vector}(對truss/cable也一併算出, 只是它們沒有塑鉸,
    不會被_find_crossing_events()用到)。"""
    result = {}
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
        L = member_L[mid]
        T = member_T[mid]
        du_local = T @ du_full[np.array(member_dofs[mid])]
        if m.member_type in ('truss', 'cable'):
            k_local = member_stiffness_local_truss(section.E, section.A, L)
        else:
            P = 0.0 if axial_forces is None else axial_forces.get(mid, 0.0)
            if hinge_states is not None and mid in hinge_states:
                k_local = full_6x6_with_hinge(section.E, section.A, section.I, L, hinge_states[mid], P=P)
            else:
                k_local = member_stiffness_local(section.E, section.I, section.A, L, P=P)
        result[mid] = k_local @ du_local
    return result


def _find_crossing_events(hinge_states, cum_forces, df_local_by_member):
    """檢查這個試探步是否有任何尚未降伏的鉸端點會被觸發。回傳
    [(ratio, member_id, end_idx), ...], ratio是這一步內事件發生的比例
    (0~1)。"""
    events = []
    for mid, hs in hinge_states.items():
        df_local = df_local_by_member[mid]
        for end_idx in (0, 1):
            if hs.yielded[end_idx]:
                continue
            M_idx = M_LOCAL_IDX[end_idx]
            M_cum = cum_forces[mid][M_idx]
            M_trial = M_cum + df_local[M_idx]
            Mp = hs.Mp[end_idx]
            if abs(M_trial) >= Mp:
                denom = abs(M_trial) - abs(M_cum)
                if denom <= 0:
                    continue   # 數值上不該發生, 跳過避免除以零
                ratio = (Mp - abs(M_cum)) / denom
                ratio = max(0.0, min(1.0, ratio))
                events.append((ratio, mid, end_idx))
    return events


def _update_theta_p(frame, hinge_states, member_T, member_dofs, member_L, u_full):
    """用目前累積的絕對位移u_full, 算每個有塑鉸的frame元素兩端"樑內部
    端點轉角跟外部節點轉角的差"(=彈簧的相對轉動), 更新
    HingeState.theta_p。只對已經降伏的端點更新(未降伏時彈簧幾乎不轉動,
    RIGID_FACTOR夠大, 直接當作0, 不用另外算);已降伏端點的彈簧轉動幾乎
    全部是塑性變形(降伏前的彈性轉動已經因為RIGID_FACTOR趨近0而可忽略),
    這是集中塑性鉸(lumped plasticity)理論標準的近似, 不是另外發明的
    公式。純粹更新"事後可查詢"的欄位, 不影響任何求解邏輯本身。"""
    from .hinge import beam_internal_rotation
    for mid, hs in hinge_states.items():
        if not any(hs.yielded):
            continue
        m = frame.members[mid]
        section = frame.sections[m.section]
        L = member_L[mid]
        T = member_T[mid]
        u_local = T @ u_full[np.array(member_dofs[mid])]
        v1, theta1, v2, theta2 = u_local[1], u_local[2], u_local[4], u_local[5]
        phi1, phi2 = beam_internal_rotation(section.E, section.I, L, hs, v1, theta1, v2, theta2)
        if hs.yielded[0]:
            hs.theta_p[0] = abs(theta1 - phi1)
        if hs.yielded[1]:
            hs.theta_p[1] = abs(theta2 - phi2)


def check_mechanism(K, prescribed_dofs, fixed_dofs, stiffness_ratio_limit=1e-8):
    """機構偵測: 算受控自由度(側推方向)的縮聚剛度(Schur complement),
    這才是「再推一單位位移, 需要多少額外力」的真實剛度。機構形成時這個
    縮聚剛度會趨近奇異, 才是「已經推不動」的正確判準(不能只看K_ff是否
    病態, 樑柱接頭轉角等自由度在側推方向以外通常還有勁度, 不會因為
    側推方向已經形成機構就跟著病態——這是portal-frame-pushover-scratch
    開發過程中修正過的錯誤, 這裡直接沿用修正後的正確版本)。

    回傳True表示「已達機構, 應停止」。"""
    n = K.shape[0]
    all_dofs = set(range(n))
    free_dofs = sorted(all_dofs - fixed_dofs - set(prescribed_dofs))
    p_dofs = list(prescribed_dofs)

    K_pp = K[np.ix_(p_dofs, p_dofs)]
    if len(free_dofs) == 0:
        K_condensed = K_pp
    else:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        K_fp = K[np.ix_(free_dofs, p_dofs)]
        K_pf = K[np.ix_(p_dofs, free_dofs)]
        try:
            K_ff_inv_K_fp = np.linalg.solve(K_ff, K_fp)
        except np.linalg.LinAlgError:
            return True
        K_condensed = K_pp - K_pf @ K_ff_inv_K_fp

    try:
        eigvals = np.linalg.eigvalsh((K_condensed + K_condensed.T) / 2.0)
    except np.linalg.LinAlgError:
        return True

    if np.any(eigvals <= 0):
        return True
    ratio = eigvals.min() / eigvals.max()
    return ratio < stiffness_ratio_limit


def apply_gravity(frame, hinge_states):
    """重力預載階段: 直接重用已驗證過的dofmanager.solve_with_hinges(),
    對frame目前的point_loads/distributed_loads做一次力控制的線性求解
    (用給定的hinge_states, 通常是全部未降伏), 取每個frame元素的局部
    端點內力當作側推階段的起始累積內力。

    回傳: cum_forces(dict member_id -> 6-vector局部端點內力, 涵蓋所有
    frame元素), result(完整的SolveResult, 供查詢重力階段本身的位移/
    反力用)。"""
    result = solve_with_hinges(frame, hinge_states)
    cum_forces = {mid: np.array(mr.end_forces_local, dtype=float)
                  for mid, mr in result.member_results.items()}
    return cum_forces, result


def _snapshot(cum_forces, hinge_states):
    """給定目前的cum_forces跟hinge_states, 回傳一份跟numpy/HingeState物件
    完全脫鉤的純Python/list快照(不會被之後的原地修改牽動)——這是逐步
    回放要顯示"當下彎矩分佈/塑鉸狀態"的資料來源。"""
    return {
        'member_forces': {mid: [float(v) for v in f] for mid, f in cum_forces.items()},
        'hinge_states': {mid: {'yielded': list(hs.yielded), 'theta_p': [float(t) for t in hs.theta_p]}
                          for mid, hs in hinge_states.items()},
    }


def run_pushover(frame, hinge_states, prescribed_dofs, direction, target_total,
                  d_nominal, base_reaction_dofs, initial_cum_forces=None,
                  use_pdelta=False, mechanism_ratio_limit=1e-8, max_steps=100000,
                  include_final_displacement=False, include_snapshots=False,
                  control_mode='displacement', geometry_update=False,
                  include_final_reactions=False, include_max_rotation=False):
    """遞增側推主迴圈, 支援位移控制(預設)或力控制。

    frame: 已定義節點/桿件/支承的Frame2D(側推的推力來自prescribed_dofs
        被強制位移或施加力增量, 不是外加point_load; 重力階段的
        point_load在apply_gravity()另外處理)
    hinge_states: {member_id: HingeState}, 通常從hinge.initial_hinge_states(frame)
        取得, 也可以是apply_gravity()跑完後(仍未降伏)的同一組物件——這個
        函式會直接原地修改傳進來的HingeState(設定yielded), 呼叫端要自己
        決定要不要事先複製一份
    prescribed_dofs: 受控自由度清單(全域dof編號, 例如[frame.dofs_of(roof_node)[0]])——
        位移控制時是"被強制位移的自由度", 力控制時是"外力施加的自由度"
    direction: 對應prescribed_dofs的方向係數(例如[1.0]表示單點控制, 或
        [1.0, 1.0]表示兩個控制點同向等量推/等量施力)
    target_total: 位移控制時是目標總側推位移量; 力控制時是目標總施加力
    d_nominal: 名目步長(沒有事件發生時每步走多少); 位移控制時是位移
        增量, 力控制時是力增量
    base_reaction_dofs: 用來加總算底剪力的自由度(例如各支承節點的水平DOF)
    initial_cum_forces: apply_gravity()算出的初始元素內力, None代表從零
        開始(無重力預載)
    use_pdelta: 是否把frame元素目前軸力組進幾何剛度矩陣(P-Delta), 軸力
        代表值取該元素node_j端(拉力為正)
    mechanism_ratio_limit: 位移控制下, 縮聚剛度最小/最大特徵值比例低於
        此值視為已達機構(力控制不用這個判斷, 見下面control_mode說明)
    max_steps: 安全閥(避免d_nominal設太小或有bug時無窮迴圈), 超過直接
        raise, 不會被誤判成"正常跑完"
    include_final_displacement: False(預設)時回傳5個值, 完全等同原本的
        呼叫方式跟既有測試不受影響。True時額外多回傳最終的絕對位移向量
        u_full_cum(全域DOF編號, 跟frame.dofs_of()對應), 供後續需要畫
        最終變形形狀(例如疊加塑鉸圈圈)的呼叫端使用。
    include_snapshots: False(預設)時不影響任何既有行為。True時額外多
        回傳history_snapshots(list, 跟history_u/history_F逐一對應,
        第0筆是u=0的初始狀態), 每筆是{'member_forces': {member_id:
        [Fx1,Fy1,M1,Fx2,Fy2,M2]}, 'hinge_states': {member_id:
        {'yielded':[bool,bool], 'theta_p':[float,float]}}}, 供"逐步
        回放"這種需要看每一步當下彎矩分佈/塑鉸狀態的功能使用。
    control_mode: 'displacement'(預設, 完全等同原本行為)或'force'。
        force模式下沒有辦法解過極限承載力之後的軟化段(切線剛度矩陣
        一旦不再正定就會直接RuntimeError被這裡包成"已達機構"提前結束,
        不是bug——這是力控制方法本身的極限), 主要用途是驗證彈性範圍
        內的結果直接對得上手算(施加已知的力, 解出轉角, 不用像位移
        控制那樣先猜位移再看對應多少力)。force模式下不會呼叫
        check_mechanism()(那是針對位移控制的"受控vs自由"DOF切分設計
        的判斷, 力控制沒有這個切分——所有非固定自由度都是要解的自由度)。
    geometry_update: False(預設, 完全等同原本行為, 全程用最初始的
        未變形幾何組裝勁度矩陣)或True(每一步用"目前累積變形後"的節點
        位置重新算桿件長度/角度, 見_member_geometry_updated())。這是
        Updated Lagrangian的基本精神(用當下真實變形的幾何建立切線
        剛度), 比預設的線性化P-Delta(幾何凍結、只疊加軸力修正項)更
        接近大變形時的真實行為, 但不是完整的co-rotational大轉角分析
        (沒有把桿件自己的剛體轉動獨立分離出來處理, 每個增量步也還是
        只解一次線性方程式、不是真正疊代到殘餘力歸零的Newton-Raphson)。
        兩種模式在變形不大時應該給出幾乎相同的結果(幾何本來就沒怎麼變),
        變形越大差異會越明顯。

    回傳: history_u(np.array, 這裡永遠是prescribed_dofs[0]這個自由度
        實際的累積位移, 不管control_mode是哪一種——位移控制下這個值
        直接等於施加的位移增量累加; 力控制下這個值是解出來的結果),
        history_F(np.array, 底剪力), event_log(list of dict),
        hinge_states(原地更新後的同一組物件), mechanism_reached(bool)
        [, u_full_cum(np.array) -- 只有include_final_displacement=True時]
        [, history_snapshots(list) -- 只有include_snapshots=True時]
        [, cum_reaction(np.array) -- 只有include_final_reactions=True時,
          全域自由度反力向量, 供畫最終狀態的結構/內力圖用]
        [, max_rotation(float) -- 只有include_max_rotation=True時, 所有
          節點裡最大的轉角絕對值(rad)——超過大約0.1rad(約5.7度)代表
          這次結果已經超出小角度假設的有效範圍, 不管有沒有開P-Delta/
          geometry_update都不可信, 見對話紀錄裡實際案例逼出來的發現]
        (以上四個都是選用的, 依序append在後面, 不管哪個組合都是同一個
        固定順序: final_displacement, snapshots, final_reactions,
        max_rotation)
    """
    direction = np.array(direction, dtype=float)
    fixed_dofs = _fixed_dof_set(frame)
    _, n_node_dof, n_extra_dof = build_dof_map(frame)
    n_total = n_node_dof + n_extra_dof
    control_dof0 = prescribed_dofs[0]

    if initial_cum_forces is None:
        cum_forces = {mid: np.zeros(6) for mid in frame.members}
    else:
        cum_forces = {mid: np.array(f, dtype=float) for mid, f in initial_cum_forces.items()}
    cum_reaction = np.zeros(n_total)

    u_control_cum = 0.0   # prescribed_dofs[0]的實際累積位移(兩種控制模式
                           # 都從du_full裡直接讀這個自由度的值, 不是假設
                           # 它等於施加的增量本身——位移控制下兩者剛好
                           # 相等, 力控制下這個值本來就是解出來的結果)
    history_u = [0.0]
    history_F = [0.0]
    event_log = []
    mechanism_reached = False
    history_snapshots = [_snapshot(cum_forces, hinge_states)] if include_snapshots else None
    u_full_cum = np.zeros(n_total)   # 累積絕對位移, 給_update_theta_p()算樑
                                      # 內部端點真正轉角用(不能只累積增量,
                                      # 因為beam_internal_rotation()的公式
                                      # 需要目前絕對的v1,theta1,v2,theta2)

    def current_axial_forces():
        return {mid: f[3] for mid, f in cum_forces.items()}   # 拉力為正(端j的Fx)

    def solve_increment(K, amount):
        if control_mode == 'force':
            return solve_force_increment(K, prescribed_dofs, direction * amount, fixed_dofs)
        return solve_displacement_increment(K, prescribed_dofs, direction * amount, fixed_dofs)

    remaining = target_total
    steps = 0
    while remaining > 1e-9:
        steps += 1
        if steps > max_steps:
            raise RuntimeError(
                f"側推超過{max_steps}步仍未達到target_total, 可能是d_nominal"
                "設太小或有其他問題, 已中止(避免無窮迴圈)。")
        d_step = min(d_nominal, remaining)
        axial = current_axial_forces() if use_pdelta else None
        K, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(
            frame, hinge_states, axial_forces=axial,
            u_full_cum=(u_full_cum if geometry_update else None))

        if control_mode == 'displacement' and check_mechanism(K, prescribed_dofs, fixed_dofs, mechanism_ratio_limit):
            mechanism_reached = True
            break

        try:
            du_full = solve_increment(K, d_step)
        except RuntimeError:
            if control_mode == 'force':
                mechanism_reached = True
                break
            raise
        df_by_member = _member_force_increments(frame, member_T, member_dofs, member_L,
                                                 hinge_states, axial, du_full)
        events = _find_crossing_events(hinge_states, cum_forces, df_by_member)

        if events:
            ratio_min = min(r for r, _, _ in events)
            d_sub = ratio_min * d_step
            du_full_sub = solve_increment(K, d_sub)
            df_by_member_sub = _member_force_increments(frame, member_T, member_dofs, member_L,
                                                          hinge_states, axial, du_full_sub)

            for mid, df in df_by_member_sub.items():
                cum_forces[mid] += df
            cum_reaction += K @ du_full_sub
            u_full_cum += du_full_sub
            u_control_cum += du_full_sub[control_dof0]
            remaining -= d_sub

            newly_yielded = []
            for r, mid, end_idx in events:
                if abs(r - ratio_min) < 1e-9:
                    hinge_states[mid].yielded[end_idx] = True
                    newly_yielded.append((mid, end_idx))

            _update_theta_p(frame, hinge_states, member_T, member_dofs, member_L, u_full_cum)

            F_base = -sum(cum_reaction[d] for d in base_reaction_dofs)
            history_u.append(u_control_cum)
            history_F.append(F_base)
            event_log.append({'u': u_control_cum, 'F': F_base, 'yielded': newly_yielded})
            if include_snapshots:
                history_snapshots.append(_snapshot(cum_forces, hinge_states))
        else:
            for mid, df in df_by_member.items():
                cum_forces[mid] += df
            cum_reaction += K @ du_full
            u_full_cum += du_full
            u_control_cum += du_full[control_dof0]
            remaining -= d_step

            _update_theta_p(frame, hinge_states, member_T, member_dofs, member_L, u_full_cum)

            F_base = -sum(cum_reaction[d] for d in base_reaction_dofs)
            history_u.append(u_control_cum)
            history_F.append(F_base)
            if include_snapshots:
                history_snapshots.append(_snapshot(cum_forces, hinge_states))

    # 轉角自由度是dofs_of()回傳的第三個(index 2), 也就是每個節點的
    # 全域dof編號裡mod 3餘2的那些——檢查有沒有任何節點轉角超過小角度
    # 假設的合理範圍(0.1rad約5.7度): 我們的彎曲公式(跟所有標準線彈性
    # 梁元素一樣)是建立在"轉角很小"這個前提上的, 一旦轉角真的變大,
    # 不管有沒有開geometry_update, 結果都已經超出模型的有效範圍——這是
    # 真實案例逼出來的發現(見對話紀錄), 不是憑空假設的邊界條件。
    max_rotation = float(np.max(np.abs(u_full_cum[2:n_node_dof:3]))) if n_node_dof >= 3 else 0.0

    result = [np.array(history_u), np.array(history_F), event_log, hinge_states, mechanism_reached]
    if include_final_displacement:
        result.append(u_full_cum)
    if include_snapshots:
        result.append(history_snapshots)
    if include_final_reactions:
        result.append(cum_reaction)
    if include_max_rotation:
        result.append(max_rotation)
    return tuple(result)


# ============================================================
# run_pushover_converged() -- 幾何平衡疊代版本, 跟run_pushover()是完全
# 獨立的兩個函式(刻意不共用主迴圈), 對話紀錄裡使用者/ChatGPT都明確建議
# 不要改掉已經驗證過的run_pushover(), 而是新增一個平行的版本, 讓兩者
# 可以互相對照。
#
# 差在哪裡(誠實記錄, 不誇大):
#   run_pushover(): 塑鉸狀態凍結的每一段區間內, 只用"這一段開始時"的
#     幾何/軸力組一次勁度矩陣、解一次線性方程式, 不管有沒有開
#     geometry_update, 都不會回頭檢查"這個答案在真正的變形終點上,
#     幾何/軸力是否還跟一開始用的一致"。
#   run_pushover_converged(): 同一段凍結區間內, 反覆疊代
#     "用目前試探的位移場重算幾何/軸力 -> 重新組K -> 重新解這段增量"
#     直到位移增量不再明顯變化(用相對範數當收斂判斷), 才承認這一段
#     真的解完。這是Picard(不動點)疊代, 不是完整的Newton-Raphson
#     (沒有算tangent的解析導數去加速收斂), 但確實會在接近極限承載力、
#     真正的幾何非線性效應顯著時, 表現出"疊代不收斂"這個明確訊號,
#     不會像run_pushover()那樣安靜地給出一條可能已經失真的曲線。
#
#   材料非線性(塑鉸降伏)這一塊完全沿用event-to-event的邏輯跟公式
#   (_find_crossing_events()等), 沒有另外設計——這裡疊代的對象只有
#   幾何/軸力這一塊, 不是"材料+幾何"完整耦合的殘餘力形式Newton-Raphson
#   (那需要把雙折線塑鉸重新設計成連續可微分的形式, 是更大的工程,
#   這裡沒有做)。
# ============================================================

def _solve_step_geom_converged(frame, hinge_states, cum_forces, u_full_cum,
                                prescribed_dofs, direction, d_amount, fixed_dofs,
                                control_mode, use_pdelta, geometry_update,
                                geom_tol, max_geom_iter):
    """對"這一步的位移/力增量大小是d_amount"這件事, 疊代到幾何/軸力
    自洽為止。回傳(converged, du_full, df_by_member, K, member_T,
    member_dofs, member_L)——K是最後一次(收斂時)那次疊代組出來的勁度
    矩陣, 呼叫端用同一個K去更新反力, 不要另外重新組一次(重組的話,
    如果用來重組的幾何/軸力狀態跟疊代收斂時不一致, 反力會對不上實際
    被接受的du_full, 這是一個容易犯的錯, 這裡直接把K傳出去避免)。

    疊代方式: 用目前試探的位移增量du_trial(初始猜測=0, 也就是先用
    這一段"開始時"的幾何/軸力算第一次)算出對應的幾何(u_full_cum+
    du_trial, 只有geometry_update=True時採用)跟軸力(cum_forces的
    軸力分量+這次試探增量算出的軸力分量, 只有use_pdelta=True時採用),
    重新組K, 重新解出新的du。如果新舊du的相對範數差距小於geom_tol,
    視為收斂;連續max_geom_iter次都沒收斂, 回傳converged=False,
    呼叫端要自己決定怎麼處理(這裡的設計是: 不收斂就視為已經到極限,
    優雅停止, 不會硬給一個不可信的答案)。
    """
    du_trial = np.zeros_like(u_full_cum)
    df_trial = None
    K = member_T = member_dofs = member_L = None
    for _ in range(max_geom_iter):
        axial = None
        if use_pdelta:
            axial = {mid: f[3] for mid, f in cum_forces.items()}
            if df_trial is not None:
                for mid, df in df_trial.items():
                    axial[mid] = axial[mid] + df[3]
        geom_u = (u_full_cum + du_trial) if geometry_update else None
        K, member_dofs, member_T, member_L = _assemble_stiffness_with_hinges(
            frame, hinge_states, axial_forces=axial, u_full_cum=geom_u)
        if control_mode == 'force':
            du_new = solve_force_increment(K, prescribed_dofs, direction * d_amount, fixed_dofs)
        else:
            du_new = solve_displacement_increment(K, prescribed_dofs, direction * d_amount, fixed_dofs)
        df_new = _member_force_increments(frame, member_T, member_dofs, member_L, hinge_states, axial, du_new)

        denom = np.linalg.norm(du_new)
        diff = np.linalg.norm(du_new - du_trial) / denom if denom > 1e-14 else 0.0
        du_trial = du_new
        df_trial = df_new
        if diff < geom_tol:
            return True, du_trial, df_trial, K, member_T, member_dofs, member_L
    return False, du_trial, df_trial, K, member_T, member_dofs, member_L


def run_pushover_converged(frame, hinge_states, prescribed_dofs, direction, target_total,
                            d_nominal, base_reaction_dofs, initial_cum_forces=None,
                            use_pdelta=True, geometry_update=True,
                            mechanism_ratio_limit=1e-8, max_steps=100000,
                            geom_tol=1e-6, max_geom_iter=30,
                            control_mode='displacement',
                            include_final_displacement=False, include_snapshots=False,
                            include_final_reactions=False, include_max_rotation=False):
    """run_pushover()的幾何平衡疊代版本, 見上面模組層級的說明註解。
    參數跟run_pushover()大致對應, 差異只列這裡沒有的/新增的:

    use_pdelta/geometry_update: 預設都是True(這個函式存在的目的就是
        正確處理這兩件事的耦合, 不像run_pushover()是預設關閉當作
        選配)。
    geom_tol: 幾何/軸力疊代的相對收斂容忍度(無因次), 預設1e-6。
    max_geom_iter: 每一段最多疊代幾次, 超過視為這一段解不出來。
    (沒有control_mode='force'時的check_mechanism()判斷——力控制在
    這個版本裡, 疊代不收斂本身就是"已經到極限"的訊號, 不需要另外用
    縮聚剛度特徵值判斷。)

    回傳: 跟run_pushover()同樣的欄位, 但mechanism_reached在這個版本
    裡的意義稍微不同: True代表某一段幾何疊代沒有在max_geom_iter內
    收斂(不是位移控制下縮聚剛度接近奇異的那種判斷)——這是這個版本
    真正的"收斂失敗"訊號, 出現時所有更晚的位移/力都不會被加進歷程,
    在那個點提前停止。
    """
    direction = np.array(direction, dtype=float)
    fixed_dofs = _fixed_dof_set(frame)
    _, n_node_dof, n_extra_dof = build_dof_map(frame)
    n_total = n_node_dof + n_extra_dof
    control_dof0 = prescribed_dofs[0]

    if initial_cum_forces is None:
        cum_forces = {mid: np.zeros(6) for mid in frame.members}
    else:
        cum_forces = {mid: np.array(f, dtype=float) for mid, f in initial_cum_forces.items()}
    cum_reaction = np.zeros(n_total)

    u_control_cum = 0.0
    history_u = [0.0]
    history_F = [0.0]
    event_log = []
    mechanism_reached = False
    history_snapshots = [_snapshot(cum_forces, hinge_states)] if include_snapshots else None
    u_full_cum = np.zeros(n_total)

    remaining = target_total
    steps = 0
    while remaining > 1e-9:
        steps += 1
        if steps > max_steps:
            raise RuntimeError(
                f"側推超過{max_steps}步仍未達到target_total, 可能是d_nominal"
                "設太小或有其他問題, 已中止(避免無窮迴圈)。")
        d_step = min(d_nominal, remaining)

        converged, du_full, df_by_member, K_step, member_T, member_dofs, member_L = _solve_step_geom_converged(
            frame, hinge_states, cum_forces, u_full_cum, prescribed_dofs, direction, d_step,
            fixed_dofs, control_mode, use_pdelta, geometry_update, geom_tol, max_geom_iter)
        if not converged:
            mechanism_reached = True
            break

        events = _find_crossing_events(hinge_states, cum_forces, df_by_member)

        if events:
            ratio_min = min(r for r, _, _ in events)
            d_sub = ratio_min * d_step
            converged_sub, du_full_sub, df_by_member_sub, K_sub, member_T_sub, member_dofs_sub, member_L_sub = \
                _solve_step_geom_converged(
                    frame, hinge_states, cum_forces, u_full_cum, prescribed_dofs, direction, d_sub,
                    fixed_dofs, control_mode, use_pdelta, geometry_update, geom_tol, max_geom_iter)
            if not converged_sub:
                mechanism_reached = True
                break

            for mid, df in df_by_member_sub.items():
                cum_forces[mid] += df
            cum_reaction += K_sub @ du_full_sub
            u_full_cum += du_full_sub
            u_control_cum += du_full_sub[control_dof0]
            remaining -= d_sub

            newly_yielded = []
            for r, mid, end_idx in events:
                if abs(r - ratio_min) < 1e-9:
                    hinge_states[mid].yielded[end_idx] = True
                    newly_yielded.append((mid, end_idx))

            _update_theta_p(frame, hinge_states, member_T_sub, member_dofs_sub, member_L_sub, u_full_cum)

            F_base = -sum(cum_reaction[d] for d in base_reaction_dofs)
            history_u.append(u_control_cum)
            history_F.append(F_base)
            event_log.append({'u': u_control_cum, 'F': F_base, 'yielded': newly_yielded})
            if include_snapshots:
                history_snapshots.append(_snapshot(cum_forces, hinge_states))
        else:
            for mid, df in df_by_member.items():
                cum_forces[mid] += df
            cum_reaction += K_step @ du_full
            u_full_cum += du_full
            u_control_cum += du_full[control_dof0]
            remaining -= d_step

            _update_theta_p(frame, hinge_states, member_T, member_dofs, member_L, u_full_cum)

            F_base = -sum(cum_reaction[d] for d in base_reaction_dofs)
            history_u.append(u_control_cum)
            history_F.append(F_base)
            if include_snapshots:
                history_snapshots.append(_snapshot(cum_forces, hinge_states))

    max_rotation = float(np.max(np.abs(u_full_cum[2:n_node_dof:3]))) if n_node_dof >= 3 else 0.0

    result = [np.array(history_u), np.array(history_F), event_log, hinge_states, mechanism_reached]
    if include_final_displacement:
        result.append(u_full_cum)
    if include_snapshots:
        result.append(history_snapshots)
    if include_final_reactions:
        result.append(cum_reaction)
    if include_max_rotation:
        result.append(max_rotation)
    return tuple(result)
