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


def _assemble_stiffness_with_hinges(frame, hinge_states, axial_forces=None):
    """組裝含鉸(+選用P-Delta)的全域勁度矩陣, 純粹供側推增量分析用,
    不含任何載重/邊界條件。

    axial_forces: 可選, {member_id: N(拉力為正)}, 用於P-Delta。

    回傳: K, member_dofs, member_T, member_L(供後續增量內力回算用)。
    """
    member_dofs, n_node_dof, n_extra_dof = build_dof_map(frame)
    n = n_node_dof + n_extra_dof
    K = np.zeros((n, n))
    member_T = {}
    member_L = {}
    for mid, m in frame.members.items():
        section = frame.sections[m.section]
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
                  include_final_displacement=False, include_snapshots=False):
    """位移控制的遞增側推主迴圈。

    frame: 已定義節點/桿件/支承的Frame2D(側推的推力來自prescribed_dofs
        被強制位移, 不是外加point_load; 重力階段的point_load在
        apply_gravity()另外處理)
    hinge_states: {member_id: HingeState}, 通常從hinge.initial_hinge_states(frame)
        取得, 也可以是apply_gravity()跑完後(仍未降伏)的同一組物件——這個
        函式會直接原地修改傳進來的HingeState(設定yielded), 呼叫端要自己
        決定要不要事先複製一份
    prescribed_dofs: 受控自由度清單(全域dof編號, 例如[frame.dofs_of(roof_node)[0]])
    direction: 對應prescribed_dofs的方向係數(例如[1.0]表示單點控制, 或
        [1.0, 1.0]表示兩個控制點同向等量推)
    target_total: 目標總側推位移量
    d_nominal: 名目步長(沒有事件發生時每步走多少)
    base_reaction_dofs: 用來加總算底剪力的自由度(例如各支承節點的水平DOF)
    initial_cum_forces: apply_gravity()算出的初始元素內力, None代表從零
        開始(無重力預載)
    use_pdelta: 是否把frame元素目前軸力組進幾何剛度矩陣(P-Delta), 軸力
        代表值取該元素node_j端(拉力為正)
    mechanism_ratio_limit: 縮聚剛度最小/最大特徵值比例低於此值視為已達機構
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

    回傳: history_u(np.array), history_F(np.array), event_log(list of dict),
        hinge_states(原地更新後的同一組物件), mechanism_reached(bool)
        [, u_full_cum(np.array) -- 只有include_final_displacement=True時]
        [, history_snapshots(list) -- 只有include_snapshots=True時, 排在
          u_full_cum後面, 不管include_final_displacement是不是True]
    """
    direction = np.array(direction, dtype=float)
    fixed_dofs = _fixed_dof_set(frame)
    _, n_node_dof, n_extra_dof = build_dof_map(frame)
    n_total = n_node_dof + n_extra_dof

    if initial_cum_forces is None:
        cum_forces = {mid: np.zeros(6) for mid in frame.members}
    else:
        cum_forces = {mid: np.array(f, dtype=float) for mid, f in initial_cum_forces.items()}
    cum_reaction = np.zeros(n_total)

    lam = 0.0
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
            frame, hinge_states, axial_forces=axial)

        if check_mechanism(K, prescribed_dofs, fixed_dofs, mechanism_ratio_limit):
            mechanism_reached = True
            break

        du_full = solve_displacement_increment(K, prescribed_dofs, direction * d_step, fixed_dofs)
        df_by_member = _member_force_increments(frame, member_T, member_dofs, member_L,
                                                 hinge_states, axial, du_full)
        events = _find_crossing_events(hinge_states, cum_forces, df_by_member)

        if events:
            ratio_min = min(r for r, _, _ in events)
            d_sub = ratio_min * d_step
            du_full_sub = solve_displacement_increment(K, prescribed_dofs, direction * d_sub, fixed_dofs)
            df_by_member_sub = _member_force_increments(frame, member_T, member_dofs, member_L,
                                                          hinge_states, axial, du_full_sub)

            for mid, df in df_by_member_sub.items():
                cum_forces[mid] += df
            cum_reaction += K @ du_full_sub
            u_full_cum += du_full_sub
            lam += d_sub
            remaining -= d_sub

            newly_yielded = []
            for r, mid, end_idx in events:
                if abs(r - ratio_min) < 1e-9:
                    hinge_states[mid].yielded[end_idx] = True
                    newly_yielded.append((mid, end_idx))

            _update_theta_p(frame, hinge_states, member_T, member_dofs, member_L, u_full_cum)

            F_base = -sum(cum_reaction[d] for d in base_reaction_dofs)
            history_u.append(lam)
            history_F.append(F_base)
            event_log.append({'u': lam, 'F': F_base, 'yielded': newly_yielded})
            if include_snapshots:
                history_snapshots.append(_snapshot(cum_forces, hinge_states))
        else:
            for mid, df in df_by_member.items():
                cum_forces[mid] += df
            cum_reaction += K @ du_full
            u_full_cum += du_full
            lam += d_step
            remaining -= d_step

            _update_theta_p(frame, hinge_states, member_T, member_dofs, member_L, u_full_cum)

            F_base = -sum(cum_reaction[d] for d in base_reaction_dofs)
            history_u.append(lam)
            history_F.append(F_base)
            if include_snapshots:
                history_snapshots.append(_snapshot(cum_forces, hinge_states))

    result = [np.array(history_u), np.array(history_F), event_log, hinge_states, mechanism_reached]
    if include_final_displacement:
        result.append(u_full_cum)
    if include_snapshots:
        result.append(history_snapshots)
    return tuple(result)
