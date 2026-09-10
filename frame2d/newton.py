"""
真正的Newton-Raphson遞增側推求解器, 用frame2d.corotational的大轉角
共旋公式——這是對話紀錄裡user明確要求的第三步(先確認要做疊代+
co-rotational, 這裡就是把兩者接起來變成完整的求解器)。

跟run_pushover()/run_pushover_converged()的關係: 三個函式完全獨立
(刻意不共用主迴圈, 見對話紀錄裡使用者/ChatGPT的建議), 適用範圍不同:
- run_pushover(): event-to-event, 材料非線性精確定位降伏事件, 幾何
  非線性用線性化P-Delta近似(或完全不用), 適合小轉角範圍的一般用途,
  速度最快。
- run_pushover_converged(): 在run_pushover()的event-to-event架構上,
  加一層"幾何/軸力疊代到自洽"——修正了"幾何更新用一次線性外推"這個
  問題, 但沒有修正"小角度假設本身在大轉角時失真"這個更根本的問題。
- run_pushover_newton()(這個檔案): 材料非線性(塑鉸降伏)跟幾何非線性
  (大轉角)都用真正的Newton-Raphson平衡疊代解, 用corotational.py的
  共旋公式(不是線性化P-Delta, 也不是"用變形後幾何硬套小角度公式"),
  轉角再大也不會失真——但相對的, 每一步都要花更多次疊代, 也是三者
  裡最慢的。

已知限制(誠實記錄):
- 不支援release端(有release_i/release_j的桿件, 直接raise清楚的
  錯誤訊息, 不是給錯誤答案)——release用額外不共用DOF體現這個既有
  機制, 目前沒有推廣到co-rotational版本, 這是刻意先縮小範圍的
  取捨, 不是隨便漏掉。
- 切線剛度矩陣用有限差分(見corotational.py開頭說明為什麼), 每次
  疊代要多算幾次內力函式, 比解析切線慢, 但保證跟殘餘力函式本身
  完全一致, 不會有"殘餘力對但切線推錯"這種難以察覺的錯誤。
- 材料降伏偵測是在Newton疊代"內部"做(每次疊代都檢查有沒有新的鉸
  超過Mp, 有的話標記降伏、標記完繼續疊代到收斂), 不是像event-to-
  event那樣先精確算出"降伏發生在這一步的哪個比例點"再切子步——這是
  一般非線性FEA軟體常見的做法(降伏偵測整合進疊代過程), 不是偷懶,
  但代表降伏發生的"確切位移點"不會像event-to-event那樣精確定位到
  很多位小數, 而是落在收斂到的那個位移點附近(差距量級是d_nominal
  這個步長, 步長切小差距就跟著變小)。
- 塑鉸彎矩用"參考狀態增量"追蹤(見_assemble_global()/corotational.
  corotational_local_forces()的說明), 修正了一個實際發生過的bug:
  塑鉸降伏瞬間R會從"近似剛接"驟降到R_post_yield, 如果直接拿總變形量
  套用新的(小很多的)R重新算彎矩, 算出來的M會遠低於Mp、瞬間"消失"
  大半, 這在物理上不合理(雙折線硬化模型, 降伏後M只會從Mp緩慢往上
  加, 不會倒退)——用"上一次成功收斂那一步"當參考點、只對"這一步
  多變形的量"套用目前的R算增量, 才能保持彎矩連續。
- 重力/桿件內部載重(distributed_loads、member_point_loads)有支援
  (見_gravity_fixed_end_forces()): 轉成固定端反力, 一部分變成等效
  節點力加進外力項, 一部分(固定端彎矩)從"變形產生的彎矩"裡扣掉,
  合起來才是塑鉸該檢查的"真正物理彎矩"——已經對照已驗證過的
  apply_gravity()跟wL^2/12公式驗證過精確一致。重力載重全程視為固定
  不變(用桿件原始幾何算一次, 不會隨pushover過程重新計算), 跟其餘
  求解器的apply_gravity()是同一種簡化, 不是另外發明規則。
"""
import numpy as np

from .corotational import corotational_kinematics, corotational_local_forces, \
    corotational_global_force, corotational_tangent_fd
from .elements import (
    member_geometry, transformation_matrix,
    fixed_end_forces_udl, fixed_end_forces_axial_udl_varying,
    fixed_end_forces_partial_udl, fixed_end_forces_axial_partial_udl,
    fixed_end_forces_point_load, fixed_end_forces_point_moment,
    fixed_end_forces_axial_point_load,
)


def _check_no_releases(frame):
    for mid, m in frame.members.items():
        if m.member_type == 'frame' and (m.release_i or m.release_j):
            raise ValueError(
                f"run_pushover_newton()目前不支援release端(桿件{mid}有"
                f"release_i={m.release_i}/release_j={m.release_j}設定)"
                f"——這是刻意先縮小範圍的取捨, 不是bug, 見本模組開頭"
                f"docstring的已知限制說明。"
            )


def _gravity_fixed_end_forces(frame):
    """把frame.distributed_loads/member_point_loads轉成"固定端反力"
    (跟dofmanager.py的_solve_once_dofmanager()同一套fixed_end_forces_*
    公式, 這裡是獨立重新寫一份呼叫, 不是共用那個函式本身——刻意不去
    動dofmanager.py, 避免任何風險影響到已經驗證過的線性求解路徑)。

    回傳(f_ext_gravity, member_fem_local):
        f_ext_gravity: n_dof長的全域"等效節點力"向量, 直接加進Newton
          求解的外力項(f_ext), 代表重力這類分布/桿件內部載重對整體
          平衡的貢獻。
        member_fem_local: {member_id: 6維局部座標固定端反力向量},
          用來修正_member_moments()算出來的"變形產生的彎矩", 讓塑鉸
          降伏判斷/顯示出來的彎矩是真正的物理彎矩(重力貢獻+變形貢獻
          兩者加總), 不是只有變形那一部分——公式跟dofmanager.py的
          end_forces_local = k_local@u_local - f_FE同一個慣例(見那邊
          第292-293行)。

    只支援frame.point_loads(直接是節點力, 不需要固定端反力處理)以外
    的distributed_loads、member_point_loads, 桁架/纜線桿件加了這些的
    話一樣明確報錯(理由同dofmanager.py的既有檢查, 這裡沿用同樣的
    錯誤訊息精神)。這裡用桿件"原始(未變形)幾何"算固定端反力+轉換
    矩陣, 重力預載全程視為固定不變, 不會隨著pushover過程重新計算
    ——跟其餘求解器的apply_gravity()是同一種簡化(重力預載定義一次,
    凍結不動), 不是另外發明新規則。
    """
    n_nodes = len(frame.nodes)
    n_dof = 3 * n_nodes
    f_ext_gravity = np.zeros(n_dof)
    member_fem_local = {mid: np.zeros(6) for mid in frame.members}
    member_T = {}
    member_L = {}
    for mid, m in frame.members.items():
        ni, nj = frame.nodes[m.node_i], frame.nodes[m.node_j]
        L, angle = member_geometry(ni, nj)
        member_T[mid] = transformation_matrix(angle)
        member_L[mid] = L

    def _add(mid, f_FE_local):
        member_fem_local[mid] = member_fem_local[mid] + f_FE_local
        dofs = list(frame.dofs_of(frame.members[mid].node_i)) + list(frame.dofs_of(frame.members[mid].node_j))
        f_ext_gravity[dofs] += member_T[mid].T @ f_FE_local

    for dl in frame.distributed_loads:
        m = frame.members[dl.member]
        if m.member_type in ('truss', 'cable'):
            raise ValueError(
                f"member {dl.member} 是{m.member_type}元素, 兩端鉸接、沒有彎曲"
                f"勁度, 不能承受分佈載重(理由同dofmanager.py既有檢查)。"
            )
        L = member_L[dl.member]
        if dl.direction == 'global_y':
            R = member_T[dl.member][0:2, 0:2]
            local_start = R @ np.array([0.0, -dl.w_start])
            local_end = R @ np.array([0.0, -dl.w_end])
            wx_start, wy_start = local_start
            wx_end, wy_end = local_end
            f_FE_local = (fixed_end_forces_udl(wy_start, wy_end, L)
                          + fixed_end_forces_axial_udl_varying(wx_start, wx_end, L))
        elif dl.direction == 'global':
            x_start = 0.0 if dl.x_start is None else dl.x_start
            x_end = L if dl.x_end is None else dl.x_end
            ang = np.radians(dl.angle_deg)
            u_global = np.array([np.cos(ang), np.sin(ang)])
            R = member_T[dl.member][0:2, 0:2]
            local_start = R @ (u_global * dl.w_start)
            local_end = R @ (u_global * dl.w_end)
            wx_start, wy_start = local_start
            wx_end, wy_end = local_end
            f_FE_local = (fixed_end_forces_partial_udl(wy_start, wy_end, x_start, x_end, L)
                          + fixed_end_forces_axial_partial_udl(wx_start, wx_end, x_start, x_end, L))
        else:
            x_start = 0.0 if dl.x_start is None else dl.x_start
            x_end = L if dl.x_end is None else dl.x_end
            if x_start <= 1e-9 and x_end >= L - 1e-9:
                f_FE_local = fixed_end_forces_udl(dl.w_start, dl.w_end, L)
            else:
                f_FE_local = fixed_end_forces_partial_udl(dl.w_start, dl.w_end, x_start, x_end, L)
        _add(dl.member, f_FE_local)

    for pl_m in frame.member_point_loads:
        m = frame.members[pl_m.member]
        if m.member_type in ('truss', 'cable'):
            raise ValueError(
                f"member {pl_m.member} 是{m.member_type}元素, 兩端鉸接、沒有彎曲"
                f"勁度, 不能承受桿件內部集中力/力矩(理由同dofmanager.py既有檢查)。"
            )
        L = member_L[pl_m.member]
        a = min(max(pl_m.a, 0.0), L)
        f_FE_local = np.zeros(6)
        if pl_m.direction == 'global':
            ang = np.radians(pl_m.angle_deg)
            u_global = np.array([np.cos(ang), np.sin(ang)])
            R = member_T[pl_m.member][0:2, 0:2]
            local_vec = R @ (u_global * pl_m.F)
            fx_local, fy_local = local_vec
        else:
            fx_local, fy_local = pl_m.fx, pl_m.fy
        if abs(fx_local) > 0:
            f_FE_local = f_FE_local + fixed_end_forces_axial_point_load(fx_local, a, L)
        if abs(fy_local) > 0:
            f_FE_local = f_FE_local + fixed_end_forces_point_load(fy_local, a, L)
        if abs(pl_m.m) > 0:
            f_FE_local = f_FE_local + fixed_end_forces_point_moment(pl_m.m, a, L)
        _add(pl_m.member, f_FE_local)

    return f_ext_gravity, member_fem_local


def _member_theta_def(frame, mid, u_full):
    """算某根桿件目前(u_full狀態下)兩端的變形轉角theta1_def/theta2_def
    (剛體轉動已扣掉)——事件偵測(有沒有哪一端剛好超過Mp)跟theta_p
    (累積塑性轉角)都要用這個, 不是原始的節點轉角。"""
    m = frame.members[mid]
    ni, nj = frame.nodes[m.node_i], frame.nodes[m.node_j]
    dofs = list(frame.dofs_of(m.node_i)) + list(frame.dofs_of(m.node_j))
    u_local = u_full[dofs]
    _, _, _, _, theta1_def, theta2_def, _ = corotational_kinematics(
        ni.x, ni.y, u_local[0], u_local[1], u_local[2],
        nj.x, nj.y, u_local[3], u_local[4], u_local[5])
    return theta1_def, theta2_def


def _member_moments(frame, mid, u_full, hinge_states, member_ref, member_fem_local=None):
    """算某根桿件目前(u_full狀態下)兩端的"真正物理彎矩"(M1,M2)——用
    member_ref(上一次成功收斂那一步的參考狀態)算變形產生的增量部分,
    見本檔案開頭docstring的說明。member_ref[mid]不存在時(這根桿件還
    沒有任何歷史, 例如第一次呼叫)視為(0,0,0,0)。

    member_fem_local: 可選, {member_id: 6維局部固定端反力向量}(來自
    _gravity_fixed_end_forces())——有給的話, 回傳的M1,M2會扣掉這根
    桿件的固定端彎矩分量(索引2,5), 讓回傳值是"變形貢獻+重力貢獻"
    加總後的真正物理彎矩, 塑鉸降伏判斷/顯示才會正確反映重力已經佔
    掉多少容量, 不是只看側推這一部分——公式跟dofmanager.py的
    end_forces_local=k_local@u_local-f_FE同一個慣例。不給(None)時
    完全等同沒有重力/桿件內部載重的情況, 維持原本行為。"""
    m = frame.members[mid]
    ni, nj = frame.nodes[m.node_i], frame.nodes[m.node_j]
    dofs = list(frame.dofs_of(m.node_i)) + list(frame.dofs_of(m.node_j))
    u_local = u_full[dofs]
    section = frame.sections[m.section]
    L0, L, beta, e, theta1_def, theta2_def, B = corotational_kinematics(
        ni.x, ni.y, u_local[0], u_local[1], u_local[2],
        nj.x, nj.y, u_local[3], u_local[4], u_local[5])
    hs = hinge_states.get(mid)
    M1_ref, M2_ref, th1_ref, th2_ref = member_ref.get(mid, (0.0, 0.0, 0.0, 0.0))
    N, M1, M2 = corotational_local_forces(
        section.E, section.A, section.I, L0, L, e, theta1_def, theta2_def, hs,
        M_ref=(M1_ref, M2_ref), theta_def_ref=(th1_ref, th2_ref))
    if member_fem_local is not None:
        fem = member_fem_local.get(mid)
        if fem is not None:
            M1 = M1 - fem[2]
            M2 = M2 - fem[5]
    return M1, M2


def _assemble_global(frame, hinge_states, u_full, n_dof, member_ref):
    """組出目前u_full狀態下, 全結構的內力向量(n_dof長)跟切線剛度矩陣
    (n_dof x n_dof), 疊加所有桿件(frame跟truss都用co-rotational公式,
    truss桿件沒有塑鉸/彎矩, 直接傳hinge_state=None且I用0讓彎矩項自然
    是0——這樣桁架桿件也能正確參與, 不用另外寫一套truss專屬邏輯)。

    member_ref: {member_id: (M1_ref,M2_ref,theta1_def_ref,theta2_def_ref)},
    每根桿件"上一次成功收斂那一步"的參考狀態, 塑鉸彎矩的路徑相依性
    (降伏後從Mp連續往上加, 不會倒退)完全靠這個維持。這一整個Newton
    "疊代"過程中(同一步內反覆試探), member_ref保持不變, 只有在一步
    真正收斂、被接受之後才會更新(見run_pushover_newton()主迴圈)。
    """
    f_int = np.zeros(n_dof)
    K_t = np.zeros((n_dof, n_dof))
    for mid, m in frame.members.items():
        ni, nj = frame.nodes[m.node_i], frame.nodes[m.node_j]
        dofs = list(frame.dofs_of(m.node_i)) + list(frame.dofs_of(m.node_j))
        u_local = u_full[dofs]
        section = frame.sections[m.section]
        E, A, I = section.E, section.A, section.I
        if m.member_type == 'truss':
            I = 0.0
        hs = hinge_states.get(mid) if (hinge_states is not None and m.member_type == 'frame') else None
        M1_ref, M2_ref, th1_ref, th2_ref = member_ref.get(mid, (0.0, 0.0, 0.0, 0.0))
        M_ref = (M1_ref, M2_ref)
        theta_def_ref = (th1_ref, th2_ref)

        f_elem = corotational_global_force(
            ni.x, ni.y, u_local[0], u_local[1], u_local[2],
            nj.x, nj.y, u_local[3], u_local[4], u_local[5], E, A, I, hs, M_ref, theta_def_ref)
        K_elem = corotational_tangent_fd(
            ni.x, ni.y, u_local[0], u_local[1], u_local[2],
            nj.x, nj.y, u_local[3], u_local[4], u_local[5], E, A, I, hs, M_ref, theta_def_ref)

        for a in range(6):
            f_int[dofs[a]] += f_elem[a]
            for b in range(6):
                K_t[dofs[a], dofs[b]] += K_elem[a, b]
    return f_int, K_t


def run_pushover_newton(frame, hinge_states, prescribed_dofs, direction, target_total,
                         d_nominal, base_reaction_dofs,
                         control_mode='displacement', tol=1e-6, max_iter=30, max_steps=100000,
                         include_final_displacement=False, include_snapshots=False):
    """真正的Newton-Raphson遞增側推。見本檔案開頭docstring。

    frame, hinge_states, prescribed_dofs, direction, target_total, d_nominal,
    base_reaction_dofs, control_mode: 跟run_pushover()的同名參數意義完全
    一致(位移控制或力控制、多點加權等)。

    tol: 殘餘力收斂容忍度——跟目前典型內力量級的相對比值;
    max_iter: 每一步最多疊代幾次。

    回傳: history_u, history_F, event_log, hinge_states, converged(bool,
    True代表整個側推過程每一步都成功收斂; False代表在某一步疊代
    max_iter次還沒收斂, 已提前停止, 後面的位移/力都不會出現在歷程裡)
    [, u_full_cum(np.array) -- 只有include_final_displacement=True時]
    [, history_snapshots(list) -- 只有include_snapshots=True時, 排在
      u_full_cum後面, 不管include_final_displacement是不是True]

    event_log的每筆紀錄跟run_pushover()格式一致: {'u':.., 'F':.., 
    'yielded': [(member_id, end_idx), ...]}——但這裡的u是"疊代收斂到
    的那一步"的位移, 不是event-to-event那種精確定位到降伏當下比例點
    的位移(見本檔案開頭"已知限制"說明)。
    """
    _check_no_releases(frame)
    direction = np.array(direction, dtype=float)
    n_nodes = len(frame.nodes)
    n_dof = 3 * n_nodes
    fixed_dofs = set()
    for s in frame.supports:
        ux_i, uy_i, rot_i = frame.dofs_of(s.node)
        if s.ux is not None:
            fixed_dofs.add(ux_i)
        if s.uy is not None:
            fixed_dofs.add(uy_i)
        if s.rot is not None:
            fixed_dofs.add(rot_i)
    free_dofs = [d for d in range(n_dof) if d not in fixed_dofs]
    control_dof0 = prescribed_dofs[0]
    if not free_dofs:
        raise ValueError(
            "這個模型所有自由度都被支承條件固定住了, 沒有任何自由度"
            "可以求解(常見情況: 控制節點本身也被設成完全固定支承,"
            "跟pushover想推動它互相矛盾), 請檢查支承條件。"
        )

    # 重力/桿件內部載重(distributed_loads、member_point_loads)轉成
    # 等效節點力(f_ext_gravity)+每根桿件的局部固定端反力
    # (member_fem_local)——全程視為固定不變(跟其餘求解器apply_gravity()
    # 同一種簡化, 不是隨pushover過程重算), 見_gravity_fixed_end_forces()
    # 說明。f_int(Newton疊代用的內力函式)只算"變形產生的部分", 重力
    # 完全走f_ext這條路; 塑鉸降伏判斷/顯示用的"真正物理彎矩"則要另外
    # 扣掉固定端彎矩才對(見_member_moments()的member_fem_local參數)。
    f_ext_gravity, member_fem_local = _gravity_fixed_end_forces(frame)

    u_full = np.zeros(n_dof)
    f_ext_cum = f_ext_gravity.copy()   # 力控制模式下, 累積施加的外力
                                        # (從重力這個固定基礎開始疊加
                                        # pushover本身的力, 不是每一步
                                        # 重算, 因為u_full是絕對量)
    # 每根桿件"上一次成功收斂那一步"的(M1,M2,theta1_def,theta2_def)
    # 參考狀態, 見本檔案開頭docstring/_assemble_global()的說明——這裡
    # 存的永遠是"變形產生的部分"(不含重力固定端彎矩), 因為
    # corotational_local_forces()的增量公式就是對這部分操作的; 重力
    # 貢獻要另外用member_fem_local在需要"真正物理彎矩"的地方(降伏
    # 判斷、顯示)另外扣回去, 不能混進member_ref裡, 否則會被
    # corotational_local_forces()的增量邏輯重複計算。一開始(u=0)全部
    # 是零。
    member_ref = {mid: (0.0, 0.0, 0.0, 0.0) for mid in frame.members}
    # 每個塑鉸端"第一次降伏當下"的theta_def, 用來算theta_p(累積塑性
    # 轉角)——這是"固定"參考點(只在剛降伏那一刻設定一次, 之後不會
    # 再變), 跟member_ref(每步都更新)是兩個不同用途的參考點, 不要
    # 混用。
    theta_def_at_yield = {mid: [None, None] for mid in hinge_states}

    history_u = [0.0]
    history_F = [0.0]
    event_log = []
    converged_all = True
    history_snapshots = None
    if include_snapshots:
        from .pushover import _snapshot
        cum_forces_display = {}
        for mid in frame.members:
            fem = member_fem_local[mid]
            cum_forces_display[mid] = np.array([0.0, 0.0, -fem[2], 0.0, 0.0, -fem[5]])
        history_snapshots = [_snapshot(cum_forces_display, hinge_states, u_full)]

    remaining = target_total
    steps = 0
    while remaining > 1e-9:
        steps += 1
        if steps > max_steps:
            raise RuntimeError(
                f"側推超過{max_steps}步仍未達到target_total, 可能是d_nominal"
                "設太小或有其他問題, 已中止(避免無窮迴圈)。")
        d_step = min(d_nominal, remaining)

        u_trial = u_full.copy()
        if control_mode == 'force':
            for i, dof in enumerate(prescribed_dofs):
                f_ext_cum[dof] += direction[i] * d_step
            f_ext = f_ext_cum
        else:
            for i, dof in enumerate(prescribed_dofs):
                u_trial[dof] += direction[i] * d_step

        step_converged = False
        step_newly_yielded = []
        for it in range(max_iter):
            f_int, K_t = _assemble_global(frame, hinge_states, u_trial, n_dof, member_ref)

            # 材料降伏偵測: 用"真正物理彎矩"(變形貢獻+重力固定端彎矩
            # 貢獻, 見_member_moments()的member_fem_local參數)檢查有
            # 沒有哪個還沒降伏的端超過Mp——這裡在疊代"內部"做, 見本
            # 檔案開頭"已知限制"說明。newly_yielded_this_iter要累加進
            # step_newly_yielded(不能只看最後一次疊代的結果, 同一步
            # 裡不同疊代次可能先後偵測到不同的鉸降伏)。
            newly_yielded_this_iter = []
            for mid, hs in hinge_states.items():
                M1, M2 = _member_moments(frame, mid, u_trial, hinge_states, member_ref, member_fem_local)
                if not hs.yielded[0] and abs(M1) >= hs.Mp[0]:
                    hs.yielded[0] = True
                    newly_yielded_this_iter.append((mid, 0))
                    th1d, _ = _member_theta_def(frame, mid, u_trial)
                    theta_def_at_yield[mid][0] = th1d
                if not hs.yielded[1] and abs(M2) >= hs.Mp[1]:
                    hs.yielded[1] = True
                    newly_yielded_this_iter.append((mid, 1))
                    _, th2d = _member_theta_def(frame, mid, u_trial)
                    theta_def_at_yield[mid][1] = th2d
            if newly_yielded_this_iter:
                step_newly_yielded.extend(newly_yielded_this_iter)
                continue   # 塑鉸狀態變了, 這次試探已經過期, 重新組裝再試一次

            if control_mode == 'force':
                resid_dofs = [d for d in free_dofs]
                residual = (f_ext - f_int)[resid_dofs]
            else:
                resid_dofs = [d for d in free_dofs if d not in prescribed_dofs]
                residual = (f_ext_gravity - f_int)[resid_dofs]

            scale = max(np.max(np.abs(f_int)), 1.0)
            if np.max(np.abs(residual)) < tol * scale:
                step_converged = True
                break

            K_sub = K_t[np.ix_(resid_dofs, resid_dofs)]
            try:
                du_sub = np.linalg.solve(K_sub, residual)
            except np.linalg.LinAlgError:
                break   # 切線奇異(通常代表已經到極限承載力附近), 視為不收斂
            for i, d in enumerate(resid_dofs):
                u_trial[d] += du_sub[i]

        if not step_converged:
            converged_all = False
            break

        u_full = u_trial
        # 這一步真正收斂了, 把每根桿件目前"變形產生的部分"
        # (M1,M2,theta1_def,theta2_def)存成新的member_ref, 供下一步
        # (以及theta_p計算)使用——這是"增量參考點只在步驟被接受時才
        # 更新"這個設計的關鍵一步, 塑鉸彎矩的連續性完全靠這裡維持。
        # 這裡刻意不傳member_fem_local(要存"變形部分"本身, 不是真正
        # 物理彎矩, 否則下一步corotational_local_forces()的增量公式
        # 會把重力貢獻重複計算進去)。
        f_int_final, _ = _assemble_global(frame, hinge_states, u_full, n_dof, member_ref)
        for mid in frame.members:
            th1d, th2d = _member_theta_def(frame, mid, u_full)
            M1_def, M2_def = _member_moments(frame, mid, u_full, hinge_states, member_ref)
            member_ref[mid] = (M1_def, M2_def, th1d, th2d)

        # 已經用懸臂樑解析解驗證過, f_int在固定支承dof上的值取負號後,
        # 才是"底剪力"這個物理量本身(標準的cum_reaction/K@u慣例)。
        F_base = -sum(f_int_final[d] for d in base_reaction_dofs)
        u_control = u_full[control_dof0]
        history_u.append(u_control)
        history_F.append(F_base)
        remaining -= d_step

        for mid, hs in hinge_states.items():
            th1d, th2d = _member_theta_def(frame, mid, u_full)
            for end_idx, th_def in [(0, th1d), (1, th2d)]:
                if hs.yielded[end_idx] and theta_def_at_yield[mid][end_idx] is not None:
                    hs.theta_p[end_idx] = max(
                        0.0, abs(th_def) - abs(theta_def_at_yield[mid][end_idx]))

        if include_snapshots:
            from .pushover import _snapshot
            cum_forces_display = {}
            for mid in frame.members:
                M1, M2 = _member_moments(frame, mid, u_full, hinge_states, member_ref, member_fem_local)
                cum_forces_display[mid] = np.array([0.0, 0.0, M1, 0.0, 0.0, M2])
            history_snapshots.append(_snapshot(cum_forces_display, hinge_states, u_full))

        if step_newly_yielded:
            event_log.append({'u': u_control, 'F': F_base, 'yielded': step_newly_yielded})

    result = [np.array(history_u), np.array(history_F), event_log, hinge_states, converged_all]
    if include_final_displacement:
        result.append(u_full)
    if include_snapshots:
        result.append(history_snapshots)
    return tuple(result)
