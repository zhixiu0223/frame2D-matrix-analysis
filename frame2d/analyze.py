"""
analyze_pushover() -- 統一入口, 純粹的dispatch層。

見ANALYSIS_ARCHITECTURE.md「規劃中」第2點的完整說明。這個檔案**不
包含任何新的物理或求解邏輯**——沒有新的力學公式、沒有新的疊代機制,
只是依material(由hinge_states本身決定, 不需要額外參數)/geometry=/
solver=這幾個獨立參數, 呼叫frame2d.pushover/frame2d.newton裡已經
各自驗證過的既有函式。

新增這一層**不會修改、也不會取代**任何一個現有函式——這是刻意的
安全網設計(對話紀錄裡達成的共識): 如果這個dispatch層本身有bug,
舊函式直接呼叫依然完全不受任何影響, 可以拿新舊兩條路徑的結果互相
對照(見tests/test_analyze_pushover.py, 每個組合都驗證過跟直接呼叫
底層函式逐位元/浮點數誤差為0一致)。
"""
from .pushover import run_pushover, run_pushover_converged
from .newton import run_pushover_newton, run_pushover_corotational_oneshot


VALID_GEOMETRY = ('linear', 'updated', 'corotational')
VALID_SOLVER = ('direct', 'picard', 'newton')


def analyze_pushover(frame, hinge_states, prescribed_dofs, direction, target_total,
                      d_nominal, base_reaction_dofs, control_mode='displacement',
                      geometry='linear', use_pdelta=False, use_local_pdelta=False,
                      solver='direct', **kwargs):
    """統一入口, 依geometry=/solver=兩個獨立參數dispatch到對應的既有
    函式——完整對照表見ANALYSIS_ARCHITECTURE.md。材料非線性(彈性 vs
    塑鉸)完全由hinge_states本身決定(空dict或None代表純彈性, 有內容
    代表塑鉸), 不需要額外的material=參數, 這件事本來就是現有函式的
    既有行為, 這裡沒有改變, 只是統一路由。

    geometry(幾何非線性怎麼處理):
      'linear'      -- 不管幾何非線性(不開P-Delta, 不開幾何更新)
      'updated'     -- 幾何更新(每步用當下變形後幾何重算L/角度, 仍是
                       小角度公式, 對應run_pushover()/run_pushover_
                       converged()的geometry_update=True)
      'corotational'-- 精確大轉角(co-rotational, 轉再多度也不失真,
                       對應run_pushover_newton()/run_pushover_
                       corotational_oneshot())

    use_pdelta: 只有geometry in ('linear','updated')時有意義——是否
      疊加global P-Delta(elements.local_geometric_stiffness()套在
      "絕對座標"v/θ上的線性化幾何剛度修正)。geometry='corotational'
      時這個參數會被忽略(不是報錯, 因為它對应的效應已經被co-rotational
      精確處理掉了, 用use_local_pdelta控制的是另一件獨立的事, 見下)。

    use_local_pdelta: 只有geometry='corotational'時有意義——是否疊加
      local P-Delta(同一條elements.local_geometric_stiffness()公式,
      但套用在co-rotational的自然座標θ_def上, 只處理"軸力對桿件自己
      局部彎曲勁度的修正", 跟大轉角本身是獨立的兩件事, 見
      ANALYSIS_ARCHITECTURE.md的說明)。geometry in ('linear','updated')
      時這個參數會被忽略。

    solver(怎麼把方程式解出來):
      'direct' -- 每一步只用"這一步開始時"的勁度矩陣解一次, 不疊代
                  (event-to-event式; geometry='corotational'時對應
                  run_pushover_corotational_oneshot(), 否則對應
                  run_pushover())
      'picard' -- 疊代到幾何/軸力自洽(不動點疊代, 只有geometry in
                  ('linear','updated')時支援, 對應run_pushover_
                  converged())
      'newton' -- 疊代到殘餘力收斂(Newton-Raphson, 只有geometry=
                  'corotational'時支援, 對應run_pushover_newton())

    不支援的geometry/solver組合會直接raise ValueError, 不是靜默退化
    成別的東西或忽略——例如geometry='corotational'配solver='picard'
    不存在(co-rotational的疊代天生就是Newton法, 沒有另外做一套Picard
    版本; 這不是"忘了實作", 是ANALYSIS_ARCHITECTURE.md裡明確記錄過
    的設計決定)。

    **kwargs: 原封不動轉傳給底下實際呼叫的函式(例如
    include_snapshots=True、include_final_displacement=True、tol=、
    max_iter=、geom_tol=、max_geom_iter=、initial_cum_forces=、
    mechanism_ratio_limit=等)——哪些kwargs有效取決於實際dispatch到
    哪個函式, 傳了該函式不認得的kwargs會直接反映成TypeError, 是刻意
    讓錯誤明顯冒出來, 不是靜默丟棄或吃掉錯誤。

    回傳: 直接回傳底下被呼叫函式的原始回傳值, 不做任何格式轉換或
    統一——不同geometry/solver組合的回傳形狀(欄位數量、順序)完全
    對應各自底層函式的簽名, 呼叫端要照底層函式的文件來解讀, 這裡
    刻意不做"統一包裝"以避免額外的轉換邏輯本身又成為新的bug來源。
    """
    if geometry not in VALID_GEOMETRY:
        raise ValueError(
            f"geometry必須是{VALID_GEOMETRY}其中之一, 實際收到'{geometry}'。"
        )
    if solver not in VALID_SOLVER:
        raise ValueError(
            f"solver必須是{VALID_SOLVER}其中之一, 實際收到'{solver}'。"
        )

    common = dict(
        frame=frame, hinge_states=hinge_states, prescribed_dofs=prescribed_dofs,
        direction=direction, target_total=target_total, d_nominal=d_nominal,
        base_reaction_dofs=base_reaction_dofs, control_mode=control_mode,
    )

    if geometry == 'corotational':
        if solver not in ('direct', 'newton'):
            raise ValueError(
                f"geometry='corotational'時solver只能是'direct'或'newton', "
                f"不支援'{solver}'——co-rotational的疊代天生就是Newton法, "
                f"沒有另外的Picard版本, 見ANALYSIS_ARCHITECTURE.md的說明。"
            )
        target_fn = run_pushover_corotational_oneshot if solver == 'direct' else run_pushover_newton
        return target_fn(**common, use_pdelta=use_local_pdelta, **kwargs)
    else:
        if solver not in ('direct', 'picard'):
            raise ValueError(
                f"geometry='{geometry}'時solver只能是'direct'或'picard', "
                f"不支援'{solver}'——Newton-Raphson目前只有搭配co-rotational"
                f"幾何的版本, 見ANALYSIS_ARCHITECTURE.md的說明。"
            )
        geometry_update = (geometry == 'updated')
        target_fn = run_pushover if solver == 'direct' else run_pushover_converged
        return target_fn(**common, use_pdelta=use_pdelta, geometry_update=geometry_update, **kwargs)
