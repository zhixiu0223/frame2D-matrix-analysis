"""
PDF 匯出: 直接重用 frame2d/plotting.py 已經驗證過的 plot_all()
(六合一總覽圖: 結構/受力/變形/N圖/V圖/M圖),不在這層重新畫任何圖,
避免跟 plotting.py 的正負號慣例(拉力側等)兩邊分別維護、互相漂移。

matplotlib 用 Agg(非互動式)backend,在沒有顯示器的伺服器環境
(Termux/Colab/Cloud Run)都能正常產生圖檔。

第二部分: 完整輸入資料表格頁(節點/斷面/桿件/支承/載重)+反力/位移
結果表格頁——之前的PDF只有六合一圖表, 看報告的人沒辦法從圖上目測
出精確的座標/斷面性質/載重數值, 沒辦法拿這份報告去手算驗證或在
別的軟體重現同一個模型。這裡補上完整的表格頁, 用matplotlib的
ax.table()畫, 跟結構圖用同一套PdfPages合併成一份多頁PDF。
"""
import base64
import io

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from frame2d import solve
from frame2d.plotting import plot_all, plot_member_fbd, plot_member_own_diagrams, plot_structure


def _fig_to_png_base64(fig, dpi=110):
    """把matplotlib figure轉成base64編碼的PNG字串, 用在「匯出前先在
    瀏覽器預覽」這個功能——不用另外存檔案再回傳URL, 直接把圖片編碼
    塞進JSON回應裡, 前端<img src="data:image/png;base64,...">就能
    直接顯示, 不用額外的靜態檔案伺服流程。"""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _force_moment_factors(units):
    """把units字典換算成(force_unit, force_factor, moment_unit, moment_factor)
    這組四元組, 直接餵給frame2d/plotting.py的plot_member_fbd()/
    plot_member_own_diagrams()/plot_all()等函式——那些函式本身不認識
    "kN"這種單位名稱, 只認係數+要印出來的字串, 這裡負責從webapi層的
    units設定轉換成它們要的格式。"""
    fu = _unit_label(units, "force", "N")
    mu = _unit_label(units, "moment", "N·m")
    return fu, 1.0 / UNIT_FACTORS["force"][fu], mu, 1.0 / UNIT_FACTORS["moment"][mu]


def _disp_factor(units):
    """回傳(disp_unit, disp_factor), 給plot_all()/plot_deformed()的
    位移標籤用。"""
    du = _unit_label(units, "disp", "m")
    return du, 1.0 / UNIT_FACTORS["disp"][du]


def _member_fbd_with_thumbnail_fig(f, result, mid, force_unit, force_factor, moment_unit, moment_factor, figsize=(12, 7)):
    """自由體圖+旁邊一張結構縮圖(目前這根桿件用橘紅色標示)合併成
    一張figure——跟瀏覽器預覽(build_fbd_previews, 分開兩張圖用CSS
    排版)不一樣, PDF裡沒有CSS可以排版, 直接用matplotlib的GridSpec
    把兩個subplot擺在同一張圖裡, 縮圖佔窄的1/4, 自由體圖佔剩下3/4。"""
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 4)
    thumb_ax = fig.add_subplot(gs[0, 0])
    plot_structure(f, ax=thumb_ax, show_node_ids=False, show_member_ids=False,
                    show_dimensions=False, highlight_member_id=mid)
    thumb_ax.set_title(f'Member {mid}', fontsize=9)
    thumb_ax.set_xlabel('')
    thumb_ax.set_ylabel('')
    fbd_ax = fig.add_subplot(gs[0, 1:])
    plot_member_fbd(f, result, mid, ax=fbd_ax, force_unit=force_unit, force_factor=force_factor,
                     moment_unit=moment_unit, moment_factor=moment_factor)
    fig.tight_layout()
    return fig


def build_fbd_previews(f, member_ids, units=None):
    """回傳每根指定桿件的自由體圖預覽(base64 PNG), 給前端在真正
    匯出PDF之前先秀出來讓使用者確認、可以個別排除不想要的桿件。
    回傳list of {"member_id":int, "image_base64":str, "thumbnail_base64":str},
    找不到的member_id直接跳過(不報錯, 讓前端自己比對哪些真的有回傳)。

    units: 前端顯示單位設定, None用SI——之前這裡沒有接units, 圖上
    的Fx/Fy/M數字永遠是SI(N/N*m), 跟使用者在「顯示設定」選的單位
    (例如kN)對不起來, 這次修正。

    thumbnail_base64: 整個結構的小縮圖, 目前這根桿件用醒目橘紅色
    粗線標示、其他桿件淡化——目的是選「全部桿件」批次預覽時, 每張
    自由體圖旁邊都能一眼看出對應到結構的哪個位置, 不用只憑桿件
    編號自己去對照。"""
    fu, ff, mu, mf = _force_moment_factors(units)
    result = solve(f)
    previews = []
    for mid in member_ids:
        if mid not in f.members:
            continue
        ax = plot_member_fbd(f, result, mid, force_unit=fu, force_factor=ff, moment_unit=mu, moment_factor=mf)
        image_b64 = _fig_to_png_base64(ax.figure)
        plt.close(ax.figure)

        thumb_fig, thumb_ax = plt.subplots(figsize=(3, 3))
        plot_structure(f, ax=thumb_ax, show_node_ids=False, show_member_ids=False,
                        show_dimensions=False, highlight_member_id=mid)
        thumb_ax.set_title(f'Member {mid}', fontsize=9)
        thumb_ax.set_xlabel(''); thumb_ax.set_ylabel('')
        thumb_fig.tight_layout()
        thumb_b64 = _fig_to_png_base64(thumb_fig, dpi=90)
        plt.close(thumb_fig)

        previews.append({"member_id": mid, "image_base64": image_b64, "thumbnail_base64": thumb_b64})
    return previews


def _table_page(title, blocks, figsize=(14, 9)):
    """blocks: [(subtitle, column_labels, rows(list of list, already stringified)), ...]
    每個block各自畫一張表, 由上到下排在同一頁, 高度依資料筆數自動分配。
    回傳一個matplotlib Figure。全部用英文, 避開matplotlib預設字型
    (DejaVu Sans)不支援中文字元的問題(這份報告要能在Termux/Colab/
    Cloud Run等不一定有裝中文字型的環境正常產生)。"""
    fig = plt.figure(figsize=figsize)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    n = len(blocks)
    # 每個block的高度依row數決定, 至少留1個row高度給空表格的情況
    row_counts = [max(len(rows), 1) for _, _, rows in blocks]
    total = sum(row_counts) + 2 * n  # 標題+表頭抓2個row高度的餘裕
    y = 0.94
    avail = 0.88
    for (subtitle, cols, rows), rc in zip(blocks, row_counts):
        h = avail * (rc + 2) / total
        ax = fig.add_axes([0.04, y - h, 0.92, h * 0.86])
        ax.axis("off")
        ax.set_title(subtitle, fontsize=10, fontweight="bold", loc="left")
        if rows:
            tbl = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center")
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8)
            tbl.scale(1, 1.3)
        else:
            ax.text(0.02, 0.5, "(none)", fontsize=9, color="#888888")
        y -= h
    return fig


def _fmt(v, nd=6):
    if v is None:
        return "free"
    if isinstance(v, float):
        return f"{v:.{nd}g}"
    return str(v)


# 跟 webapi/static/index.html 的 UNIT_FACTORS 完全對應的Python版本
# ——PDF匯出時如果前端有帶units參數過來, 用這裡的係數把SI值換算成
# 使用者畫面上實際看到的單位, 讓PDF報告的數字/單位標籤跟使用者當下
# 的顯示設定一致, 不用另外猜這份報告到底是哪個單位基準。
UNIT_FACTORS = {
    "E": {"Pa": 1, "MPa": 1e6, "GPa": 1e9},
    "I": {"m4": 1, "mm4": 1e-12, "cm4": 1e-8},
    "A": {"m2": 1, "mm2": 1e-6, "cm2": 1e-4},
    "disp": {"m": 1, "cm": 1e-2, "mm": 1e-3},
    "force": {"N": 1, "kN": 1e3},
    "moment": {"N·m": 1, "kN·m": 1e3},
    # 動力分析單位(基準: kg, kg·m², kg/m³, m/s², m/s), 對應index.html新增的選單
    "mass": {"kg": 1, "t": 1e3},
    "inertia": {"kg·m²": 1, "t·m²": 1e3},
    "density": {"kg/m³": 1, "t/m³": 1e3},
    "accel": {"m/s²": 1, "g": 9.80665, "gal": 1e-2},
    "vel": {"m/s": 1, "cm/s": 1e-2, "mm/s": 1e-3},
}


def _unit_label(units, kind, default):
    """units是前端傳過來的{"E":"GPa",...}字典(可能是None或缺某個key),
    回傳實際要用的單位字串; 沒有指定的話用SI(default)。"""
    if not units or kind not in units or units[kind] not in UNIT_FACTORS[kind]:
        return default
    return units[kind]


def _from_si(units, kind, si_val, default_unit):
    """把SI值換算成units指定的顯示單位(跟前端fromSI()同一套公式);
    units沒指定該kind時, 直接用SI(等於除以1, 不變)。"""
    unit = _unit_label(units, kind, default_unit)
    return si_val / UNIT_FACTORS[kind][unit]


def build_input_data_pages(f, units=None):
    """把Frame2D f 的所有輸入資料(節點/斷面/桿件/支承/三種載重)排成
    表格頁, 回傳figure的list——目的是讓報告本身就能完全重現模型或
    手算驗證, 不用只能從六合一圖上目測。表格標籤(欄位名稱)全部用
    英文, 跟plot_all()既有的慣例一致——這不是隨便選的, 是刻意避開
    matplotlib預設字型(DejaVu Sans)不支援中文字元的問題: 這份報告
    要能在Termux/Colab/Cloud Run等沒有裝中文字型的環境正常產生,
    不能假設部署環境一定有CJK字型可用。

    units: 前端「顯示設定」目前選的單位({"E":"GPa","force":"kN",...}),
    None的話用SI——數值跟欄位標題裡的單位字串會一起換算, 讓報告
    上看到的數字/單位標籤永遠一致(不會發生"標題寫kN, 數字卻是N"
    這種對不起來的情況)。"""
    eu = _unit_label(units, "E", "Pa")
    iu = _unit_label(units, "I", "m4")
    au = _unit_label(units, "A", "m2")
    fu = _unit_label(units, "force", "N")
    mu = _unit_label(units, "moment", "N·m")

    nodes = sorted(f.nodes.values(), key=lambda n: n.id)
    sections = sorted(f.sections.values(), key=lambda s: s.name)
    members = sorted(f.members.values(), key=lambda m: m.id)
    supports = sorted(f.supports, key=lambda s: s.node)
    point_loads = sorted(f.point_loads, key=lambda p: p.node)
    dist_loads = sorted(f.distributed_loads, key=lambda d: d.member)
    member_point_loads = sorted(f.member_point_loads, key=lambda p: p.member)

    node_rows = [[str(n.id), _fmt(n.x), _fmt(n.y)] for n in nodes]  # 座標固定用公尺(跟前端一致, 目前沒有座標單位設定)
    section_rows = [[
        s.name, _fmt(_from_si(units, "E", s.E, "Pa")),
        _fmt(_from_si(units, "I", s.I, "m4")), _fmt(_from_si(units, "A", s.A, "m2")),
    ] for s in sections]
    page1 = _table_page(f"Input Data (1/3): Node Coordinates / Section Properties", [
        ("Node Coordinates (m)", ["Node ID", "x (m)", "y (m)"], node_rows),
        (f"Section Properties (E in {eu}, I in {iu}, A in {au})",
         ["Section Name", f"E ({eu})", f"I ({iu})", f"A ({au})"], section_rows),
    ])

    member_rows = [[
        str(m.id), str(m.node_i), str(m.node_j), m.section, m.member_type,
        "yes" if m.release_i else "no", "yes" if m.release_j else "no",
        _fmt(_from_si(units, "moment", m.Mp_i, "N·m")) if m.Mp_i is not None else "-",
        _fmt(_from_si(units, "moment", m.R_post_yield_i, "N·m")) if m.R_post_yield_i is not None else "-",
        _fmt(_from_si(units, "moment", m.Mp_j, "N·m")) if m.Mp_j is not None else "-",
        _fmt(_from_si(units, "moment", m.R_post_yield_j, "N·m")) if m.R_post_yield_j is not None else "-",
    ] for m in members]
    support_rows = [[str(s.node), _fmt(s.ux), _fmt(s.uy), _fmt(s.rot)] for s in supports]
    page2 = _table_page("Input Data (2/3): Member Connectivity / Support Conditions", [
        (f"Member Connectivity (plastic hinge capacity: Mp/R in {mu}, blank = never yields -- see Pushover section)",
         ["Member ID", "Node i", "Node j", "Section", "Type", "Release i", "Release j",
          f"Mp_i ({mu})", f"R_i ({mu}/rad)", f"Mp_j ({mu})", f"R_j ({mu}/rad)"], member_rows),
        ("Support Conditions (m, rad; blank = free, 0 = fixed, nonzero = prescribed displacement)",
         ["Node ID", "ux (m)", "uy (m)", "rot (rad)"], support_rows),
    ])

    pl_rows = [[
        str(p.node), _fmt(_from_si(units, "force", p.fx, "N")),
        _fmt(_from_si(units, "force", p.fy, "N")), _fmt(_from_si(units, "moment", p.m, "N·m")),
    ] for p in point_loads]
    dl_rows = [[
        str(d.member), _fmt(_from_si(units, "force", d.w_start, "N")),
        _fmt(_from_si(units, "force", d.w_end, "N")) if d.w_end is not None else "= w1",
        _fmt(d.x_start) if d.x_start is not None else "0",
        _fmt(d.x_end) if d.x_end is not None else "full length",
        d.direction, _fmt(d.angle_deg) if d.angle_deg is not None else "-",
    ] for d in dist_loads]
    mpl_rows = [[
        str(p.member), _fmt(p.a), p.direction,
        _fmt(_from_si(units, "force", p.fx, "N")) if p.direction == "local" else "-",
        _fmt(_from_si(units, "force", p.fy, "N")) if p.direction == "local" else "-",
        _fmt(_from_si(units, "force", p.F, "N")) if p.direction == "global" else "-",
        _fmt(p.angle_deg) if p.direction == "global" else "-",
        _fmt(_from_si(units, "moment", p.m, "N·m")),
    ] for p in member_point_loads]
    page3 = _table_page("Input Data (3/3): Nodal Point Loads / Distributed Loads / Member Point Loads", [
        (f"Nodal Point Loads (fx/fy in {fu}, m in {mu})", ["Node ID", f"fx ({fu})", f"fy ({fu})", f"m ({mu})"], pl_rows),
        (f"Distributed Loads (w in {fu}/m, measured along member length)",
         ["Member ID", "w1 (start)", "w2 (end)", "x start (m)", "x end (m)", "direction", "angle (deg)"], dl_rows),
        (f"Member Internal Point Loads (a = distance from node i along member axis, in m; forces in {fu}, moment in {mu})",
         ["Member ID", "a (m)", "direction", f"fx ({fu})", f"fy ({fu})", f"F ({fu})", "angle (deg)", f"m ({mu})"], mpl_rows),
    ])
    return [page1, page2, page3]


def build_mass_data_page(f, units=None):
    """質量設定(斷面密度ρ + 節點集中質量)表格頁, 給D2以後的動力分析報告用(模態/反應譜/
    循環/非線性地震)——D0/D1(線性/P-Delta/Pushover)不需要質量矩陣, 沒有這一頁, 這裡只在
    build_modal_pdf_report()等四個新函式裡呼叫。跟 build_input_data_pages() 用同一套
    _table_page()/_fmt()/_from_si()機制, 保持格式一致。沒有任何ρ或節點質量時仍然回傳一頁
    (兩個表都是空的), 讓報告頁碼固定、呼叫端不用另外判斷要不要插這一頁。"""
    du_mass = _unit_label(units, "mass", "kg")
    du_inertia = _unit_label(units, "inertia", "kg·m²")
    du_density = _unit_label(units, "density", "kg/m³")

    sections = sorted(f.sections.values(), key=lambda s: s.name)
    rho_rows = [[s.name, _fmt(_from_si(units, "density", s.rho, "kg/m³"))] for s in sections if s.rho]
    nodes_with_mass = sorted(f.node_masses, key=lambda nm: nm.node)
    mass_rows = [[
        str(nm.node),
        _fmt(_from_si(units, "mass", nm.mx, "kg")) if nm.mx else "-",
        _fmt(_from_si(units, "mass", nm.my, "kg")) if nm.my else "-",
        _fmt(_from_si(units, "inertia", nm.Iz, "kg·m²")) if nm.Iz else "-",
    ] for nm in nodes_with_mass]
    return _table_page(f"Input Data: Mass (dynamic analysis)", [
        (f"Section Density (ρ in {du_density})", ["Section Name", f"ρ ({du_density})"], rho_rows),
        (f"Nodal Lumped Mass (translational in {du_mass}, rotational inertia in {du_inertia})",
         ["Node ID", f"mx ({du_mass})", f"my ({du_mass})", f"Iz ({du_inertia})"], mass_rows),
    ])


def build_result_data_page(f, result, units=None):
    """反力/位移結果表格頁, 讓報告本身就能對照手算或別的軟體算出的
    反力/位移數值, 不用只能從圖上目測峰值。units同build_input_data_pages。"""
    du = _unit_label(units, "disp", "m")
    fu = _unit_label(units, "force", "N")
    mu = _unit_label(units, "moment", "N·m")
    node_rows = []
    for n in sorted(f.nodes.values(), key=lambda n: n.id):
        ux, uy, rot = f.dofs_of(n.id)
        node_rows.append([
            str(n.id),
            _fmt(_from_si(units, "disp", result.displacements[ux], "m")),
            _fmt(_from_si(units, "disp", result.displacements[uy], "m")),
            _fmt(result.displacements[rot]),
            _fmt(_from_si(units, "force", result.reactions[ux], "N")),
            _fmt(_from_si(units, "force", result.reactions[uy], "N")),
            _fmt(_from_si(units, "moment", result.reactions[rot], "N·m")),
        ])
    return _table_page(f"Solve Results: Node Displacements / Reactions", [
        (f"Node Displacements ({du}, rad) & Reactions ({fu}, {mu}; reaction is nonzero only where that DOF is restrained)",
         ["Node ID", f"ux ({du})", f"uy ({du})", "rot (rad)", f"Rx ({fu})", f"Ry ({fu})", f"M ({mu})"], node_rows),
    ])


def build_fbd_images_archive(f, member_ids, units=None) -> bytes:
    """回傳一個zip檔的bytes, 裡面是每根指定桿件的自由體圖PNG(跟
    PDF裡的版面一樣: 自由體圖+旁邊的結構縮圖合併成一張圖, 不是
    分開兩個檔案)。單一桿件時, 呼叫端(main.py/server.py)可以決定
    直接回傳裡面那張PNG而不包zip——這裡固定產生zip, 由呼叫端依
    member_ids長度決定要不要解開, 邏輯不用寫兩次。"""
    import zipfile
    fu, ff, mu, mf = _force_moment_factors(units)
    result = solve(f)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for mid in member_ids:
            if mid not in f.members:
                continue
            fig = _member_fbd_with_thumbnail_fig(f, result, mid, fu, ff, mu, mf)
            img_buf = io.BytesIO()
            fig.savefig(img_buf, format="png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            zf.writestr(f"member_{mid}_fbd.png", img_buf.getvalue())
    return buf.getvalue()


def _pushover_final_solve_result(f, cum_forces, u_full_cum, cum_reaction):
    """從Pushover最後一步的累積內力/位移/反力, 組出一個標準的SolveResult,
    直接餵給plot_all()/plot_diagram()/plot_deformed()——這幾個函式本來
    就是拿(frame, result)當參數, 不管result是"線性solve()算出來的"還是
    "Pushover最後一步的狀態", 只要end_forces_local跟displacements的
    符號慣例一致就能正常運作。已經確認過cum_forces跟MemberResult.
    end_forces_local是同一種原始慣例(都是k_local@u_local直接算出來的,
    沒有另外做符號校正), 不需要額外轉換。

    L/angle用未變形的原始幾何(跟frame.members當初定義的一樣), 不是
    Pushover過程中(如果開了geometry_update)每一步重算的變形後幾何——
    這裡刻意畫的是"最終那一步的內力/位移狀態長什麼樣子", 用原始幾何
    畫圖比較符合"疊在結構原始外形上看內力分布"這個報告用途的直覺,
    不是要重現geometry_update內部計算時用的每一步幾何。
    """
    from frame2d.result import SolveResult, MemberResult
    from frame2d.elements import member_geometry
    member_results = {}
    for mid, m in f.members.items():
        ni, nj = f.nodes[m.node_i], f.nodes[m.node_j]
        L, angle = member_geometry(ni, nj)
        member_results[mid] = MemberResult(
            member_id=mid, L=L, angle=angle,
            end_forces_local=np.array(cum_forces[mid], dtype=float), slack=False,
        )
    return SolveResult(displacements=u_full_cum, reactions=cum_reaction,
                        member_results=member_results, frame=f)


def build_pushover_capacity_curve_page(history_u, history_F, event_log, mechanism_reached,
                                        du_factor=1.0, du_unit="m", fu_factor=1.0, fu_unit="N",
                                        solver="event_to_event", max_rotation=None,
                                        figsize=(14, 9)):
    """容量曲線(控制點位移 vs 底剪力)圖表頁, 降伏事件用紅點標出。"""
    u_disp = np.array(history_u) * du_factor
    F_disp = np.array(history_F) * fu_factor
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(u_disp, F_disp, color="#16a34a", linewidth=2)
    if event_log:
        ev_u = [ev["u"] * du_factor for ev in event_log]
        ev_F = [ev["F"] * fu_factor for ev in event_log]
        ax.scatter(ev_u, ev_F, color="#dc2626", zorder=5, s=30, label="Yield event")
        ax.legend(loc="lower right")
    ax.set_xlabel(f"Control point displacement ({du_unit})")
    ax.set_ylabel(f"Base shear ({fu_unit})")
    solver_label = "Newton-ish geometric-equilibrium-iteration solver" if solver == "converged" else "event-to-event solver (single linear solve per segment)"
    title = (f"Capacity Curve [{solver_label}] -- final: u={u_disp[-1]:.3f}{du_unit}, "
             f"F={F_disp[-1]:.3f}{fu_unit}, {len(event_log)} yield events")
    if mechanism_reached:
        reason = "geometric iteration failed to converge" if solver == "converged" else "mechanism reached (condensed stiffness near-singular)"
        title += f" [STOPPED EARLY: {reason}]"
    ax.set_title(title, fontsize=10, fontweight="bold")
    if max_rotation is not None and max_rotation > 0.0873:   # 0.0873 rad ≈ 5度
        import math
        ax.text(0.02, 0.98,
                f"⚠ max nodal rotation = {math.degrees(max_rotation):.1f}° -- exceeds small-angle\n"
                f"assumption range (~5°). Results beyond this point may not be physically valid\n"
                f"regardless of solver or convergence status.",
                transform=ax.transAxes, fontsize=9, color="#dc2626", va="top",
                bbox=dict(boxstyle="round", facecolor="#fef2f2", edgecolor="#dc2626"))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def build_pushover_hinge_diagram(f, hinge_summary, figsize=(11, 8.5)):
    """結構圖疊塑鉸圈圈, 跟網頁上「結構/荷載」分頁看到的畫法對應:
    藍色虛線空心圈=有設定塑鉸容量但這一步(這裡固定是最終狀態)還沒
    降伏, 紅色實心圈=已經降伏(半徑隨累積塑性轉角θp變大)。讓看報告的
    人不用自己對照桿件編號表格, 一眼就能看出「哪裡形成了塑鉸」。
    """
    fig, ax = plt.subplots(figsize=figsize)
    plot_structure(f, ax=ax, show_node_ids=True, show_member_ids=False, show_dimensions=False)

    for mid_str, hs in hinge_summary.items():
        mid = int(mid_str)
        if mid not in f.members:
            continue
        m = f.members[mid]
        ni, nj = f.nodes[m.node_i], f.nodes[m.node_j]
        dx, dy = nj.x - ni.x, nj.y - ni.y
        L = (dx ** 2 + dy ** 2) ** 0.5 or 1.0
        ux, uy = dx / L, dy / L
        offset = min(0.5, L * 0.18)
        ends = [(0, ni.x + ux * offset, ni.y + uy * offset, "i"),
                (1, nj.x - ux * offset, nj.y - uy * offset, "j")]
        for end_idx, ex, ey, label in ends:
            if hs["Mp"][end_idx] is None:
                continue
            yielded = hs["yielded"][end_idx]
            theta_p = hs["theta_p"][end_idx]
            r = min(0.08 + theta_p * 3.5, 0.30) if yielded else 0.08   # 資料座標(m), 跟畫布比例配合過
            if yielded:
                ax.add_patch(plt.Circle((ex, ey), r, facecolor="#dc2626", edgecolor="#dc2626",
                                         alpha=0.35, zorder=6))
            else:
                ax.add_patch(plt.Circle((ex, ey), r, facecolor="#2563eb", edgecolor="#2563eb",
                                         alpha=0.15, linestyle="--", linewidth=1.2, zorder=6))

    n_yielded = sum(1 for hs in hinge_summary.values() for e in (0, 1) if hs["yielded"][e])
    n_defined = sum(1 for hs in hinge_summary.values() for e in (0, 1) if hs["Mp"][e] is not None)
    ax.set_title(
        f"Final Plastic Hinge Locations (blue dashed = defined, not yet yielded; "
        f"red filled = yielded, {n_yielded}/{n_defined} hinges yielded)",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig


def build_pushover_hinge_pages(event_log, hinge_summary, du_factor=1.0, du_unit="m",
                                fu_factor=1.0, fu_unit="N", mu_factor=1.0, mu_unit="N·m",
                                figsize=(14, 9)):
    """降伏事件清單 + 最終塑鉸狀態兩張表格頁——讓看報告的人可以查到
    「哪根桿件哪一端、在容量曲線的哪個位置降伏」, 不用只能從圖上目測
    紅點對應到哪裡。"""
    event_rows = []
    for i, ev in enumerate(event_log):
        yielded_str = ", ".join(f"M{mid} end {'i' if e == 0 else 'j'}" for mid, e in ev["yielded"])
        event_rows.append([str(i + 1), _fmt(ev["u"] * du_factor, 4), _fmt(ev["F"] * fu_factor, 3), yielded_str])
    page1 = _table_page("Pushover Results (1/2): Yield Event Sequence", [
        (f"Yield events in order (displacement in {du_unit}, base shear in {fu_unit})",
         ["#", f"Displacement ({du_unit})", f"Base Shear ({fu_unit})", "Yielded member end(s)"], event_rows),
    ], figsize=figsize)

    hinge_rows = []
    for mid_str in sorted(hinge_summary.keys(), key=lambda x: int(x)):
        hs = hinge_summary[mid_str]
        for end_idx, end_label in [(0, "i"), (1, "j")]:
            if hs["Mp"][end_idx] is None:
                continue
            hinge_rows.append([
                mid_str, end_label, _fmt(hs["Mp"][end_idx] * mu_factor, 3),
                "yes" if hs["yielded"][end_idx] else "no",
                _fmt(hs["theta_p"][end_idx], 5),
                hs["performance_level"][end_idx] or "-",
            ])
    page2 = _table_page("Pushover Results (2/2): Final Hinge State", [
        (f"Final state of every defined plastic hinge (Mp in {mu_unit})",
         ["Member", "End", f"Mp ({mu_unit})", "Yielded", "Cumulative θp (rad)", "Performance Level"], hinge_rows),
    ], figsize=figsize)
    return [page1, page2]


def build_pushover_pdf_report(f, pushover_result, units=None) -> bytes:
    """Pushover結果的PDF報告: 完整輸入資料(含塑鉸容量) + 容量曲線圖 +
    降伏事件清單 + 最終塑鉸狀態表 + 最後一步的結構/變形/N/V/M總覽圖
    (直接重用plot_all(), 見_pushover_final_solve_result())——讓看報告
    的人不用另外拿到原始檔案就能重新檢驗整個側推過程跟結果。

    pushover_result: dict, 至少要有 history_u/history_F/event_log/
        hinge_summary/mechanism_reached/cum_forces/u_full_cum/cum_reaction
        這幾個key(webapi/main.py的_export_pushover_pdf()負責準備這個
        dict, 直接來自run_pushover()的回傳值)。

    已知限制(誠實記錄): 最後一步的N/V/M/變形圖用的是"最終那一步的內力
    狀態", 不是像容量曲線那樣涵蓋整個歷程——想看中間某一步的狀態,
    現在只有網頁上的「Pushover回放」看得到, PDF匯出目前只有最終狀態。
    """
    fu, ff, mu, mf = _force_moment_factors(units)
    du, df = _disp_factor(units)

    result = _pushover_final_solve_result(
        f, pushover_result["cum_forces"], pushover_result["u_full_cum"], pushover_result["cum_reaction"])

    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        fig = plot_all(f, result, figsize=(14, 9), force_factor=ff, force_unit=fu,
                       moment_factor=mf, moment_unit=mu, disp_factor=df, disp_unit=du)
        fig.suptitle("Pushover -- Final Step Overview (structure / deformed shape / N / V / M)",
                     fontsize=12, fontweight="bold", y=1.03)
        fig.subplots_adjust(top=0.90)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        curve_fig = build_pushover_capacity_curve_page(
            pushover_result["history_u"], pushover_result["history_F"], pushover_result["event_log"],
            pushover_result["mechanism_reached"], du_factor=df, du_unit=du, fu_factor=ff, fu_unit=fu,
            solver=pushover_result.get("solver", "event_to_event"),
            max_rotation=pushover_result.get("max_rotation"))
        pdf.savefig(curve_fig)
        plt.close(curve_fig)

        hinge_diagram_fig = build_pushover_hinge_diagram(f, pushover_result["hinge_summary"])
        pdf.savefig(hinge_diagram_fig)
        plt.close(hinge_diagram_fig)

        for page_fig in build_pushover_hinge_pages(
                pushover_result["event_log"], pushover_result["hinge_summary"],
                du_factor=df, du_unit=du, fu_factor=ff, fu_unit=fu, mu_factor=mf, mu_unit=mu):
            pdf.savefig(page_fig)
            plt.close(page_fig)

        for page_fig in build_input_data_pages(f, units):
            pdf.savefig(page_fig)
            plt.close(page_fig)

    return buf.getvalue()


def build_pdf_report(f, units=None, member_ids=None, include_member_diagrams=True) -> bytes:
    """f: 已經建好的 frame2d.Frame2D。units: 前端目前顯示單位設定,
    None的話用SI。member_ids: 如果有指定(選一個/多個/全選桿件),
    每根桿件會多附加自由體圖(含旁邊的結構縮圖, 驗證Fx/Fy/M平衡)頁面。
    include_member_diagrams: True(預設, 維持原本行為)時, 每根桿件
    的自由體圖後面還會再附一頁那根桿件自己的N/V/M/變形圖; 設成
    False可以只要自由體圖(含縮圖)、跳過N/V/M/變形這幾張, 讓PDF
    更精簡——例如想要「結構裡每一根桿件的自由體圖總覽」但不需要
    看每根的內力分布細節時用。回傳 PDF 的 bytes(六合一總覽圖 +
    完整輸入資料表格 + 求解結果表格 + 選擇性的桿件自由體圖頁),
    讓報告本身就能完全重現模型或手算/跨軟體驗證, 不用只能從圖上
    目測。"""
    fu, ff, mu, mf = _force_moment_factors(units)
    du, df = _disp_factor(units)
    result = solve(f)
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        fig = plot_all(f, result, figsize=(14, 9), force_factor=ff, force_unit=fu,
                       moment_factor=mf, moment_unit=mu, disp_factor=df, disp_unit=du)
        pdf.savefig(fig)
        plt.close(fig)

        for page_fig in build_input_data_pages(f, units):
            pdf.savefig(page_fig)
            plt.close(page_fig)

        result_fig = build_result_data_page(f, result, units)
        pdf.savefig(result_fig)
        plt.close(result_fig)

        if member_ids:
            for mid in member_ids:
                if mid not in f.members:
                    continue
                fbd_fig = _member_fbd_with_thumbnail_fig(f, result, mid, fu, ff, mu, mf)
                pdf.savefig(fbd_fig)
                plt.close(fbd_fig)

                if include_member_diagrams:
                    own_fig = plot_member_own_diagrams(f, result, mid, force_unit=fu, force_factor=ff,
                                                        moment_unit=mu, moment_factor=mf)
                    pdf.savefig(own_fig)
                    plt.close(own_fig)
    return buf.getvalue()


def _mode_shape_fig(f, modal_result, mode_indices, ncols=2):
    """畫一頁最多6個模態振型小圖(灰色虛線=未變形結構, 橘色實線=變形後形狀, 座標已經是
    modal_to_dict()算好的絕對座標, 這裡不重算)。mode_indices是這一頁要畫的模態編號
    (modal_result['modes']裡的index, 1起算)清單。"""
    n = len(mode_indices)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4.2 * nrows), squeeze=False)
    modes_by_index = {m["index"]: m for m in modal_result["modes"]}
    for k, idx in enumerate(mode_indices):
        ax = axes[k // ncols][k % ncols]
        m = modes_by_index[idx]
        for member in f.members.values():
            ni, nj = f.nodes[member.node_i], f.nodes[member.node_j]
            ax.plot([ni.x, nj.x], [ni.y, nj.y], color="#bbbbbb", lw=1.2, ls="--", zorder=1)
        for mid_str, curve in m["curves"].items():
            ax.plot(curve["X"], curve["Y"], color="#ea580c", lw=1.8, zorder=2)
        ax.set_title(f"Mode {idx}: T={m['period']:.4f}s, f={m['frequency']:.4f}Hz "
                     f"(Γx={m['gamma_x']:.3f}, Γy={m['gamma_y']:.3f})", fontsize=10)
        ax.set_aspect("equal", adjustable="datalim")
        ax.tick_params(labelsize=8)
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    fig.suptitle("Mode Shapes (scaled for visibility, not true amplitude)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def build_modal_pdf_report(f, modal_result, units=None) -> bytes:
    """模態分析的PDF報告: 週期/頻率/參與係數/有效質量比表 + 振型圖(每頁最多6個, 直接用
    modal_to_dict()已經算好的變形曲線座標, 不重算) + 完整輸入資料(含質量設定)。

    modal_result: dict, `frame2d.modal.modal_to_dict()`的回傳值(webapi/main.py的
    /export/pdf(modal)路徑負責準備, 重新跑一次eigen()拿到)。
    """
    mu_ = _unit_label(units, "mass", "kg")

    mode_rows = [[
        str(m["index"]), f"{m['period']:.5f}", f"{m['frequency']:.4f}", f"{m['omega']:.4f}",
        f"{m['gamma_x']:.4f}", f"{(m['ratio_x'] or 0) * 100:.1f}%" if m["ratio_x"] is not None else "-",
        f"{(m['cum_total_x'] or 0) * 100:.1f}%" if m["cum_total_x"] is not None else "-",
        f"{m['gamma_y']:.4f}", f"{(m['ratio_y'] or 0) * 100:.1f}%" if m["ratio_y"] is not None else "-",
        f"{(m['cum_total_y'] or 0) * 100:.1f}%" if m["cum_total_y"] is not None else "-",
    ] for m in modal_result["modes"]]
    table_fig = _table_page(
        f"Modal Analysis Results ({modal_result['mass_kind']} mass, {modal_result['n_modes']} modes; "
        f"free mass x={_from_si(units,'mass',modal_result['mass_free']['x'],'kg'):.4g}{mu_}, "
        f"y={_from_si(units,'mass',modal_result['mass_free']['y'],'kg'):.4g}{mu_})",
        [("", ["Mode", "T (s)", "f (Hz)", "ω (rad/s)", "Γx", "Ratio x", "Cum(total) x",
              "Γy", "Ratio y", "Cum(total) y"], mode_rows)], figsize=(14, 3 + 0.35 * len(mode_rows)))

    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        pdf.savefig(table_fig, bbox_inches="tight")
        plt.close(table_fig)

        all_idx = [m["index"] for m in modal_result["modes"]]
        for i in range(0, len(all_idx), 6):
            fig = _mode_shape_fig(f, modal_result, all_idx[i:i + 6])
            pdf.savefig(fig)
            plt.close(fig)

        mass_fig = build_mass_data_page(f, units)
        pdf.savefig(mass_fig, bbox_inches="tight")
        plt.close(mass_fig)

        for page_fig in build_input_data_pages(f, units):
            pdf.savefig(page_fig)
            plt.close(page_fig)

    return buf.getvalue()


def _rsa_spectrum_curve_fig(rsa_result, units=None):
    """反應譜曲線圖(週期-加速度), 疊上每個模態實際落在譜上的位置(紅點+模態編號標籤)。
    T軸固定用秒(跟前端一致, 沒有另外的週期單位選項), Sa軸依units換算。"""
    accu = _unit_label(units, "accel", "m/s²")
    T = np.array(rsa_result["curve"]["T"])
    Sa = np.array([_from_si(units, "accel", v, "m/s²") for v in rsa_result["curve"]["Sa"]])
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(T, Sa, color="#0891b2", lw=1.6)
    for m in rsa_result["modes"]:
        sa_disp = _from_si(units, "accel", m["sa"], "m/s²")
        ax.plot(m["period"], sa_disp, "o", color="#dc2626", ms=6, zorder=3)
        ax.annotate(f"Mode {m['index']}", (m["period"], sa_disp), textcoords="offset points",
                    xytext=(6, 6), fontsize=8, color="#dc2626")
    ax.set_xlabel("Period T (s)")
    ax.set_ylabel(f"Sa ({accu})")
    ax.set_title(f"Response Spectrum ({rsa_result['direction']}-direction, {rsa_result['combine']}, "
                 f"damping={rsa_result['damping']})")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def build_rsa_result_data_page(rsa_result, units=None):
    """反應譜結果表格頁: 各模態的組合貢獻表 + 節點位移表 + 桿件端點內力表。跟
    build_result_data_page() 不一樣, 這裡的資料來自 `rsa_to_dict()` 的 dict, 不是
    frame2d.postprocess 的 solve result 物件——反應譜結果本身就是組合過的(SRSS/CQC),
    沒有對應的單一SolveResult可以重用。"""
    du = _unit_label(units, "disp", "m")
    fu = _unit_label(units, "force", "N")
    mu = _unit_label(units, "moment", "N·m")
    accu = _unit_label(units, "accel", "m/s²")

    mode_rows = [[
        str(m["index"]), f"{m['period']:.5f}", f"{m['frequency']:.4f}",
        _fmt(_from_si(units, "accel", m["sa"], "m/s²")), f"{m['gamma']:.4f}",
        f"{m['eff_ratio'] * 100:.1f}%" if m["eff_ratio"] is not None else "-",
        f"{m['cum_ratio_total'] * 100:.1f}%" if m["cum_ratio_total"] is not None else "-",
        _fmt(_from_si(units, "force", m["modal_base_shear"], "N")),
    ] for m in rsa_result["modes"]]
    page1 = _table_page(
        f"Response Spectrum Results ({rsa_result['direction']}-dir, {rsa_result['combine']}, "
        f"damping={rsa_result['damping']}, {rsa_result['mass_kind']} mass, {rsa_result['n_modes']} modes): "
        f"Base Shear = {_fmt(_from_si(units,'force',rsa_result['base_shear'],'N'))} {fu}, "
        f"Cumulative Mass Ratio = {rsa_result['cum_ratio_total'] * 100:.1f}%",
        [("", ["Mode", "T (s)", "f (Hz)", f"Sa ({accu})", "Γ", "Eff. Ratio", "Cum(total)",
              f"Modal Base Shear ({fu})"], mode_rows)], figsize=(14, 3 + 0.35 * len(mode_rows)))

    node_rows = [[
        nid, _fmt(_from_si(units, "disp", v["ux"], "m")), _fmt(_from_si(units, "disp", v["uy"], "m")),
    ] for nid, v in sorted(rsa_result["nodes"].items(), key=lambda kv: int(kv[0]))]
    member_rows = [[
        mid, _fmt(_from_si(units, "force", v["Fx_i"], "N")), _fmt(_from_si(units, "force", v["Fy_i"], "N")),
        _fmt(_from_si(units, "moment", v["M_i"], "N·m")), _fmt(_from_si(units, "force", v["Fx_j"], "N")),
        _fmt(_from_si(units, "force", v["Fy_j"], "N")), _fmt(_from_si(units, "moment", v["M_j"], "N·m")),
    ] for mid, v in sorted(rsa_result["members"].items(), key=lambda kv: int(kv[0]))]
    page2 = _table_page("Response Spectrum Results: Node Displacements / Member End Forces", [
        (f"Node Displacements ({du}, SRSS/CQC combined -- always non-negative-ish magnitude, not a signed value)",
         ["Node ID", f"ux ({du})", f"uy ({du})"], node_rows),
        (f"Member End Forces ({fu}, {mu}, combined)",
         ["Member ID", f"Nᵢ ({fu})", f"Vᵢ ({fu})", f"Mᵢ ({mu})", f"Nⱼ ({fu})", f"Vⱼ ({fu})", f"Mⱼ ({mu})"], member_rows),
    ])
    return [page1, page2]


def build_rsa_pdf_report(f, rsa_result, units=None) -> bytes:
    """反應譜分析的PDF報告: 反應譜曲線圖(疊模態點) + 模態組合表(含基底剪力) + 節點位移/
    桿件端點內力表 + 完整輸入資料(含質量設定)。

    rsa_result: dict, `frame2d.spectrum.rsa_to_dict()`的回傳值(webapi/main.py的
    /export/pdf(rsa)路徑負責準備, 重新跑一次spectrum_analysis()拿到)。
    """
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        curve_fig = _rsa_spectrum_curve_fig(rsa_result, units)
        pdf.savefig(curve_fig, bbox_inches="tight")
        plt.close(curve_fig)

        for page_fig in build_rsa_result_data_page(rsa_result, units):
            pdf.savefig(page_fig, bbox_inches="tight")
            plt.close(page_fig)

        mass_fig = build_mass_data_page(f, units)
        pdf.savefig(mass_fig, bbox_inches="tight")
        plt.close(mass_fig)

        for page_fig in build_input_data_pages(f, units):
            pdf.savefig(page_fig)
            plt.close(page_fig)

    return buf.getvalue()


def _hinge_label_en(h):
    """`cyclic_to_dict()`/`seismic_to_dict()`的塑鉸label是中文(例如"M0 i端"), matplotlib預設
    字型(DejaVu Sans)不支援中文字元, 直接塞進圖表會變成缺字框——跟 build_input_data_pages()
    刻意避開中文的慣例一致, 這裡改用「Member 0, end i」這種英文格式重新組字串, 不能直接沿用
    h['label']。"""
    return f"Member {h['member']}, end {'i' if h['end'] == 0 else 'j'}"


def _cyclic_loop_fig(cyclic_result, units=None):
    """整體遲滯迴圈圖(控制節點位移 vs 力), 黃色三角形標記每個降伏事件的位置(跟網頁
    drawCyclicLoop()的視覺慣例一致: 降伏事件是路徑上的轉折點)。"""
    du = _unit_label(units, "disp", "m")
    fu = _unit_label(units, "force", "N")
    u = np.array([_from_si(units, "disp", v, "m") for v in cyclic_result["u"]])
    F = np.array([_from_si(units, "force", v, "N") for v in cyclic_result["F"]])
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.plot(u, F, color="#7c3aed", lw=1.4, zorder=1)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.axvline(0, color="#999999", lw=0.6)
    yield_events = [e for e in cyclic_result["events"] if e["kind"] == "yield"]
    if yield_events:
        yu = [_from_si(units, "disp", e["u"], "m") for e in yield_events]
        yf = [_from_si(units, "force", e["F"], "N") for e in yield_events]
        ax.scatter(yu, yf, marker="^", color="#f59e0b", s=45, zorder=3, label="Yield event")
        ax.legend(loc="best", fontsize=9)
    ax.set_xlabel(f"Displacement ({du})")
    ax.set_ylabel(f"Force ({fu})")
    ax.set_title(f"Cyclic Hysteresis Loop (control node(s) {cyclic_result['control_nodes']}, "
                 f"{cyclic_result['direction']}-direction)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _cyclic_hinge_fig(cyclic_result, hinges_subset, units=None, ncols=2):
    """一頁最多6個塑鉸的 M-θp 小圖(跟 _mode_shape_fig() 同一種每頁最多6個的分頁慣例)。"""
    mu = _unit_label(units, "moment", "N·m")
    n = len(hinges_subset)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4.0 * nrows), squeeze=False)
    for k, h in enumerate(hinges_subset):
        ax = axes[k // ncols][k % ncols]
        theta_p = np.array(h["theta_p"]) * 1000
        M = np.array([_from_si(units, "moment", v, "N·m") for v in h["M"]])
        ax.plot(theta_p, M, color="#ea580c", lw=1.4)
        ax.axhline(0, color="#999999", lw=0.6)
        ax.axvline(0, color="#999999", lw=0.6)
        ax.set_title(f"{_hinge_label_en(h)} (yielded {h['n_yield']}x)", fontsize=10)
        ax.set_xlabel("θp (mrad)", fontsize=9)
        ax.set_ylabel(f"M ({mu})", fontsize=9)
        ax.tick_params(labelsize=8)
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    fig.suptitle("Hinge Moment - Plastic Rotation (M-θp)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def build_cyclic_result_data_page(cyclic_result, units=None):
    """循環分析的每一級幅值耗能表 + 塑鉸降伏順序表。"""
    du = _unit_label(units, "disp", "m")
    fu = _unit_label(units, "force", "N")
    mu = _unit_label(units, "moment", "N·m")

    loop_rows = [[
        _fmt(_from_si(units, "disp", loop["amplitude"], "m")),
        _fmt(_from_si(units, "force", loop["F_max"], "N")) if loop["F_max"] is not None else "-",
        _fmt(_from_si(units, "moment", loop["energy"], "N·m")) if loop["energy"] is not None else "-",
        f"{loop['xi_eq']:.4f}" if loop["xi_eq"] is not None else "-",
    ] for loop in cyclic_result["loops"]]
    hinge_rows = [[
        str(i + 1), _hinge_label_en(h), _fmt(_from_si(units, "disp", h["first_yield"]["u"], "m")),
        _fmt(_from_si(units, "force", h["first_yield"]["F"], "N")), str(h["n_yield"]),
        _fmt(_from_si(units, "moment", h["work"], "N·m")),
    ] for i, h in enumerate(cyclic_result["hinges"])]
    return _table_page(
        f"Cyclic Analysis Results ({cyclic_result['n_steps']} steps, "
        f"total plastic work = {_fmt(_from_si(units,'moment',cyclic_result['total_plastic_work'],'N·m'))} {mu})",
        [
            (f"Energy per Complete Loop at Each Amplitude ({du}, {fu}, {mu})",
             ["Amplitude", "F_max", "Energy", "ξ_eq"], loop_rows),
            (f"Hinge Yield Sequence ({du}, {fu}, {mu})",
             ["#", "Hinge", "Yield u", "Yield F", "N. Yields", "Plastic Work"], hinge_rows),
        ])


def build_cyclic_pdf_report(f, cyclic_result, units=None) -> bytes:
    """循環(遲滯)分析的PDF報告: 整體遲滯迴圈圖 + 每級幅值耗能表/塑鉸降伏順序表 + 塑鉸M-θp
    小圖(每頁最多6個)+ 完整輸入資料(含質量設定)。

    cyclic_result: dict, `frame2d.cyclic.cyclic_to_dict()`的回傳值(webapi/main.py的
    /export/pdf(cyclic)路徑負責準備, 重新跑一次cyclic_analysis()拿到)。
    """
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        loop_fig = _cyclic_loop_fig(cyclic_result, units)
        pdf.savefig(loop_fig, bbox_inches="tight")
        plt.close(loop_fig)

        result_fig = build_cyclic_result_data_page(cyclic_result, units)
        pdf.savefig(result_fig, bbox_inches="tight")
        plt.close(result_fig)

        hinges = cyclic_result["hinges"]
        for i in range(0, len(hinges), 6):
            fig = _cyclic_hinge_fig(cyclic_result, hinges[i:i + 6], units)
            pdf.savefig(fig)
            plt.close(fig)

        mass_fig = build_mass_data_page(f, units)
        pdf.savefig(mass_fig, bbox_inches="tight")
        plt.close(mass_fig)

        for page_fig in build_input_data_pages(f, units):
            pdf.savefig(page_fig)
            plt.close(page_fig)

    return buf.getvalue()


def _seismic_history_fig(seismic_result, units=None):
    """地震時程圖: 地面加速度(g)+ 控制節點相對位移雙圖(上下疊, 跟網頁drawSeismicHistory()
    同一種版面)。"""
    du = _unit_label(units, "disp", "m")
    t = np.array(seismic_result["t"])
    ag_g = np.array(seismic_result["ground_motion"]) / 9.80665
    disp = np.array([_from_si(units, "disp", v, "m") for v in seismic_result["control_disp"]])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8))
    ax1.plot(t, ag_g, color="#64748b", lw=1.0)
    ax1.axhline(0, color="#999999", lw=0.5)
    ax1.set_ylabel("Ground accel. (g)")
    ax1.set_title("Ground Acceleration ag(t)")
    ax1.grid(alpha=0.3)
    ax2.plot(t, disp, color="#0891b2", lw=1.1)
    ax2.axhline(0, color="#999999", lw=0.5)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel(f"Rel. disp ({du})")
    ax2.set_title(f"Control Node {seismic_result['control_node']} Relative Displacement")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _seismic_hysteresis_fig(seismic_result, units=None):
    """全域遲滯圖: 控制節點位移 vs 有效慣性力合力(不是嚴格的支承反力, 見
    `seismic.seismic_to_dict()` 的說明)。"""
    du = _unit_label(units, "disp", "m")
    fu = _unit_label(units, "force", "N")
    disp = np.array([_from_si(units, "disp", v, "m") for v in seismic_result["control_disp"]])
    force = np.array([_from_si(units, "force", v, "N") for v in seismic_result["equivalent_force"]])
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.plot(disp, force, color="#7c3aed", lw=1.0)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.set_xlabel(f"Displacement ({du})")
    ax.set_ylabel(f"Equivalent inertial force ({fu})")
    ax.set_title(f"Displacement vs Equivalent Inertial Force (control node {seismic_result['control_node']}, "
                 "NOT a strict support reaction -- see report notes)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _seismic_energy_fig(seismic_result, units=None):
    """能量平衡圖: 動能+阻尼耗能+桿件內力作功疊加(stackplot), 對照外力作功(虛線)——虛線應該
    精確蓋在疊加區域的頂端, 這張圖本身就是驗證這組結果最直接的方式。"""
    mu = _unit_label(units, "moment", "N·m")
    t = np.array(seismic_result["energy"]["t"])
    KE = np.array([_from_si(units, "moment", v, "N·m") for v in seismic_result["energy"]["KE"]])
    Wdamp = np.array([_from_si(units, "moment", v, "N·m") for v in seismic_result["energy"]["Wdamp"]])
    Wint = np.array([_from_si(units, "moment", v, "N·m") for v in seismic_result["energy"]["Wint"]])
    Wext = np.array([_from_si(units, "moment", v, "N·m") for v in seismic_result["energy"]["Wext"]])
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.stackplot(t, KE, Wdamp, Wint, labels=["Kinetic", "Damping dissipated", "Member internal work"],
                colors=["#0891b2", "#64748b", "#ea580c"], alpha=0.85)
    ax.plot(t, Wext, color="#111827", lw=1.4, ls="--", label="External (effective EQ force) work")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(f"Cumulative energy ({mu})")
    ax.set_title("Energy Balance: KE + Damping + Member Internal Work vs External Work")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def build_seismic_result_data_page(seismic_result, units=None):
    """非線性地震分析的摘要頁: 分析設定/結果摘要表(含能量平衡自我檢核)+ 降伏塑鉸表。"""
    du = _unit_label(units, "disp", "m")
    fu = _unit_label(units, "force", "N")
    mu = _unit_label(units, "moment", "N·m")
    r = seismic_result
    n_yield = sum(1 for e in r["events"] if e["kind"] == "yield")
    n_unload = sum(1 for e in r["events"] if e["kind"] == "unload")
    e_final = r["energy"]["Wext"][-1]
    e_check = r["energy"]["KE"][-1] + r["energy"]["Wdamp"][-1] + r["energy"]["Wint"][-1]

    summary_rows = [
        ["Control node / direction", f"{r['control_node']} / {r['direction']}"],
        ["Mass kind / gravity preload", f"{r['mass_kind']} / {'yes' if r['apply_gravity_loads'] else 'no'}"],
        ["T1 / ζ_target", f"{r['period1']:.4f} s / {r['zeta_target']}"],
        ["Rayleigh α / β", f"{r['alpha']:.5f} / {r['beta']:.3e}"],
        ["Δt / steps", f"{r['dt']:.5g} s / {r['n_steps']}"],
        [f"Peak displacement ({du})", _fmt(_from_si(units, "disp", r["peak_displacement"], "m"))],
        ["Yield / unload events", f"{n_yield} / {n_unload}"],
        [f"Energy balance: external work ({mu})", _fmt(_from_si(units, "moment", e_final, "N·m"))],
        [f"Energy balance: KE+damping+internal ({mu}, should match above exactly)",
         _fmt(_from_si(units, "moment", e_check, "N·m"))],
    ]
    hinge_rows = [[_hinge_label_en(h), str(h["n_yield"]),
                  _fmt(_from_si(units, "moment", h["work_final"], "N·m"))] for h in r["hinges"]]
    return _table_page(
        "Nonlinear Seismic Analysis Results",
        [("Summary", ["Item", "Value"], summary_rows),
         (f"Yielded Hinges (plastic work in {mu})", ["Hinge", "N. Yields", "Plastic Work"], hinge_rows)],
        figsize=(14, 6 + 0.3 * max(1, len(hinge_rows))))


def build_seismic_pdf_report(f, seismic_result, units=None) -> bytes:
    """非線性地震反應分析的PDF報告: 地震時程圖(地面加速度+位移)+ 全域遲滯圖 + 能量平衡圖 +
    摘要表(含能量平衡自我檢核)+ 塑鉸M-θp小圖(每頁最多6個, 重用 `_cyclic_hinge_fig()`——
    seismic_to_dict()跟cyclic_to_dict()的塑鉸dict形狀相容, 都有member/end/M/theta_p/
    n_yield)+ 質量設定頁 + 完整輸入資料頁。

    seismic_result: dict, `frame2d.seismic.seismic_to_dict()`的回傳值(webapi/main.py的
    /export/pdf(seismic)路徑負責準備, 重新跑一次nonlinear_seismic_web_analysis()拿到)。
    """
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        for fig in (_seismic_history_fig(seismic_result, units),
                   _seismic_hysteresis_fig(seismic_result, units),
                   _seismic_energy_fig(seismic_result, units)):
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        result_fig = build_seismic_result_data_page(seismic_result, units)
        pdf.savefig(result_fig, bbox_inches="tight")
        plt.close(result_fig)

        hinges = seismic_result["hinges"]
        for i in range(0, len(hinges), 6):
            hinge_fig = _cyclic_hinge_fig(seismic_result, hinges[i:i + 6], units)
            pdf.savefig(hinge_fig)
            plt.close(hinge_fig)

        mass_fig = build_mass_data_page(f, units)
        pdf.savefig(mass_fig, bbox_inches="tight")
        plt.close(mass_fig)

        for page_fig in build_input_data_pages(f, units):
            pdf.savefig(page_fig)
            plt.close(page_fig)

    return buf.getvalue()
