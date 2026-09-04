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

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from frame2d import solve
from frame2d.plotting import plot_all, plot_member_fbd, plot_member_own_diagrams


def _fig_to_png_base64(fig, dpi=110):
    """把matplotlib figure轉成base64編碼的PNG字串, 用在「匯出前先在
    瀏覽器預覽」這個功能——不用另外存檔案再回傳URL, 直接把圖片編碼
    塞進JSON回應裡, 前端<img src="data:image/png;base64,...">就能
    直接顯示, 不用額外的靜態檔案伺服流程。"""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def build_fbd_previews(f, member_ids):
    """回傳每根指定桿件的自由體圖預覽(base64 PNG), 給前端在真正
    匯出PDF之前先秀出來讓使用者確認、可以個別排除不想要的桿件。
    回傳list of {"member_id":int, "image_base64":str}, 找不到的
    member_id直接跳過(不報錯, 讓前端自己比對哪些真的有回傳)。"""
    result = solve(f)
    previews = []
    for mid in member_ids:
        if mid not in f.members:
            continue
        ax = plot_member_fbd(f, result, mid)
        previews.append({"member_id": mid, "image_base64": _fig_to_png_base64(ax.figure)})
        plt.close(ax.figure)
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
    ] for m in members]
    support_rows = [[str(s.node), _fmt(s.ux), _fmt(s.uy), _fmt(s.rot)] for s in supports]
    page2 = _table_page("Input Data (2/3): Member Connectivity / Support Conditions", [
        ("Member Connectivity", ["Member ID", "Node i", "Node j", "Section", "Type", "Release i", "Release j"], member_rows),
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


def build_pdf_report(f, units=None, member_ids=None) -> bytes:
    """f: 已經建好的 frame2d.Frame2D。units: 前端目前顯示單位設定,
    None的話用SI。member_ids: 如果有指定(選一個/多個/全選桿件),
    每根桿件會多兩頁: 自由體圖(驗證Fx/Fy/M平衡)+ 這根桿件自己的
    N/V/M/變形圖, 方便把特定桿件抓出來單獨檢查或設計。回傳 PDF 的
    bytes(六合一總覽圖 + 完整輸入資料表格 + 求解結果表格 + 選擇性
    的桿件自由體圖頁), 讓報告本身就能完全重現模型或手算/跨軟體
    驗證, 不用只能從圖上目測。"""
    result = solve(f)
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        fig = plot_all(f, result, figsize=(14, 9))
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
                fbd_ax = plot_member_fbd(f, result, mid)
                pdf.savefig(fbd_ax.figure)
                plt.close(fbd_ax.figure)

                own_fig = plot_member_own_diagrams(f, result, mid)
                pdf.savefig(own_fig)
                plt.close(own_fig)
    return buf.getvalue()
