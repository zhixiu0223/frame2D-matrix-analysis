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
import io

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from frame2d import solve
from frame2d.plotting import plot_all


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


def build_input_data_pages(f):
    """把Frame2D f 的所有輸入資料(節點/斷面/桿件/支承/三種載重)排成
    表格頁, 回傳figure的list——目的是讓報告本身就能完全重現模型或
    手算驗證, 不用只能從六合一圖上目測。表格標籤全部用英文, 跟
    plot_all()既有的慣例一致——這不是隨便選的, 是刻意避開matplotlib
    預設字型(DejaVu Sans)不支援中文字元的問題: 這份報告要能在Termux/
    Colab/Cloud Run等沒有裝中文字型的環境正常產生, 不能假設部署環境
    一定有CJK字型可用。"""
    nodes = sorted(f.nodes.values(), key=lambda n: n.id)
    sections = sorted(f.sections.values(), key=lambda s: s.name)
    members = sorted(f.members.values(), key=lambda m: m.id)
    supports = sorted(f.supports, key=lambda s: s.node)
    point_loads = sorted(f.point_loads, key=lambda p: p.node)
    dist_loads = sorted(f.distributed_loads, key=lambda d: d.member)
    member_point_loads = sorted(f.member_point_loads, key=lambda p: p.member)

    node_rows = [[str(n.id), _fmt(n.x), _fmt(n.y)] for n in nodes]
    section_rows = [[s.name, _fmt(s.E), _fmt(s.I), _fmt(s.A)] for s in sections]
    page1 = _table_page("Input Data (1/3): Node Coordinates / Section Properties (SI units: m, Pa, m^4, m^2)", [
        ("Node Coordinates", ["Node ID", "x (m)", "y (m)"], node_rows),
        ("Section Properties", ["Section Name", "E (Pa)", "I (m^4)", "A (m^2)"], section_rows),
    ])

    member_rows = [[
        str(m.id), str(m.node_i), str(m.node_j), m.section, m.member_type,
        "yes" if m.release_i else "no", "yes" if m.release_j else "no",
    ] for m in members]
    support_rows = [[str(s.node), _fmt(s.ux), _fmt(s.uy), _fmt(s.rot)] for s in supports]
    page2 = _table_page("Input Data (2/3): Member Connectivity / Support Conditions", [
        ("Member Connectivity", ["Member ID", "Node i", "Node j", "Section", "Type", "Release i", "Release j"], member_rows),
        ("Support Conditions (blank = free, 0 = fixed, nonzero = prescribed displacement)",
         ["Node ID", "ux (m)", "uy (m)", "rot (rad)"], support_rows),
    ])

    pl_rows = [[str(p.node), _fmt(p.fx), _fmt(p.fy), _fmt(p.m)] for p in point_loads]
    dl_rows = [[
        str(d.member), _fmt(d.w_start), _fmt(d.w_end) if d.w_end is not None else "= w1",
        _fmt(d.x_start) if d.x_start is not None else "0",
        _fmt(d.x_end) if d.x_end is not None else "full length",
        d.direction, _fmt(d.angle_deg) if d.angle_deg is not None else "-",
    ] for d in dist_loads]
    mpl_rows = [[
        str(p.member), _fmt(p.a), p.direction,
        _fmt(p.fx) if p.direction == "local" else "-",
        _fmt(p.fy) if p.direction == "local" else "-",
        _fmt(p.F) if p.direction == "global" else "-",
        _fmt(p.angle_deg) if p.direction == "global" else "-",
        _fmt(p.m),
    ] for p in member_point_loads]
    page3 = _table_page("Input Data (3/3): Nodal Point Loads / Distributed Loads / Member Point Loads", [
        ("Nodal Point Loads", ["Node ID", "fx (N)", "fy (N)", "m (N*m)"], pl_rows),
        ("Distributed Loads (w measured along member length)",
         ["Member ID", "w1 (start)", "w2 (end)", "x start (m)", "x end (m)", "direction", "angle (deg)"], dl_rows),
        ("Member Internal Point Loads (a = distance from node i along member axis)",
         ["Member ID", "a (m)", "direction", "fx (N)", "fy (N)", "F (N)", "angle (deg)", "m (N*m)"], mpl_rows),
    ])
    return [page1, page2, page3]


def build_result_data_page(f, result):
    """反力/位移結果表格頁, 讓報告本身就能對照手算或別的軟體算出的
    反力/位移數值, 不用只能從圖上目測峰值。"""
    node_rows = []
    for n in sorted(f.nodes.values(), key=lambda n: n.id):
        ux, uy, rot = f.dofs_of(n.id)
        node_rows.append([
            str(n.id), _fmt(result.displacements[ux]), _fmt(result.displacements[uy]),
            _fmt(result.displacements[rot]), _fmt(result.reactions[ux]),
            _fmt(result.reactions[uy]), _fmt(result.reactions[rot]),
        ])
    return _table_page("Solve Results: Node Displacements / Reactions (SI units: m, rad, N, N*m)", [
        ("Node Displacements & Reactions (reaction is nonzero only where that DOF is restrained)",
         ["Node ID", "ux (m)", "uy (m)", "rot (rad)", "Rx (N)", "Ry (N)", "M (N*m)"], node_rows),
    ])


def build_pdf_report(f) -> bytes:
    """f: 已經建好的 frame2d.Frame2D。回傳 PDF 的 bytes(多頁: 六合一
    總覽圖 + 完整輸入資料表格 + 求解結果表格), 讓報告本身就能完全
    重現模型或手算/跨軟體驗證, 不用只能從圖上目測。"""
    result = solve(f)
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        fig = plot_all(f, result, figsize=(14, 9))
        pdf.savefig(fig)
        plt.close(fig)

        for page_fig in build_input_data_pages(f):
            pdf.savefig(page_fig)
            plt.close(page_fig)

        result_fig = build_result_data_page(f, result)
        pdf.savefig(result_fig)
        plt.close(result_fig)
    return buf.getvalue()
