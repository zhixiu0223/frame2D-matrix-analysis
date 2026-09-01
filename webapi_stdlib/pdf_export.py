"""
PDF 匯出: 直接重用 frame2d/plotting.py 已經驗證過的 plot_all()
(六合一總覽圖: 結構/受力/變形/N圖/V圖/M圖),不在這層重新畫任何圖,
避免跟 plotting.py 的正負號慣例(拉力側等)兩邊分別維護、互相漂移。

matplotlib 用 Agg(非互動式)backend,在沒有顯示器的伺服器環境
(Termux/Colab/Cloud Run)都能正常產生圖檔。
"""
import io

import matplotlib
matplotlib.use("Agg")

from frame2d import solve
from frame2d.plotting import plot_all


def build_pdf_report(f) -> bytes:
    """f: 已經建好的 frame2d.Frame2D。回傳 PDF 的 bytes。"""
    result = solve(f)
    fig = plot_all(f, result, figsize=(14, 9))
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return buf.getvalue()
