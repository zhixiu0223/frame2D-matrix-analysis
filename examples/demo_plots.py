"""
Demo: 畫出六合一視覺化 (結構圖/受力圖/軸力圖/剪力圖/彎矩圖/變形圖)

在 Pydroid3 上直接執行, plt.show() 會跳出內建繪圖檢視畫面。
在 Termux / 無顯示器環境, 改用 plt.savefig(...) (下面已經示範)。
"""
import matplotlib
import matplotlib.pyplot as plt
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all

# ---- Case-04: 單跨門型鋼架, 水平力+梁均佈載重 ----
# 斷面用 sd_framework 自己的 EI_numeric=15000.0 慣例(見tests/test_case04_vs_sdframework.py),
# 已經跟 sd_framework/anastruct 驗證過吻合。E=1.0,I=EI直接當EI用,A給近似剛體的大數字。
E, I, A = 1.0, 15000.0, 1e8
H, L, P = 4.0, 6.0, 12.0

f = Frame2D()
f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, L, H).add_node(3, L, 0)
f.add_section('sec', E=E, I=I, A=A)
f.add_member(0, node_i=0, node_j=1, section='sec')
f.add_member(1, node_i=1, node_j=2, section='sec')
f.add_member(2, node_i=3, node_j=2, section='sec')
f.fix(0)
f.fix(3)
f.point_load(1, fx=P)
f.distributed_load(1, w=-24.0)   # Case-04.5的梁均佈載重, 已驗證過的數值

result = solve(f)
fig = plot_all(f, result)
fig.savefig('case04_six_views.png', dpi=120)
print("已存檔: case04_six_views.png")

# 若在有GUI的環境(Pydroid3等), 取消下一行的註解可直接跳出檢視畫面:
# plt.show()
