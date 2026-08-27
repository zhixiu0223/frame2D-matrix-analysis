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
E, I, A = 200e6, 8e-5, 1e-2
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
f.distributed_load(1, w=-3.0)

result = solve(f)
fig = plot_all(f, result)
fig.savefig('case04_six_views.png', dpi=120)
print("已存檔: case04_six_views.png")

# 若在有GUI的環境(Pydroid3等), 取消下一行的註解可直接跳出檢視畫面:
# plt.show()
