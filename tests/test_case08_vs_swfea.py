"""
驗證案例3: 兩層兩跨側移鋼架 (九節點/十桿件), 對照第三方工具 SW FEA (swframe2d,
Android app) 的 PDF 分析報告, 純點載重版本(Report.pdf) + 含均佈載重版本(Report-udl.pdf)。

幾何/斷面/載重數值直接取自使用者上傳的 SW FEA .frame (SQLite) 檔案:
  E=150 GPa, I=1e-6 m^4, A=0.05 m^2 (全桿件同斷面)
  水平點載重: node3(1F, x=0) Fx=+15kN(向右), node6(roof, x=0) Fx=+10kN(向右)
  均佈載重(僅UDL版本): 1F梁(F3,F4) 10 kN/m 向下, 屋頂梁(F8,F9) 12 kN/m 向下

**方向慣例: 全部直接吻合, 不需要任何翻轉或取絕對值。**
開發過程中一度誤判方向(見git歷史/對話記錄), 後來根據SW FEA app畫面上的箭頭
(不是對話框裡轉盤圖示的猜測)重新確認: 點載重180°=沿全域+x方向(向右), 均佈
載重方向就是直覺的「向下」——這跟frame2d、跟sd_framework/anastruct都是同一套
標準工程慣例(Y向上為正, 右手座標系), 完全不需要任何特例或翻轉。之前記錄的
"SW FEA內部方向不一致"是誤判, 一部分原因是最初PDF反力數字的Ry/M正負號被
抄錯, 另一部分是把點載重方向看反了; 這裡已用位移(dX,dY, 全部精確吻合到
小數點後3位, 比反力更不容易受慣例混淆)重新交叉確認過, 兩個版本(純點載重/
含UDL)的Rx,Ry,M三個反力分量也都直接吻合, 不需要任何取負號的特例處理。
"""
from frame2d import Frame2D, solve

E, I, A = 150e9, 1e-6, 0.05

nodes = {0: (0, 0), 1: (5, 0), 2: (12, 0),
         3: (0, 4), 4: (5, 4), 5: (12, 4),
         6: (0, 7.5), 7: (5, 7.5), 8: (12, 7.5)}
members = [(0, 0, 3), (1, 1, 4), (2, 5, 2), (3, 3, 4), (4, 4, 5),
           (5, 6, 3), (6, 7, 4), (7, 8, 5), (8, 6, 7), (9, 7, 8)]


def _build(with_udl):
    f = Frame2D()
    for nid, (x, y) in nodes.items():
        f.add_node(nid, x, y)
    f.add_section('sec', E=E, I=I, A=A)
    for mid, ni, nj in members:
        f.add_member(mid, node_i=ni, node_j=nj, section='sec')
    for n in [0, 1, 2]:
        f.fix(n)
    f.point_load(6, fx=10000)
    f.point_load(3, fx=15000)
    if with_udl:
        f.distributed_load(8, w=-12000)
        f.distributed_load(9, w=-12000)
        f.distributed_load(3, w=-10000)
        f.distributed_load(4, w=-10000)
    return f


def _check(f, result, sw_reactions, label):
    print(f"=== {label} ===")
    print(f"{'node':<6}{'Rx(FEM)':>12}{'Rx(SW)':>12}{'Ry(FEM)':>12}{'Ry(SW)':>12}{'M(FEM)':>12}{'M(SW)':>12}")
    for n in [0, 1, 2]:
        ux, uy, rot = f.dofs_of(n)
        Rx, Ry, M = result.reactions[ux] / 1000, result.reactions[uy] / 1000, result.reactions[rot] / 1000
        ref = sw_reactions[n]
        print(f"{n:<6}{Rx:>12.3f}{ref['Rx']:>12.3f}{Ry:>12.3f}{ref['Ry']:>12.3f}{M:>12.3f}{ref['M']:>12.3f}")
        for fem_val, ref_val, comp in [(Rx, ref['Rx'], 'Rx'), (Ry, ref['Ry'], 'Ry'), (M, ref['M'], 'M')]:
            assert abs(fem_val - ref_val) < 0.02, f"{label} node{n} {comp} 不吻合: FEM={fem_val}, SW={ref_val}"


# ---- 純點載重版本 (Case-08-two-story-two-bay_Report.pdf) ----
sw_no_udl = {
    0: dict(Rx=-8.058, Ry=-8.538, M=19.621),
    1: dict(Rx=-9.437, Ry=3.917, M=21.460),
    2: dict(Rx=-7.505, Ry=4.621, M=18.885),
}
f1 = _build(with_udl=False)
r1 = solve(f1)
_check(f1, r1, sw_no_udl, "純點載重版本")

# ---- 含均佈載重版本 (Case-08-two-story-two-bay-udl_Report.pdf) ----
sw_udl = {
    0: dict(Rx=-6.229, Ry=41.157, M=16.692),
    1: dict(Rx=-7.609, Ry=144.801, M=18.532),
    2: dict(Rx=-11.161, Ry=78.042, M=23.268),
}
f2 = _build(with_udl=True)
r2 = solve(f2)
_check(f2, r2, sw_udl, "含均佈載重版本")

print()
print("PASS: Case-08 兩個版本(純點載重+含UDL), 共18個反力分量, 全數跟SW FEA直接吻合")
print("      (不需要任何取負號/取絕對值的特例, 三方(frame2d/sd_framework/SW FEA)方向完全一致)")
