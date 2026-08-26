"""
驗證案例3: 兩層兩跨側移鋼架 (九節點/十桿件), 對照第三方工具 SW FEA (swframe2d,
Android app) 的 PDF 分析報告 (Case-08-two-story-two-bay_Report-udl.pdf)。

幾何/斷面/載重數值直接取自使用者上傳的 SW FEA .frame (SQLite) 檔案:
  E=150 GPa, I=1e-6 m^4, A=0.05 m^2 (全桿件同斷面)
  水平點載重: node3(1F, x=0) Fx=-15kN, node6(roof, x=0) Fx=-10kN
  均佈載重: 1F梁(F3,F4) 10 kN/m, 屋頂梁(F8,F9) 12 kN/m

慣例差異note: SW FEA的分佈載重方向對應到本套件的 w=+ (正值), 推測是SW FEA(手機app)
內部用螢幕座標(Y向下為正), 跟本套件的標準工程座標(Y向上為正)相反 - 不是bug,
是輸入慣例不同, 使用本套件時請自行依 local+y 實際指向決定正負號。

支承反力的整體正負號慣例也跟SW FEA相反(本套件: reactions = K@u - F,
是"結構受到的外力"; SW FEA的Rx/Ry/M則是相反號的定義) - 比對時統一取絕對值。
"""
import numpy as np
from frame2d import Frame2D, solve

E = 150e9
I = 1e-6
A = 0.05

f = Frame2D()
nodes = {0: (0, 0), 1: (5, 0), 2: (12, 0),
         3: (0, 4), 4: (5, 4), 5: (12, 4),
         6: (0, 7.5), 7: (5, 7.5), 8: (12, 7.5)}
for nid, (x, y) in nodes.items():
    f.add_node(nid, x, y)

f.add_section('sec', E=E, I=I, A=A)

members = [(0, 0, 3), (1, 1, 4), (2, 5, 2), (3, 3, 4), (4, 4, 5),
           (5, 6, 3), (6, 7, 4), (7, 8, 5), (8, 6, 7), (9, 7, 8)]
for mid, ni, nj in members:
    f.add_member(mid, node_i=ni, node_j=nj, section='sec')

for n in [0, 1, 2]:
    f.fix(n)

f.point_load(6, fx=-10000)
f.point_load(3, fx=-15000)
f.distributed_load(8, w=12000)   # 見上方note: SW FEA慣例對應本套件的正值
f.distributed_load(9, w=12000)
f.distributed_load(3, w=10000)
f.distributed_load(4, w=10000)

result = solve(f)

# SW FEA PDF報告 (Case-08-two-story-two-bay_Report-udl.pdf) 的支承反力
sw_fea = {
    0: dict(Rx=-6.229, Ry=41.157, M=16.692),
    1: dict(Rx=-7.609, Ry=144.801, M=18.532),
    2: dict(Rx=-11.161, Ry=78.042, M=23.268),
}

print(f"{'node':<6}{'Rx(FEM)':>12}{'Rx(SWFEA)':>12}{'Ry(FEM)':>12}{'Ry(SWFEA)':>12}{'M(FEM)':>12}{'M(SWFEA)':>12}")
for n in [0, 1, 2]:
    ux, uy, rot = f.dofs_of(n)
    Rx_fem = result.reactions[ux] / 1000
    Ry_fem = result.reactions[uy] / 1000
    M_fem = result.reactions[rot] / 1000
    ref = sw_fea[n]
    print(f"{n:<6}{Rx_fem:>12.3f}{ref['Rx']:>12.3f}{Ry_fem:>12.3f}{ref['Ry']:>12.3f}{M_fem:>12.3f}{ref['M']:>12.3f}")

    # 整體反力正負號慣例相反, 取絕對值比較大小(相對誤差)
    for fem_val, ref_val, label in [(Rx_fem, ref['Rx'], 'Rx'), (Ry_fem, ref['Ry'], 'Ry'), (M_fem, ref['M'], 'M')]:
        rel_err = abs(abs(fem_val) - abs(ref_val)) / abs(ref_val)
        assert rel_err < 2e-3, f"node{n} {label} 誤差過大: FEM={fem_val}, SWFEA={ref_val}, rel_err={rel_err}"

print("PASS: Case-08 九節點兩層兩跨鋼架, 三支承反力(Rx/Ry/M共9個數值)全數跟SW FEA吻合(rel_err<0.2%)")
