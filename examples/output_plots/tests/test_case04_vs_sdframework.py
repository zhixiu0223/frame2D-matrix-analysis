"""
驗證案例4: 側移單跨剛架 Case-04 / Case-04.5, 對照使用者自己的
slope_deflection_framework repo 的 notebook 結果
(Case-04-sway-frame-wget.ipynb / Case-04_5-sway-frame-with-udl-wget.ipynb),
該repo內部已經用 sd_framework(傾角變位法) 跟 anastruct 互相交叉驗證過。

真實參數(從notebook挖出來, 不是隨手編的示範值):
  H=4.0, L=6.0, P=12.0 (B點水平力), w=24.0 kN/m (僅Case-04.5, 梁上均佈載重向下)
  EI=15000.0, EA=1e8 (anastruct: SystemElements(EA=1e8, EI=15000.0))

慣例note: frame2d的整體反力/彎矩正負號跟sd_framework相反(統一取絕對值比較),
但分佈載重方向不用像SW FEA那次一樣轉換 —— 這裡向下直接用負值(w=-24)就對了,
因為anastruct/sd_framework跟frame2d都是標準工程座標(Y向上為正), 不像SW FEA
(Android app)用螢幕座標。
"""
import numpy as np
from frame2d import Frame2D, solve


def build_case04(w=0.0):
    H, L, P = 4.0, 6.0, 12.0
    f = Frame2D()
    f.add_node(0, 0, 0)   # A (base)
    f.add_node(1, 0, H)   # B (top-left)
    f.add_node(2, L, H)   # C (top-right)
    f.add_node(3, L, 0)   # D (base)
    f.add_section('sec', E=1.0, I=15000.0, A=1e8)
    f.add_member(0, node_i=0, node_j=1, section='sec')  # AB
    f.add_member(1, node_i=1, node_j=2, section='sec')  # BC (梁)
    f.add_member(2, node_i=3, node_j=2, section='sec')  # DC
    f.fix(0)
    f.fix(3)
    f.point_load(1, fx=P)
    if w != 0.0:
        f.distributed_load(1, w=w)
    return f


def check(label, fem_val, ref_val, tol=2e-2):
    rel_err = abs(abs(fem_val) - abs(ref_val)) / max(abs(ref_val), 1e-9)
    print(f"{label:<8} FEM={fem_val:>10.4f}  ref={ref_val:>10.4f}  rel_err={rel_err:.2e}")
    assert rel_err < tol, f"{label} 誤差過大"


print("=== Case-04 (無UDL) ===")
f = build_case04(w=0.0)
r = solve(f)
Rx_A = r.reactions[f.dofs_of(0)[0]]
Rx_D = r.reactions[f.dofs_of(3)[0]]
check("H_A", Rx_A, -6.0)
check("H_D", Rx_D, -6.0)
print("PASS: Case-04 水平反力吻合 notebook (H_A=H_D=6.0, 對稱斷面下P對半分)\n")

print("=== Case-04.5 (含梁UDL w=24 kN/m 向下) ===")
f = build_case04(w=-24.0)
r = solve(f)

# notebook反力交叉檢查值 (由端彎矩推算): H_A=14.25, H_D=-26.25
Rx_A = r.reactions[f.dofs_of(0)[0]]
Rx_D = r.reactions[f.dofs_of(3)[0]]
check("H_A", Rx_A, 14.25)
check("H_D", Rx_D, -26.25)

# notebook 六個桿端彎矩交叉檢查值 (cell 11, 已通過該repo自己的anastruct驗證)
ab = r.member_results[0].end_forces_local  # AB: [.., .., M_A, .., .., M_B]
bc = r.member_results[1].end_forces_local  # BC: [.., .., M_B, .., .., M_C]
dc = r.member_results[2].end_forces_local  # DC: [.., .., M_D, .., .., M_C]

check("M_AB", ab[2], 12.6)
check("M_BA", ab[5], 44.4)
check("M_BC", bc[2], -44.4)
check("M_CB", bc[5], 63.6)
check("M_DC", dc[2], -41.4)
check("M_CD", dc[5], -63.6)

print("PASS: Case-04.5 水平反力+六個桿端彎矩, 全數跟 sd_framework/anastruct 吻合")
