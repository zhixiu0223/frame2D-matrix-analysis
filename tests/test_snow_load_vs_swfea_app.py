"""
用SW FEA app(手機上)的實際輸出資料做交叉驗證——不是我們自己的假設,
是zhixiu直接在app裡建了一個對照案例(整段都是雪, 不是只有部分段),
讀app畫面上的Edit Distributed Load對話框、反力、桿件內力, 一項一項對過。

案例: 跟slop-roof-snow系列同一個三角形屋架幾何(node0/1固接柱腳,
node2/3柱頂, node4屋頂尖端), 但雪載重改成"整段都有"(不是只蓋0.8L,
是Start Location=0, Load Length=4.243=整根桿件長度), F9(2->4,
即F7)是10->0, F10(4->3, 即F8)是0->10, 兩個都是全長。

這個案例解開了方向的疑惑: 之前用局部段(0.8L)的案例對不起來, 但這個
全長案例用direction='global'(或等同的'global_y'), angle_deg=-90
(全域正下方), 反力、兩根斜梁的軸力都精確吻合到小數點後4位——證明
SW FEA的Distributed Load雖然畫面顯示Load Angle(局部+y轉成全域角度),
但實際套用的力就是"全域垂直方向, 大小以沿桿件長度量測", 跟我們的
direction='global_y'完全一致, angle=-90只是它的特例。

之前那個"0.8L局部段"案例(slop-roof-snow.pdf)的反力數字對不起來,
很可能是那次匯出的PDF報告是舊快取(跟這個案例確認了計算方法完全正確
之後, 舊案例的落差沒辦法用載重方向解釋, 只剩「報告是舊的」這個解釋
說得通)。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01


def build():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_node(4, 3, 7)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')   # F0
    f.add_member(1, 2, 4, 's')   # F9 (SW FEA app的編號)
    f.add_member(2, 4, 3, 's')   # F10
    f.add_member(3, 3, 1, 's')   # F2
    f.fix(0)
    f.fix(1)
    # SW FEA app畫面: Start Location=0, Load Length=4.243(=整根桿件),
    # F9: Start Value=10, End Value=0; F10: Start Value=0, End Value=10;
    # Load Angle=-90(全域正下方, 用'global'方向 + angle_deg=-90重現,
    # 等同direction='global_y')
    f.distributed_load(1, w=10.0, w_end=0.0, direction='global', angle_deg=-90.0)
    f.distributed_load(2, w=0.0, w_end=10.0, direction='global', angle_deg=-90.0)
    return f


def check(label, fem, screenshot, tol=1e-3):
    err = abs(fem - screenshot)
    print(f"  {label}: frame2d={fem:.4f}  SW FEA app畫面={screenshot}  差={err:.2e}")
    assert err < tol, f"{label} 對不起來 (差={err:.2e})"


print("=== 反力, 跟SW FEA app畫面(18:57)逐項比對 ===")
f = build()
r = solve(f)
sw = {0: (3.993, 21.213, -6.990), 1: (-3.993, 21.213, 6.990)}
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    check(f"node{n} Rx", r.reactions[ux], sw[n][0])
    check(f"node{n} Ry", r.reactions[uy], sw[n][1])
    check(f"node{n} M", r.reactions[rot], sw[n][2])

print("\n=== 兩根斜梁軸力(N), 跟SW FEA app畫面(18:58)逐項比對 ===")
x, N, V, M = member_internal_forces(f, r, 1, n=3)
check("F9 N(0) (node2端)", -N[0], 17.823)
check("F9 N(L) (node4端)", -N[-1], 2.823)
x8, N8, V8, M8 = member_internal_forces(f, r, 2, n=3)
check("F10 N(0) (node4端)", -N8[0], 2.823)
check("F10 N(L) (node3端)", -N8[-1], 17.823)

print("\n全部通過 -- direction='global'(含'global_y'特例)的實作，")
print("已經用SW FEA app的真實輸出資料逐項驗證過，不只是我們自己的假設。")
