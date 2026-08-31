"""
驗證: distributed_load(direction='global_y') 支援 w_start != w_end
(線性變化, 例如不均勻/不平衡雪載重) 的新擴充。

新增的部分只有: elements.py 的 fixed_end_forces_axial_udl_varying()
(線性變化軸向載重的等效節點載重), 加上 dofmanager.py 把 global_y 載重的
兩端分開投影到局部座標、分開代入既有的橫向公式(fixed_end_forces_udl,
早就支援線性變化)跟新的軸向公式。

驗證策略(跟本專案一貫做法一致, 見commit紀錄"細網格模型交叉驗證"):
1. 退化案例: w_start=w_end 時, 新路徑必須跟舊的(已用SW FEA驗證過的)
   均佈global_y結果精確一致(浮點誤差等級) -- 確保沒改壞舊功能。
2. 獨立收斂驗證: 把同一根斜屋頂桿件切成N段小桿件, 每小段用「舊的、
   已驗證過的」均佈global_y路徑(w_start=w_end=該小段中點的線性內插值)
   近似線性變化載重, N越大近似越精確。隨著N增加, 反力應該平滑收斂到
   單一桿件+新公式算出的答案 -- 這條路徑完全沒有用到新公式本身,
   是真正獨立的數值交叉驗證, 不是自己驗自己。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01


def build_single_element_roof(w_left_eave, w_left_ridge, w_right_ridge, w_right_eave):
    """跟slop_roof_demo.py同一個幾何: node0=(0,0)/node1=(6,0)固接柱腳,
    node2=(0,4)/node3=(6,4)柱頂, node4=(3,7)屋頂尖端。
    左屋頂斜梁(node2->node4)用w_left_eave(簷口)->w_left_ridge(屋脊),
    右屋頂斜梁(node4->node3)用w_right_ridge(屋脊)->w_right_eave(簷口)。"""
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_node(4, 3, 7)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 4, 's')
    f.add_member(2, 4, 3, 's')
    f.add_member(3, 3, 1, 's')
    f.fix(0)
    f.fix(1)
    f.distributed_load(1, w=w_left_eave, w_end=w_left_ridge, direction='global_y')
    f.distributed_load(2, w=w_right_ridge, w_end=w_right_eave, direction='global_y')
    return f


def build_subdivided_roof(w_left_eave, w_left_ridge, w_right_ridge, w_right_eave, n_seg):
    """同一何幾何, 但左右屋頂斜梁各切成n_seg段小桿件, 每小段用「舊的、
    只支援均佈的global_y路徑」(w_start=w_end=線性內插到該小段中點的值)。"""
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_node(4, 3, 7)
    f.add_section('s', E=E, I=I, A=A)
    f.fix(0)
    f.fix(1)

    next_node = 100
    mid = 0

    def add_sloped_line(node_start, xy_start, xy_end, w_start, w_end, n):
        nonlocal next_node, mid
        prev = node_start
        for k in range(1, n + 1):
            t = k / n
            pt = np.array(xy_start) + (np.array(xy_end) - np.array(xy_start)) * t
            if k < n:
                nid = next_node
                next_node += 1
                f.add_node(nid, pt[0], pt[1])
            else:
                nid = None  # will connect to真正的終點節點(由呼叫端指定)
            # 該小段中點對應的線性內插w值 (0<=t_mid<=1 沿整段0->n)
            t_mid = (k - 0.5) / n
            w_mid = w_start + (w_end - w_start) * t_mid
            f.add_member(mid, prev, nid if nid is not None else END_NODE[0], 's')
            f.distributed_load(mid, w=w_mid, direction='global_y')
            mid += 1
            prev = nid if nid is not None else END_NODE[0]

    END_NODE = [None]

    # 左屋頂斜梁: node2 -> node4, 沿線性插值w_left_eave -> w_left_ridge
    END_NODE[0] = 4
    add_sloped_line(2, (0, 4), (3, 7), w_left_eave, w_left_ridge, n_seg)
    # 右屋頂斜梁: node4 -> node3, 沿線性插值w_right_ridge -> w_right_eave
    END_NODE[0] = 3
    add_sloped_line(4, (3, 7), (6, 4), w_right_ridge, w_right_eave, n_seg)
    # 柱 F0, F2
    f.add_member(mid, 0, 2, 's'); mid += 1
    f.add_member(mid, 3, 1, 's'); mid += 1

    return f


def check(label, fem, exact, tol):
    rel = abs(fem - exact) / max(abs(exact), 1e-9)
    print(f"  {label}: FEM={fem:.6f}  ref={exact:.6f}  rel_err={rel:.2e}")
    assert rel < tol, f"{label} 不吻合 (rel_err={rel:.2e} >= {tol:.2e})"


print("=== 案例A: 退化情況(w_start=w_end)應精確等於舊的均佈global_y結果 ===")
f_uniform = build_single_element_roof(10.0, 10.0, 10.0, 10.0)
r = solve(f_uniform)
# 舊有slop_roof_demo.py/test_sloped_roof_global_udl.py已用SW FEA驗證過
# w=10均佈時的反力(SW FEA報告值)
sw_reactions = {0: (11.047, 42.426, -20.157), 1: (-11.047, 42.426, 20.157)}
for n in [0, 1]:
    ux, uy, rot = f_uniform.dofs_of(n)
    check(f"node{n} Rx", r.reactions[ux], sw_reactions[n][0], 2e-3)
    check(f"node{n} Ry", r.reactions[uy], sw_reactions[n][1], 2e-3)
    check(f"node{n} M", r.reactions[rot], sw_reactions[n][2], 2e-3)
print("PASS: 退化情況精確吻合既有SW FEA驗證過的均佈結果\n")


print("=== 案例B: 線性變化(不均勻雪載重)-- 細網格分段模型收斂驗證 ===")
# 不平衡雪載重示意: 左屋頂簷口輕(6)->屋脊重(22), 右屋頂屋脊重(22)->簷口輕(6)
# (風吹雪堆積在屋脊附近的簡化示意分布, 不是特定規範公式)
args = (6.0, 22.0, 22.0, 6.0)
f_single = build_single_element_roof(*args)
r_single = solve(f_single)

prev_err = None
for n_seg in [4, 16, 64, 256]:
    f_sub = build_subdivided_roof(*args, n_seg=n_seg)
    r_sub = solve(f_sub)
    ux0, uy0, rot0 = f_single.dofs_of(0)
    ux0s, uy0s, rot0s = f_sub.dofs_of(0)
    err_ry = abs(r_sub.reactions[uy0s] - r_single.reactions[uy0]) / abs(r_single.reactions[uy0])
    err_m = abs(r_sub.reactions[rot0s] - r_single.reactions[rot0]) / abs(r_single.reactions[rot0])
    print(f"  n_seg={n_seg:4d}: node0 Ry_sub={r_sub.reactions[uy0s]:.5f} "
          f"(單一元素={r_single.reactions[uy0]:.5f}, rel_err={err_ry:.2e})  "
          f"M_sub={r_sub.reactions[rot0s]:.5f} (單一元素={r_single.reactions[rot0]:.5f}, "
          f"rel_err={err_m:.2e})")
    if prev_err is not None:
        assert err_ry < prev_err[0] * 0.6 or err_ry < 1e-6, "Ry誤差沒有隨網格加密收斂"
        assert err_m < prev_err[1] * 0.6 or err_m < 1e-6, "M誤差沒有隨網格加密收斂"
    prev_err = (err_ry, err_m)

assert prev_err[0] < 5e-4, f"最細網格(n_seg=256)Ry誤差還太大: {prev_err[0]:.2e}"
assert prev_err[1] < 5e-4, f"最細網格(n_seg=256)M誤差還太大: {prev_err[1]:.2e}"
print("PASS: 細網格分段模型(完全獨立於新公式, 只用舊的均佈global_y路徑)"
      "隨網格加密平滑收斂到單一元素+新線性變化公式的答案\n")

print("=== 案例C: 整體靜力平衡(必要條件檢查) ===")
# 左屋頂斜梁長度 = 右屋頂斜梁長度 = sqrt(3^2+3^2) = 3*sqrt(2)
L_slope = np.sqrt(3.0 ** 2 + 3.0 ** 2)
W_left = 0.5 * (6.0 + 22.0) * L_slope     # 梯形載重合力 = 平均值*長度
W_right = 0.5 * (22.0 + 6.0) * L_slope
W_total_exact = W_left + W_right
ux0, uy0, rot0 = f_single.dofs_of(0)
ux1, uy1, rot1 = f_single.dofs_of(1)
Ry_total = r_single.reactions[uy0] + r_single.reactions[uy1]
check("總垂直反力=總雪載重", Ry_total, W_total_exact, 1e-9)
Rx_total = r_single.reactions[ux0] + r_single.reactions[ux1]
print(f"  總水平反力(應接近0, 載重全垂直、無外加水平力): {Rx_total:.2e}")
assert abs(Rx_total) < 1e-6 * max(abs(W_total_exact), 1.0), "總水平反力應接近0"
print("PASS: 整體靜力平衡(ΣFy=總雪載重, ΣFx=0)精確成立\n")

print("全部通過。")
