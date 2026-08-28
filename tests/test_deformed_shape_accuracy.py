"""
驗證案例12: member_deformed_shape() 的撓度曲線精確度

背景: 舊版用單一桿件兩端內插一條Hermite三次多項式畫變形曲線, 桿件內部有
集中力/局部段均佈載重時, 因為真實撓度曲線在那個位置有更高階的不連續,
單一三次多項式描述不出來, 造成明顯低估(實測case: 低估約23%, 位置也偏)。
新版做法: 直接對已經驗證過的 M(x) 除以EI做兩次數值積分, 用兩端FEM算出的
真實節點位移釘住積分常數, 求解器本身的精度直接帶到畫圖這一步, 不需要
另外維護一個「近似畫圖版」+「精確驗算版」。

驗證分兩層:
案例A: 簡支梁均佈載重(已有解析解 5wL⁴/384EI), 確認新方法在均佈載重情況
       依然精確(不會因為改寫而退步)。
案例B: 簡支梁內部集中力(重現開發過程中發現的舊方法23%低估案例), 對照
       「直接在載重點插入節點求解」的精確節點解, 新方法應該吻合。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_deformed_shape

E, I, A = 200e6, 8e-5, 1e-2


def check(label, fem, exact, tol=1e-4):
    rel = abs(fem - exact) / abs(exact)
    print(f"  {label}: FEM={fem:.6f}  解析解={exact:.6f}  rel_err={rel:.2e}")
    assert rel < tol, f"{label} 不吻合"


# ---- 案例A: 簡支梁均佈載重, 新方法應該依然精確吻合解析解 ----
print("=== 案例A: 簡支梁均佈載重, 新的積分法跨中撓度對照解析解 ===")
L, w = 6.0, -5.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.distributed_load(0, w=w)
r = solve(f)

X, Y = member_deformed_shape(f, r, 0, scale=1.0, n=201)
i_mid = np.argmin(np.abs(X - L / 2))
v_mid_exact = 5 * w * L**4 / (384 * E * I)
check("跨中撓度", Y[i_mid], v_mid_exact, tol=1e-3)
print("PASS: 均佈載重情況新方法依然精確吻合解析解\n")


# ---- 案例B: 簡支梁內部集中力(重現舊方法23%低估的case), 對照節點分割精確解 ----
print("=== 案例B: 簡支梁內部集中力, 新方法對照節點分割精確解 ===")
L2, P, a = 10.0, 30.0, 3.5
f2 = Frame2D()
f2.add_node(0, 0, 0)
f2.add_node(1, L2, 0)
f2.add_section('s', E=E, I=I, A=A)
f2.add_member(0, 0, 1, 's')
f2.pin(0)
f2.roller_y(1)
f2.member_point_load(0, a=a, fy=-P)
r2 = solve(f2)
X2, Y2 = member_deformed_shape(f2, r2, 0, scale=1.0, n=201)
i_max = np.argmax(np.abs(Y2))
v_max_new = Y2[i_max]
x_max_new = X2[i_max]
print(f"  新方法: 最大撓度={v_max_new*1000:.4f}mm  位置x={x_max_new:.3f}")

# 節點分割精確解(在最大值附近x=4.5直接放節點求解, 該點位移是FEM節點解無內插誤差)
x_check = 4.5
f2_exact = Frame2D()
f2_exact.add_node(0, 0, 0)
f2_exact.add_node(1, x_check, 0)
f2_exact.add_node(2, L2, 0)
f2_exact.add_section('s', E=E, I=I, A=A)
f2_exact.add_member(0, 0, 1, 's')
f2_exact.add_member(1, 1, 2, 's')
f2_exact.pin(0)
f2_exact.roller_y(2)
f2_exact.member_point_load(0, a=a, fy=-P)   # 載重仍在a=3.5(member0內部), 不在新節點上
r2_exact = solve(f2_exact)
v_exact_at_x = r2_exact.displacements[f2_exact.dofs_of(1)[1]]
print(f"  節點分割精確解(x={x_check}): {v_exact_at_x*1000:.4f}mm")

# 用新方法在同一個x=4.5位置取值來比較(不是比較各自的argmax, 避免取樣網格差異)
i_45 = np.argmin(np.abs(X2 - x_check))
check(f"x={x_check}處撓度(新方法 vs 節點分割精確解)", Y2[i_45], v_exact_at_x, tol=1e-3)
print(f"  (對照: 舊Hermite內插法在這個case只算出約-26.73mm, 明顯低估)")
print("PASS: 桿件內部集中力情況, 新方法精確吻合節點分割精確解, 不再低估\n")


# ---- 案例C: 局部段均佈載重, 對照200段細網格模型(完全獨立的第二種算法) ----
print("=== 案例C: 局部段均佈載重, 積分法對照200段細網格模型 ===")
L3, w3, c3, d3 = 12.0, -20.0, 3.0, 7.0
f3 = Frame2D()
f3.add_node(0, 0, 0)
f3.add_node(1, L3, 0)
f3.add_section('s', E=E, I=I, A=A)
f3.add_member(0, 0, 1, 's')
f3.pin(0)
f3.roller_y(1)
f3.distributed_load(0, w=w3, x_start=c3, x_end=d3)
r3 = solve(f3)
X3, Y3 = member_deformed_shape(f3, r3, 0, scale=1.0, n=1001)
i_max3 = np.argmax(np.abs(Y3))
v_max_integration = Y3[i_max3]

# 200段細網格(每段各自承受對應比例的均佈載重, 節點解本身不需要任何內插)
n_seg = 200
f3b = Frame2D()
xs = np.linspace(0, L3, n_seg + 1)
for i, x in enumerate(xs):
    f3b.add_node(i, x, 0)
f3b.add_section('s', E=E, I=I, A=A)
for i in range(n_seg):
    f3b.add_member(i, i, i + 1, 's')
f3b.pin(0)
f3b.roller_y(n_seg)
for i in range(n_seg):
    x0, x1 = xs[i], xs[i + 1]
    lo, hi = max(x0, c3), min(x1, d3)
    if hi > lo:
        f3b.distributed_load(i, w=w3, x_start=lo - x0, x_end=hi - x0)
r3b = solve(f3b)
disp3b = np.array([r3b.displacements[f3b.dofs_of(i)[1]] for i in range(n_seg + 1)])
i_max3b = np.argmax(np.abs(disp3b))
v_max_mesh = disp3b[i_max3b]

print(f"  積分法(單一桿件): 最大撓度={v_max_integration*1000:.4f}mm  位置x={X3[i_max3]:.3f}")
print(f"  200段細網格(節點解): 最大撓度={v_max_mesh*1000:.4f}mm  位置x={xs[i_max3b]:.3f}")
rel_diff = abs(v_max_integration - v_max_mesh) / abs(v_max_mesh)
assert rel_diff < 1e-3, f"積分法跟細網格模型誤差過大: {rel_diff:.2e}"
print(f"  相對誤差: {rel_diff:.2e}")
print("PASS: 局部段均佈載重情況, 積分法精確吻合完全獨立的細網格模型\n")

print("PASS: member_deformed_shape() 撓度曲線精確度驗證完成")