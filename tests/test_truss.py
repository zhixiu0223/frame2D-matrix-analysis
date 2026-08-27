"""
驗證案例6: 桁架(truss)元素

案例A: 單一桁架桿件, 軸向點載重, 解析解 (伸長量 = PL/EA, 軸力 = P)
案例B: 對稱雙桿桁架(A字撐架), 節點2受垂直向下載重P, 用「節點法」
       (method of joints) 手算靜力平衡獨立驗證軸力大小與正負號(受壓)
"""
import numpy as np
from frame2d import Frame2D, solve

# ---- 案例A: 單桿軸力 ----
print("=== 案例A: 單一桁架桿件 ===")
E, A, L, P = 200e6, 0.01, 5.0, 50.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=1.0, A=A)   # I對truss無意義, 隨便給
f.add_truss(0, 0, 1, 's')
f.pin(0)   # 桁架節點沒有轉角勁度, 用pin(鉸接符號)比fix(固定端符號)更符合物理意涵
f.point_load(1, fx=P)
r = solve(f)

elong_exact = P * L / (E * A)
elong_fem = r.displacements[f.dofs_of(1)[0]]
print(f"伸長量: FEM={elong_fem:.8e}  解析解={elong_exact:.8e}")
assert abs((elong_fem - elong_exact) / elong_exact) < 1e-10, "伸長量不符解析解"

N_fem = r.member_results[0].end_forces_local[3]   # Fx2, 應該=+P(拉力)
print(f"軸力(Fx2端): FEM={N_fem:.6f}  應為 +{P}(拉力)")
assert abs(N_fem - P) < 1e-8, "軸力不符"
print("PASS: 單桿軸力完全吻合解析解 (伸長量+軸力)\n")


# ---- 案例B: 對稱雙桿桁架(A字撐架), 節點法獨立驗證 ----
print("=== 案例B: 對稱雙桿桁架 (A-frame) ===")
# 節點0=(0,0)、節點1=(4,0) 為地面錨點, 節點2=(2,3) 承受垂直載重P(向下)
E, A = 200e6, 0.01
P = 10.0

f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, 4, 0)
f.add_node(2, 2, 3)
f.add_section('s', E=E, I=1.0, A=A)
f.add_truss(0, 0, 2, 's')   # 桿0: 節點0-節點2
f.add_truss(1, 1, 2, 's')   # 桿1: 節點1-節點2
f.pin(0)
f.pin(1)
f.point_load(2, fy=-P)
r = solve(f)

# 節點法手算(獨立於frame2d的計算路徑):
# 桿件方向都是 (2,3)<->(0,0)/(4,0), 長度 sqrt(2^2+3^2)=sqrt(13),
# 垂直方向的投影比例 = 3/sqrt(13)。對稱結構, 兩桿內力大小相等。
# 節點2垂直方向平衡: 2*|N|*(3/sqrt(13)) = P (只有垂直分量需要平衡, 水平分量因對稱互相抵銷)
# 結構像A字撐架撐住向下的載重, 直覺上兩根桿應該是"受壓"(壓力把節點2撐住不往下垮),
# 不是"受拉"(受拉的話反而會被拉向下方錨點, 沒有東西能提供向上的反力平衡外加載重)
L_bar = np.hypot(2, 3)
N_exact_magnitude = P * L_bar / (2 * 3)   # = P*sqrt(13)/6
N_exact = -N_exact_magnitude   # 受壓, 取負號(tension-positive慣例下應為負)

for mid in [0, 1]:
    N_fem = r.member_results[mid].end_forces_local[3]
    print(f"桿{mid}軸力: FEM={N_fem:.6f}  節點法手算={N_exact:.6f}")
    assert abs(N_fem - N_exact) < 1e-6, f"桿{mid}軸力不符節點法手算"

print("PASS: 對稱雙桿桁架軸力(大小+受壓方向)完全吻合節點法獨立手算\n")

# ---- 額外鎖住 postprocess.member_internal_forces() 的N(x)正負號慣例 ----
# (這裡曾經抓到一個真的bug: N原本用+Fx1, 剛好是"壓力為正", 跟工程慣例相反,
#  導致這個A字撐架的壓力桿在軸力圖上顯示成正值、看起來像在受拉, 已修正成-Fx1)
from frame2d.postprocess import member_internal_forces
x, N, V, M = member_internal_forces(f, r, 0, n=3)
assert np.all(N < 0), f"A字撐架應該是壓力(N<0), 得到 {N}"
assert abs(N[0] - N_exact) < 1e-6, f"N(x)大小應為{N_exact}, 得到{N[0]}"
print(f"PASS: postprocess軸力圖 N(x)={N[0]:.4f} 拉力為正慣例正確(壓力顯示為負)\n")

print("PASS: truss元素(軸力桁架)驗證完成")
