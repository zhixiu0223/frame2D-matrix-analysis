"""
驗證: point_load(F, angle_deg) 跟 member_point_load(direction='global',
F, angle_deg) -- 集中力的角度便利介面

背景: distributed_load已經有direction='global'(任意角度)的便利介面,
但point_load/member_point_load還要自己手動分解fx/fy, 這裡補齊同一套
模式:
  - point_load: fx/fy本來就是全域座標, 不用分解, 只是加F+angle_deg
    這個便利寫法(算出fx,fy後直接加總, 兩者等價)。
  - member_point_load: fx/fy是該桿件自己的局部座標(跟distributed_load
    的w一樣), direction='global'時用F+angle_deg指定全域方向, 依桿件
    角度自動分解成局部fx/fy(跟distributed_load(direction='global')
    同一套旋轉邏輯)。

驗證策略: 退化案例對照, 不是自己驗自己——手動算出局部分量(用跟
dofmanager.py同一套旋轉公式手算), 跟新介面的direction='global'自動
分解結果比較, 必須完全一致(不是"接近", 是精確到浮點誤差等級, 因為
理論上是同一個計算, 只是誰來做旋轉的差別)。
"""
import math
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01


def base_roof():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_node(4, 3, 7)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 4, 's')   # F7: 45度斜梁
    f.add_member(2, 4, 3, 's')
    f.add_member(3, 3, 1, 's')
    f.fix(0)
    f.fix(1)
    return f


print("=== 案例A: point_load(F, angle_deg) vs 手動填fx/fy, 應完全一致 ===")
F, ang = 10.0, -90.0
fA1 = base_roof()
fA1.point_load(4, fx=F * math.cos(math.radians(ang)), fy=F * math.sin(math.radians(ang)))
rA1 = solve(fA1)
fA2 = base_roof()
fA2.point_load(4, F=F, angle_deg=ang)
rA2 = solve(fA2)
diff = max(abs(rA1.reactions[i] - rA2.reactions[i]) for i in fA1.dofs_of(0))
print(f"  差異: {diff:.2e}")
assert diff < 1e-12
print("  PASS\n")


print("=== 案例B: member_point_load(direction='global') vs 手動分解局部fx/fy ===")
angle_F7 = math.radians(45)   # F7桿件自己的角度
F, a_glob = 10.0, math.radians(-90)   # 全域垂直向下
gx, gy = F * math.cos(a_glob), F * math.sin(a_glob)
c, s = math.cos(angle_F7), math.sin(angle_F7)
# 跟dofmanager.py同一套旋轉: local = [[c,s],[-s,c]] @ global
fx_local_manual = c * gx + s * gy
fy_local_manual = -s * gx + c * gy

fB1 = base_roof()
fB1.member_point_load(1, a=2.0, fx=fx_local_manual, fy=fy_local_manual)
rB1 = solve(fB1)
fB2 = base_roof()
fB2.member_point_load(1, a=2.0, direction='global', F=F, angle_deg=-90.0)
rB2 = solve(fB2)
diff = max(abs(rB1.reactions[i] - rB2.reactions[i]) for i in fB1.dofs_of(0))
print(f"  反力差異: {diff:.2e}")
assert diff < 1e-9

x1, N1, V1, M1 = member_internal_forces(fB1, rB1, 1, n=21)
x2, N2, V2, M2 = member_internal_forces(fB2, rB2, 1, n=21)
assert np.max(np.abs(N1 - N2)) < 1e-9
assert np.max(np.abs(V1 - V2)) < 1e-9
assert np.max(np.abs(M1 - M2)) < 1e-9
print("  N/V/M內力圖也完全一致")
print("  PASS\n")


print("=== 案例C: 錯誤處理 -- direction='global'缺angle_deg, 或跟fx/fy混用, 都要報錯 ===")
try:
    f = base_roof()
    f.member_point_load(1, a=2.0, direction='global', F=10.0)  # 缺angle_deg
    raise AssertionError("應該要報錯")
except ValueError:
    print("  缺angle_deg: 正確報錯")

try:
    f = base_roof()
    f.member_point_load(1, a=2.0, direction='global', F=10.0, angle_deg=-90, fx=5.0)
    raise AssertionError("應該要報錯")
except ValueError:
    print("  跟fx混用: 正確報錯")
print("  PASS\n")

print("PASS: 集中力角度便利介面驗證完成")
