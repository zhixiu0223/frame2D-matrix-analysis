"""
驗證: distributed_load(direction='global') -- 全域任意角度均佈載重
(global_y的推廣版, 支援任意角度+局部段+線性變化的任意組合)。

驗證策略: 退化案例對照, 不是自己驗自己。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01


def build_roof(direction_kwargs_left, direction_kwargs_right):
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
    f.distributed_load(1, **direction_kwargs_left)
    f.distributed_load(2, **direction_kwargs_right)
    return f


def check(label, fem, exact, tol=1e-8):
    rel = abs(fem - exact) / max(abs(exact), 1e-9)
    print(f"  {label}: FEM={fem:.6f}  ref={exact:.6f}  rel_err={rel:.2e}")
    assert rel < tol, f"{label} 不吻合"


print("=== 案例A: angle_deg=-90(全域直下) 應精確等於既有direction='global_y' ===")
f_gy = build_roof(dict(w=10.0, direction='global_y'),
                   dict(w=10.0, direction='global_y'))
f_ga = build_roof(dict(w=10.0, direction='global', angle_deg=-90.0),
                   dict(w=10.0, direction='global', angle_deg=-90.0))
r_gy = solve(f_gy)
r_ga = solve(f_ga)
for n in [0, 1]:
    ux, uy, rot = f_gy.dofs_of(n)
    check(f"node{n} Rx", r_ga.reactions[ux], r_gy.reactions[ux])
    check(f"node{n} Ry", r_ga.reactions[uy], r_gy.reactions[uy])
    check(f"node{n} M", r_ga.reactions[rot], r_gy.reactions[rot])
print("PASS\n")

print("=== 案例B: 局部段+angle_deg=-90 應精確等於'local'方向在水平桿件上的結果 ===")
L = 12.0
f_local = Frame2D()
f_local.add_node(0, 0, 0)
f_local.add_node(1, L, 0)
f_local.add_section('s', E=E, I=I, A=A)
f_local.add_member(0, 0, 1, 's')
f_local.pin(0)
f_local.roller_y(1)
f_local.distributed_load(0, w=-20.0, x_start=3.0, x_end=7.0)
r_local = solve(f_local)

f_global = Frame2D()
f_global.add_node(0, 0, 0)
f_global.add_node(1, L, 0)
f_global.add_section('s', E=E, I=I, A=A)
f_global.add_member(0, 0, 1, 's')
f_global.pin(0)
f_global.roller_y(1)
f_global.distributed_load(0, w=20.0, x_start=3.0, x_end=7.0,
                           direction='global', angle_deg=-90.0)
r_global = solve(f_global)

ux0, uy0, _ = f_local.dofs_of(0)
ux1, uy1, _ = f_local.dofs_of(1)
check("node0 Ry", r_global.reactions[uy0], r_local.reactions[uy0])
check("node1 Ry", r_global.reactions[uy1], r_local.reactions[uy1])
print("PASS\n")

print("全部通過。")
