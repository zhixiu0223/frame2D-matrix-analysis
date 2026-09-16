"""
驗證frame2d.model.Frame2D.equal_dof()(跟OpenSeesPy的equalDOF同一個
概念)——用高勁度彈簧懲罰法實現, 見frame2d/dofmanager.py的
_apply_equal_dof()說明。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A, L = 200e9, 8e-5, 1e-2, 4.0


def two_cantilever_columns():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_node(2, 5, 0)
    f.add_node(3, 5, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec')
    f.add_member(1, node_i=2, node_j=3, section='sec')
    f.fix(0)
    f.fix(2)
    return f


print("=== 案例1: equal_dof(ux)讓兩根柱子頂端水平位移幾乎相等 ===")
f1 = two_cantilever_columns()
f1.equal_dof(1, 3, ux=True)
f1.point_load(1, fx=100e3, fy=0, m=0)
r1 = solve(f1)
u1x = r1.displacements[f1.dofs_of(1)[0]]
u3x = r1.displacements[f1.dofs_of(3)[0]]
rel_err = abs(u1x - u3x) / abs(u1x)
assert rel_err < 1e-6, f"equal_dof(ux)綁定的兩個自由度應該幾乎相等, 實際相對誤差={rel_err}"
print(f"PASS: u1x={u1x:.6f}, u3x={u3x:.6f}, 相對誤差={rel_err:.2e}\n")


print("=== 案例2: 物理正確性——兩柱並聯共同承擔側向力, 位移應該是單柱獨自承擔時的一半 ===")
f2 = two_cantilever_columns()
f2.point_load(1, fx=100e3, fy=0, m=0)
r2 = solve(f2)
u1x_alone = r2.displacements[f2.dofs_of(1)[0]]
u3x_alone = r2.displacements[f2.dofs_of(3)[0]]
F_theory = 100e3 * L**3 / (3 * E * I)
assert abs(u1x_alone - F_theory) / F_theory < 1e-9, "沒有equal_dof時應該精確對得上懸臂樑解析解"
assert abs(u3x_alone) < 1e-9, "柱2沒有equal_dof時應該完全不受柱1的力影響"
rel_err2 = abs(u1x - u1x_alone / 2) / (u1x_alone / 2)
assert rel_err2 < 1e-6, (
    f"兩柱equal_dof(ux)並聯共同承擔側向力, 位移應該是單柱獨自承擔時的一半"
    f"(兩根相同柱子並聯, 側向勁度加倍), 實際: 並聯u1x={u1x}, 單柱u1x_alone/2={u1x_alone/2}, "
    f"相對誤差={rel_err2}"
)
print(f"PASS: 單柱獨自承擔u1x_alone={u1x_alone:.6f}, 並聯後u1x={u1x:.6f}(≈一半), "
      f"相對誤差={rel_err2:.2e}\n")


print("=== 案例3: 沒有指定的方向(uy/rot)應該保持各自獨立, 不受equal_dof影響 ===")
f3 = two_cantilever_columns()
f3.equal_dof(1, 3, ux=True)   # 只綁ux, 不綁uy/rot
f3.point_load(1, fx=100e3, fy=-50e3, m=30e3)
f3.point_load(3, fx=0, fy=-200e3, m=0)   # 柱2另外單獨加一個完全不同的垂直力
r3 = solve(f3)
uy1 = r3.displacements[f3.dofs_of(1)[1]]
uy3 = r3.displacements[f3.dofs_of(3)[1]]
assert abs(uy1 - uy3) / max(abs(uy1), abs(uy3), 1e-30) > 0.5, (
    f"uy方向沒有被equal_dof綁定, 兩柱垂直位移應該各自獨立、明顯不同"
    f"(柱2額外加了完全不同的垂直力), 實際: uy1={uy1}, uy3={uy3}"
)
print(f"PASS: uy1={uy1:.6e}, uy3={uy3:.6e}(明顯不同, 沒有被錯誤綁定)\n")


print("=== 案例4: 至少要指定一個方向, 不然應該明確拒絕 ===")
f4 = two_cantilever_columns()
try:
    f4.equal_dof(1, 3)   # 沒有指定任何ux/uy/rot=True
    assert False, "沒有指定任何方向應該要raise ValueError"
except ValueError as e:
    assert 'equal_dof' in str(e)
print("PASS: 沒有指定任何方向正確拒絕\n")


print("=== 案例5: 沒有使用equal_dof時(預設空list), 完全不影響既有結果 ===")
f5 = two_cantilever_columns()
f5.point_load(1, fx=100e3, fy=0, m=0)
r5 = solve(f5)
assert np.array_equal(r5.displacements, r2.displacements), (
    "沒有equal_dof時(空list), 結果應該跟完全沒有這個功能時逐位元一致"
)
print("PASS: 空的equal_dofs列表對結果完全沒有影響\n")


print("=== 案例6: solve_pdelta()透過同一個_solve_once_dofmanager()核心,")
print("    equal_dof自動也支援, 不用另外改 ===")
from frame2d.dofmanager import solve_pdelta

f6 = two_cantilever_columns()
f6.equal_dof(1, 3, ux=True)
f6.point_load(1, fx=100e3, fy=-500e3, m=0)
r6 = solve_pdelta(f6)
u1x6 = r6.displacements[f6.dofs_of(1)[0]]
u3x6 = r6.displacements[f6.dofs_of(3)[0]]
rel_err6 = abs(u1x6 - u3x6) / abs(u1x6)
assert rel_err6 < 1e-6, f"solve_pdelta()也應該支援equal_dof, 實際相對誤差={rel_err6}"
print(f"PASS: solve_pdelta() u1x={u1x6:.6f}, u3x={u3x6:.6f}, 相對誤差={rel_err6:.2e}\n")

print("PASS: frame2d equal_dof所有案例通過(含solve_pdelta)")
