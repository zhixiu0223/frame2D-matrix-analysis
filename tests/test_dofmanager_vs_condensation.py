"""
驗證案例17: DOFManager版本(frame2d/dofmanager.py) 跟主要求解器的
靜力凝縮(static condensation)版本交叉驗證

這兩套是完全獨立寫的程式碼(DOFManager版本不呼叫solve.py的任何組裝/求解
邏輯), 用來驗證"內部鉸接"這個功能的兩種不同實作方式(靜力凝縮 vs 真的
給每個release端獨立DOF)算出來的答案是否一致——理論上必須一致, 因為兩者
在數學上是等價的(高斯消去法, 只是消去時機不同)。

案例A/B: Gerber梁跟門型鋼架(跟test_element_release.py、
test_element_release_vs_swfea.py同一組模型), 用兩套實作分別求解比對。

案例C: 展示DOFManager版本的額外能力——桿件內部集中力可以直接加在有鉸接
的桿件上(靜力凝縮版本目前這個組合會直接報錯, 因為沒有推導release專屬
的固定端反力公式)。用「節點分割」的方式獨立驗證這個新案例的正確性
(在點載重位置直接放節點, 拆成兩段, 這樣就不用靠release專屬公式, 靜力
凝縮版本也能算, 兩邊再對一次)。
"""
import numpy as np
from frame2d import Frame2D, solve_condensation
from frame2d.dofmanager import solve_dofmanager

E, I, A = 200e6, 8e-5, 0.01


def check(label, a, b, atol=1e-6):
    err = abs(a - b)
    print(f"  {label}: 靜力凝縮={a:.6f}  DOFManager={b:.6f}  差={err:.2e}")
    assert err < atol, f"{label} 不吻合"


# ---- 案例A: Gerber梁 ----
print("=== 案例A: Gerber梁, 兩套實作交叉驗證 ===")
L1, L2, w = 6.0, 5.0, -10.0
fA = Frame2D()
fA.add_node(0, 0, 0)
fA.add_node(1, L1, 0)
fA.add_node(2, L1 + L2, 0)
fA.add_section('s', E=E, I=I, A=A)
fA.add_member(0, 0, 1, 's', release_j=True)
fA.add_member(1, 1, 2, 's')
fA.pin(0)
fA.roller_y(1)
fA.roller_y(2)
fA.distributed_load(0, w=w)

rA_static = solve_condensation(fA)
rA_dof = solve_dofmanager(fA)
for n in [0, 1, 2]:
    ux, uy, rot = fA.dofs_of(n)
    check(f"node{n} Ry", rA_static.reactions[uy], rA_dof.reactions[uy])
for mid in [0, 1]:
    diff = np.max(np.abs(rA_static.member_results[mid].end_forces_local - rA_dof.member_results[mid].end_forces_local))
    print(f"  member{mid} 端力最大差: {diff:.2e}")
    assert diff < 1e-6
print("PASS\n")


# ---- 案例B: 門型鋼架(跟test_element_release_vs_swfea.py同一個模型) ----
print("=== 案例B: 門型鋼架, 兩套實作交叉驗證 ===")
fB = Frame2D()
fB.add_node(0, 0, 0)
fB.add_node(1, 6, 0)
fB.add_node(2, 0, 4)
fB.add_node(3, 6, 4)
fB.add_section('s', E=E, I=I, A=A)
fB.add_member(0, 0, 2, 's')
fB.add_member(1, 2, 3, 's', release_j=True)
fB.add_member(2, 3, 1, 's')
fB.fix(0)
fB.fix(1)
fB.point_load(2, fx=10.0)

rB_static = solve_condensation(fB)
rB_dof = solve_dofmanager(fB)
for n in [0, 1]:
    ux, uy, rot = fB.dofs_of(n)
    check(f"node{n} Rx", rB_static.reactions[ux], rB_dof.reactions[ux])
    check(f"node{n} Ry", rB_static.reactions[uy], rB_dof.reactions[uy])
    check(f"node{n} M", rB_static.reactions[rot], rB_dof.reactions[rot])
max_diff = max(
    np.max(np.abs(rB_static.member_results[mid].end_forces_local - rB_dof.member_results[mid].end_forces_local))
    for mid in [0, 1, 2]
)
print(f"  三桿件端力最大差: {max_diff:.2e}")
assert max_diff < 1e-6
print("PASS\n")


# ---- 案例C: DOFManager的額外能力 -- 桿件內部集中力直接加在鉸接桿件上
#      (靜力凝縮版本這個組合會報錯), 用節點分割法獨立驗證正確性 ----
print("=== 案例C: DOFManager額外能力(release桿件上的桿件內部點載重) ===")
L, P, a = 8.0, 20.0, 3.0

fC = Frame2D()
fC.add_node(0, 0, 0)
fC.add_node(1, L, 0)
fC.add_section('s', E=E, I=I, A=A)
fC.add_member(0, 0, 1, 's', release_j=True)
fC.fix(0)
fC.pin(1)
fC.member_point_load(0, a=a, fy=-P)

# 確認靜力凝縮版本(參考實作)這個組合確實會報錯(還沒支援)
try:
    solve_condensation(fC)
    raise AssertionError("預期應該報錯(靜力凝縮版本還不支援這個組合), 但沒有報錯")
except ValueError:
    print("  確認: 靜力凝縮版本(參考實作)這個組合會報錯(還沒支援), 符合預期")

rC_dof = solve_dofmanager(fC)

# 獨立驗證: 節點分割版本(在a=3處放真實節點, 拆成兩段, 這樣不需要release
# 專屬固定端反力公式, 靜力凝縮版本也能算, 用來獨立核對DOFManager的新能力)
fC_split = Frame2D()
fC_split.add_node(0, 0, 0)
fC_split.add_node(1, a, 0)
fC_split.add_node(2, L, 0)
fC_split.add_section('s', E=E, I=I, A=A)
fC_split.add_member(0, 0, 1, 's')
fC_split.add_member(1, 1, 2, 's', release_j=True)
fC_split.fix(0)
fC_split.pin(2)
fC_split.point_load(1, fy=-P)
rC_split = solve_condensation(fC_split)

ux0, uy0, rot0 = fC.dofs_of(0)
check("node0 Ry (release桿件內部點載重 vs 節點分割獨立驗證)",
      rC_dof.reactions[uy0], rC_split.reactions[fC_split.dofs_of(0)[1]])
check("node0 M", rC_dof.reactions[rot0], rC_split.reactions[fC_split.dofs_of(0)[2]])
print("PASS: DOFManager能處理靜力凝縮版本目前還不支援的組合, 且答案通過獨立驗證\n")

# 額外確認: 現在的主要求解器 frame2d.solve (= solve_dofmanager) 也能直接算這個組合
from frame2d import solve as solve_production
rC_production = solve_production(fC)
check("node0 Ry (主要求解器solve() vs 節點分割獨立驗證)",
      rC_production.reactions[uy0], rC_split.reactions[fC_split.dofs_of(0)[1]])
print("PASS: 主要求解器frame2d.solve()現在直接支援這個組合, 不用另外呼叫solve_dofmanager\n")

print("PASS: DOFManager vs 靜力凝縮, 交叉驗證+額外能力展示, 全數完成")


# ---- 案例D: 不連續/不從0開始的node_id, 確認dofs_of()真的解耦了 ----
print()
print("=== 案例D: 不連續node_id(100,250,999), 兩套實作依然一致 ===")
fD = Frame2D()
fD.add_node(100, 0, 0)
fD.add_node(250, L1, 0)
fD.add_node(999, L1 + L2, 0)
fD.add_section('s', E=E, I=I, A=A)
fD.add_member(0, 100, 250, 's', release_j=True)
fD.add_member(1, 250, 999, 's')
fD.pin(100)
fD.roller_y(250)
fD.roller_y(999)
fD.distributed_load(0, w=w)
rD_static = solve_condensation(fD)
rD_dof = solve_dofmanager(fD)
assert fD.n_dof() == 9, f"3個節點應該只佔用9個DOF, 得到{fD.n_dof()}"
for nid in [100, 250, 999]:
    ux, uy, rot = fD.dofs_of(nid)
    check(f"node{nid} Ry(不連續id)", rD_static.reactions[uy], rD_dof.reactions[uy])
print("PASS: 不連續node_id下, DOF數量正確(9個, 不是300個), 兩套實作依然一致")
