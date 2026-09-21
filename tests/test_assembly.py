"""
驗證案例: frame2d.assembly.assemble_K() -- 動力分析D0

assemble_K()只是`pushover._assemble_stiffness_with_hinges()`的公開薄包裝
(見assembly.py說明), 沒有修改任何既有求解器。但這代表專案裡有「兩份
獨立寫的K組裝」(這裡呼叫的那份, 跟`dofmanager._solve_once_dofmanager()`
內嵌的那份)。這支測試守住的是「兩份不會悄悄分歧」, 分兩層:

層1(不依賴任何求解器): 剛體運動檢核
  沒有支承、沒有release的連通結構, K必須剛好有3個零特徵值(2個平移+1個
  轉動剛體模態), 而且K乘上剛體位移向量必須是零。這是結構力學最基本的
  性質, 完全不需要跟另一份組裝比較。

層2(跟既有求解器交叉比對): 用assemble_K()自己組K、自己劃分邊界條件、
  自己解(只含節點集中力), 對照 solve() / _solve_once_dofmanager()(帶
  P-Delta軸力) / solve_with_hinges() 的位移結果。涵蓋: 一般剛架、桿件
  release(額外DOF)、桁架(無轉角勁度的DOF)、equalDOF懲罰法、P-Delta
  幾何勁度、含塑鉸勁度。
"""
import numpy as np

from frame2d import Frame2D, solve
from frame2d.assembly import assemble_K
from frame2d.dofmanager import _solve_once_dofmanager, solve_with_hinges
from frame2d.hinge import HingeState

E, I, A = 200e6, 8e-5, 1e-2


def manual_static_solve(frame, asm):
    """用assemble_K()回傳的K, 自己做邊界條件劃分並求解(只含節點集中力,
    支承只接受0.0/None)。刻意不呼叫frame2d任何求解函式。"""
    n = asm.n_dof
    F = np.zeros(n)
    for pl in frame.point_loads:
        ux, uy, rot = frame.dofs_of(pl.node)
        F[ux] += pl.fx
        F[uy] += pl.fy
        F[rot] += pl.m
    fixed = np.zeros(n, dtype=bool)
    for s in frame.supports:
        for dof, val in zip(frame.dofs_of(s.node), (s.ux, s.uy, s.rot)):
            if val is not None:
                assert val == 0.0, "manual_static_solve只接受0.0/None的支承"
                fixed[dof] = True
    inactive = (~fixed) & np.isclose(np.diag(asm.K), 0.0)
    free = np.where((~fixed) & (~inactive))[0]
    u = np.zeros(n)
    u[free] = np.linalg.solve(asm.K[np.ix_(free, free)], F[free])
    return u[:asm.n_node_dof]


def max_rel_diff(a, b):
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


# ------------------------------------------------------------------
# 模型建構
# ------------------------------------------------------------------
def portal_frame():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 'sec').add_member(1, 1, 2, 'sec').add_member(2, 3, 2, 'sec')
    f.fix(0).fix(3)
    f.point_load(1, fx=12.0)
    f.point_load(2, fy=-20.0, m=3.0)
    return f


def portal_frame_with_release():
    f = portal_frame()
    f.members[1].release_j = True     # 梁右端鉸接(額外DOF)
    f.members[0].release_j = True     # 柱頂鉸接(第二個額外DOF)
    return f


def truss_triangle():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 4, 0).add_node(2, 2, 3)
    f.add_section('t', E=E, I=I, A=A)
    f.add_truss(0, 0, 1, 't').add_truss(1, 1, 2, 't').add_truss(2, 0, 2, 't')
    f.pin(0).roller_y(1)
    f.point_load(2, fx=5.0, fy=-10.0)
    return f


def portal_frame_equal_dof():
    f = portal_frame()
    f.equal_dof(master_node=1, slave_node=2, ux=True)   # 剛性樓板: 水平位移綁定
    return f


def free_free_frame():
    """沒有支承、沒有release: 連通的門型剛架, 用來檢查剛體模態。"""
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 'sec').add_member(1, 1, 2, 'sec').add_member(2, 3, 2, 'sec')
    return f


# ------------------------------------------------------------------
# 層1: 剛體運動檢核(不依賴任何求解器)
# ------------------------------------------------------------------
print("=== 層1: 無支承剛架的剛體運動檢核 ===")
f = free_free_frame()
asm = assemble_K(f)
K = asm.K
assert K.shape == (12, 12) and asm.n_node_dof == 12 and asm.n_extra_dof == 0
assert np.allclose(K, K.T, rtol=0, atol=1e-9 * np.max(np.abs(K))), "K不對稱"

eig = np.linalg.eigvalsh(K)
scale = np.max(np.abs(eig))
n_zero = int(np.sum(np.abs(eig) < 1e-9 * scale))
print(f"  零特徵值數量 = {n_zero} (預期3), 最小的4個特徵值/最大值 = {eig[:4] / scale}")
assert n_zero == 3, f"剛體模態數應為3, 實際{n_zero}"
assert eig[0] > -1e-9 * scale, "K不是半正定(有負特徵值)"

n = asm.n_node_dof
rx, ry, rr = np.zeros(n), np.zeros(n), np.zeros(n)
for nid, node in f.nodes.items():
    ux, uy, rot = f.dofs_of(nid)
    rx[ux] = 1.0
    ry[uy] = 1.0
    rr[ux], rr[uy], rr[rot] = -node.y, node.x, 1.0    # 繞原點的剛體轉動
for name, r in (("x平移", rx), ("y平移", ry), ("繞原點轉動", rr)):
    resid = np.max(np.abs(K @ r)) / np.max(np.abs(K))
    print(f"  K @ 剛體{name} 的最大殘差(相對) = {resid:.2e}")
    assert resid < 1e-12, f"剛體{name}不在K的零空間裡"

# ------------------------------------------------------------------
# 層2: 跟既有求解器交叉比對
# ------------------------------------------------------------------
print("=== 層2a: assemble_K自解 vs solve() ===")
cases = [
    ("一般門型剛架", portal_frame()),
    ("含2個release(額外DOF)", portal_frame_with_release()),
    ("桁架(無轉角勁度DOF)", truss_triangle()),
    ("equalDOF懲罰法", portal_frame_equal_dof()),
]
for name, frame in cases:
    a = assemble_K(frame)
    u_mine = manual_static_solve(frame, a)
    u_ref = solve(frame).displacements
    d = max_rel_diff(u_mine, u_ref)
    print(f"  {name}: n_dof={a.n_dof} (額外{a.n_extra_dof}), 位移最大相對差 = {d:.2e}")
    assert d < 1e-10, f"{name}: assemble_K自解的位移跟solve()不一致"

f = portal_frame_with_release()
assert assemble_K(f).n_extra_dof == 2 and assemble_K(f).K.shape == (14, 14), "release額外DOF數量不對"

print("=== 層2b: P-Delta幾何勁度(axial_forces) vs _solve_once_dofmanager(member_axial=) ===")
frame = portal_frame()
axial = {0: -80.0, 1: 5.0, 2: -60.0}     # 任意給定的軸力(拉力為正), 兩邊吃同一組
a = assemble_K(frame, axial_forces=axial)
u_mine = manual_static_solve(frame, a)
u_ref = _solve_once_dofmanager(frame, set(), member_axial=axial).displacements
d = max_rel_diff(u_mine, u_ref)
a0 = assemble_K(frame)
assert not np.allclose(a.K, a0.K), "軸力沒有改變K, 幾何勁度沒有疊加進去"
print(f"  位移最大相對差 = {d:.2e}")
assert d < 1e-10, "P-Delta: assemble_K自解的位移跟_solve_once_dofmanager不一致"

print("=== 層2c: 含塑鉸勁度(hinge_states) vs solve_with_hinges() ===")
frame = portal_frame()
hs = {0: HingeState(Mp1=50.0, Mp2=float('inf'), R_post_yield_1=500.0, R_post_yield_2=0.0),
      2: HingeState(Mp1=float('inf'), Mp2=40.0, R_post_yield_1=0.0, R_post_yield_2=300.0)}
hs[0].yielded[0] = True     # 柱0的i端已降伏, 用降伏後的軟化旋轉彈簧
hs[2].yielded[1] = True     # 柱2的j端已降伏
a = assemble_K(frame, hinge_states=hs)
u_mine = manual_static_solve(frame, a)
u_ref = solve_with_hinges(frame, hs).displacements
u_elastic = solve(frame).displacements
d = max_rel_diff(u_mine, u_ref)
print(f"  位移最大相對差 = {d:.2e}; 降伏後位移比純彈性大 {max_rel_diff(u_ref, u_elastic):.1%}")
assert d < 1e-10, "hinge: assemble_K自解的位移跟solve_with_hinges不一致"
assert max_rel_diff(u_ref, u_elastic) > 1e-3, "塑鉸沒有讓結構變軟, 測試案例本身沒有測到東西"

print("\n全部通過: assemble_K() 與既有求解器內的組裝一致, 且通過剛體運動檢核。")
