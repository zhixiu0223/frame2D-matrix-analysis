"""
驗證案例5: slope_deflection_framework 的 Case-01 ~ Case-08 全部案例
(共11個: 01,02,03,04,04.5,05,06,06.5,07,07.5,08), 對照該repo自己的
sd_framework.py 本體直接算出的正確答案(不是讀notebook/PDF圖片轉錄,
是實際執行 SlopeDeflectionSolver._solve_core() 拿到的數字)。

幾何/支承/載重全部從 samples/model_*.py 的 draw_geometry() 讀出來,
節點座標、支承類型、載重施加位置都跟原始model一致。
"""
from frame2d import Frame2D, solve

_EI_DEFAULT = dict(E=1.0, I=15000.0, A=1e8)  # 全部案例統一斷面(對應各model的EI_numeric=15000)


def _check(label, fem, ref, tol=2e-2):
    rel = abs(abs(fem) - abs(ref)) / max(abs(ref), 1e-9)
    status = 'OK' if rel < tol else 'FAIL'
    print(f"  [{status}] {label:<8} FEM={fem:>10.4f}  ref={ref:>10.4f}  rel_err={rel:.2e}")
    assert rel < tol, f"{label} 誤差過大: FEM={fem}, ref={ref}"


def test_case01():
    print("=== Case-01: propped cantilever (L=6, w=20) ===")
    L, w = 6.0, 20.0
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', **_EI_DEFAULT)
    f.add_member(0, 0, 1, 's')
    f.fix(0)
    f.roller_y(1)
    f.distributed_load(0, w=-w)
    r = solve(f)
    _check('R_A', r.reactions[f.dofs_of(0)[1]], 75.0)
    _check('R_B', r.reactions[f.dofs_of(1)[1]], 45.0)
    _check('M_AB', -r.member_results[0].end_forces_local[2], -90.0)


def test_case02():
    print("=== Case-02: two-span beam (L1=5,L2=6, w1=15,w2=20) ===")
    L1, L2, w1, w2 = 5.0, 6.0, 15.0, 20.0
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L1, 0)
    f.add_node(2, L1 + L2, 0)
    f.add_section('s', **_EI_DEFAULT)
    f.add_member(0, 0, 1, 's')
    f.add_member(1, 1, 2, 's')
    f.fix(0)
    f.roller_y(1)
    f.roller_y(2)
    f.distributed_load(0, w=-w1)
    f.distributed_load(1, w=-w2)
    r = solve(f)
    _check('R_A', r.reactions[f.dofs_of(0)[1]], 26.654)
    _check('R_B', r.reactions[f.dofs_of(1)[1]], 119.58)
    _check('R_C', r.reactions[f.dofs_of(2)[1]], 48.766)
    _check('M_AB', -r.member_results[0].end_forces_local[2], -13.173)
    _check('M_BA', -r.member_results[0].end_forces_local[5], 67.404)


def test_case03():
    print("=== Case-03: no-sway frame (H=4,L=6, w=24) ===")
    H, L, w = 4.0, 6.0, 24.0
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, H)
    f.add_node(2, L, H)
    f.add_node(3, L, 0)
    f.add_section('s', **_EI_DEFAULT)
    f.add_member(0, 0, 1, 's')
    f.add_member(1, 1, 2, 's')
    f.add_member(2, 3, 2, 's')
    f.fix(0)
    f.fix(3)
    f.distributed_load(1, w=-w)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], 20.25)
    _check('H_D', r.reactions[f.dofs_of(3)[0]], -20.25)
    _check('M_AB', -r.member_results[0].end_forces_local[2], 27.0)
    _check('M_BA', -r.member_results[0].end_forces_local[5], 54.0)


def _build_sway_frame(H, L, P=0.0, w=0.0):
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, H)
    f.add_node(2, L, H)
    f.add_node(3, L, 0)
    f.add_section('s', **_EI_DEFAULT)
    f.add_member(0, 0, 1, 's')
    f.add_member(1, 1, 2, 's')
    f.add_member(2, 3, 2, 's')
    f.fix(0)
    f.fix(3)
    if P:
        f.point_load(1, fx=P)
    if w:
        f.distributed_load(1, w=-w)
    return f


def test_case04():
    print("=== Case-04: sway frame, no UDL (H=4,L=6,P=12) ===")
    f = _build_sway_frame(4.0, 6.0, P=12.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], -6.0)
    _check('H_D', r.reactions[f.dofs_of(3)[0]], -6.0)


def test_case04_5():
    print("=== Case-04.5: sway frame + UDL (H=4,L=6,P=12,w=24) ===")
    f = _build_sway_frame(4.0, 6.0, P=12.0, w=24.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], 14.25)
    _check('H_D', r.reactions[f.dofs_of(3)[0]], -26.25)
    ab = r.member_results[0].end_forces_local
    bc = r.member_results[1].end_forces_local
    dc = r.member_results[2].end_forces_local
    _check('M_AB', -ab[2], 12.6)
    _check('M_BA', -ab[5], 44.4)
    _check('M_BC', -bc[2], -44.4)
    _check('M_CB', -bc[5], 63.6)
    _check('M_DC', -dc[2], -41.4)
    _check('M_CD', -dc[5], -63.6)


def _build_two_story(H1, H2, L, w1=0.0, w2=0.0, P1=0.0, P2=0.0):
    Htot = H1 + H2
    f = Frame2D()
    f.add_node(0, 0, 0)       # A base-left
    f.add_node(1, 0, H1)      # B 1F-left
    f.add_node(2, 0, Htot)    # C roof-left
    f.add_node(3, L, Htot)    # D roof-right
    f.add_node(4, L, H1)      # E 1F-right
    f.add_node(5, L, 0)       # F base-right
    f.add_section('s', **_EI_DEFAULT)
    f.add_member(0, 0, 1, 's')  # AB
    f.add_member(1, 1, 2, 's')  # BC
    f.add_member(2, 2, 3, 's')  # CD (roof梁)
    f.add_member(3, 3, 4, 's')  # DE
    f.add_member(4, 4, 5, 's')  # EF
    f.add_member(5, 1, 4, 's')  # BE (1F梁)
    f.fix(0)
    f.fix(5)
    if w1:
        f.distributed_load(5, w=-w1)
    if w2:
        f.distributed_load(2, w=-w2)
    if P1:
        f.point_load(1, fx=P1)
    if P2:
        f.point_load(2, fx=P2)
    return f


def test_case05():
    print("=== Case-05: two-story frame, no lateral load (H1=4,H2=3.5,L=6,w1=24,w2=18) ===")
    f = _build_two_story(4.0, 3.5, 6.0, w1=24.0, w2=18.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], 8.497)
    _check('H_F', r.reactions[f.dofs_of(5)[0]], -8.497)
    _check('M_AB', -r.member_results[0].end_forces_local[2], 11.33)
    _check('M_BE', -r.member_results[5].end_forces_local[2], -64.447)


def test_case06():
    print("=== Case-06: two-story sway, no UDL (P1=15,P2=10) ===")
    f = _build_two_story(4.0, 3.5, 6.0, P1=15.0, P2=10.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], -12.5)
    _check('H_F', r.reactions[f.dofs_of(5)[0]], -12.5)
    _check('M_AB', -r.member_results[0].end_forces_local[2], -31.061)
    _check('M_BE', -r.member_results[5].end_forces_local[2], 24.245)


def test_case06_5():
    print("=== Case-06.5: two-story sway + UDL (P1=15,P2=10,w1=24,w2=18) ===")
    f = _build_two_story(4.0, 3.5, 6.0, w1=24.0, w2=18.0, P1=15.0, P2=10.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], -4.003)
    _check('H_F', r.reactions[f.dofs_of(5)[0]], -20.997)
    _check('M_AB', -r.member_results[0].end_forces_local[2], -19.732)
    _check('M_EB', -r.member_results[5].end_forces_local[5], 88.692)


def _build_two_bay(H, L1, L2, w1=0.0, w2=0.0, P=0.0):
    f = Frame2D()
    f.add_node(0, 0, 0)         # A base-left
    f.add_node(1, 0, H)         # B top-left
    f.add_node(2, L1, H)        # F top-mid
    f.add_node(3, L1, 0)        # E base-mid
    f.add_node(4, L1 + L2, H)   # H top-right
    f.add_node(5, L1 + L2, 0)   # G base-right
    f.add_section('s', **_EI_DEFAULT)
    f.add_member(0, 0, 1, 's')  # AB
    f.add_member(1, 3, 2, 's')  # EF
    f.add_member(2, 5, 4, 's')  # GH
    f.add_member(3, 1, 2, 's')  # BF
    f.add_member(4, 2, 4, 's')  # FH
    f.fix(0)
    f.fix(3)
    f.fix(5)
    if w1:
        f.distributed_load(3, w=-w1)
    if w2:
        f.distributed_load(4, w=-w2)
    if P:
        f.point_load(1, fx=P)
    return f


def test_case07():
    print("=== Case-07: two-bay frame, no lateral load (H=4,L1=5,L2=7,w1=20,w2=15) ===")
    f = _build_two_bay(4.0, 5.0, 7.0, w1=20.0, w2=15.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], 9.075)
    _check('H_E', r.reactions[f.dofs_of(3)[0]], 5.132)
    _check('H_G', r.reactions[f.dofs_of(5)[0]], -14.207)
    _check('M_AB', -r.member_results[0].end_forces_local[2], 13.371)
    _check('M_FB', -r.member_results[3].end_forces_local[5], 56.196)


def test_case07_5():
    print("=== Case-07.5: two-bay frame + P (H=4,L1=5,L2=7,w1=20,w2=15,P=10) ===")
    f = _build_two_bay(4.0, 5.0, 7.0, w1=20.0, w2=15.0, P=10.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], 5.915)
    _check('H_E', r.reactions[f.dofs_of(3)[0]], 1.107)
    _check('H_G', r.reactions[f.dofs_of(5)[0]], -17.022)
    _check('M_AB', -r.member_results[0].end_forces_local[2], 5.801)
    _check('M_GH', -r.member_results[2].end_forces_local[2], -24.781)
    _check('M_HG', -r.member_results[2].end_forces_local[5], -43.305)


def test_case08():
    print("=== Case-08: two-story two-bay, point loads only (H1=4,H2=3.5,L1=5,L2=7,P1=15,P2=10) ===")
    f = Frame2D()
    nodes = {0: (0, 0), 1: (5, 0), 2: (12, 0),
             3: (0, 4), 4: (5, 4), 5: (12, 4),
             6: (0, 7.5), 7: (5, 7.5), 8: (12, 7.5)}
    for nid, (x, y) in nodes.items():
        f.add_node(nid, x, y)
    f.add_section('s', **_EI_DEFAULT)
    members = [(0, 0, 3), (1, 1, 4), (2, 5, 2), (3, 3, 4), (4, 4, 5),
               (5, 6, 3), (6, 7, 4), (7, 8, 5), (8, 6, 7), (9, 7, 8)]
    for mid, ni, nj in members:
        f.add_member(mid, ni, nj, 's')
    for n in [0, 1, 2]:
        f.fix(n)
    f.point_load(6, fx=10.0)
    f.point_load(3, fx=15.0)
    r = solve(f)
    _check('H_A', r.reactions[f.dofs_of(0)[0]], -8.057)
    _check('H_E', r.reactions[f.dofs_of(1)[0]], -9.437)
    _check('H_I', r.reactions[f.dofs_of(2)[0]], -7.505)


if __name__ == '__main__':
    test_case01()
    test_case02()
    test_case03()
    test_case04()
    test_case04_5()
    test_case05()
    test_case06()
    test_case06_5()
    test_case07()
    test_case07_5()
    test_case08()
    print()
    print("PASS: slope_deflection_framework Case-01~08 (共11個案例) 全數跟 sd_framework.py 本體吻合")
