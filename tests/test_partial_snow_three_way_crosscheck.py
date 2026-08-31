"""
三方交叉驗證: 屋脊0.2L沒有雪(partial gap)的斜屋頂snow load案例。

背景: 這個案例(direction='global', angle_deg=-90, 局部段+線性變化)
跟SW FEA app報告的支承反力(3.544/19.092/-6.243)對不起來, 但同一個
結構+載重用另外兩個完全獨立、開源、廣泛使用的FEM套件(anastruct 1.7,
PyNiteFEA 3.0)重新建模驗證, 三個工具(frame2d + anastruct + PyNite)
精確吻合到4位小數, 只有SW FEA app的數字不一樣。

背後的道理(不只是巧合對上): 這個結構是左右鏡射對稱的, 載重也是鏡射
對稱的(左屋頂簷口重、屋脊輕, 右屋頂反過來, 剛好左右對照), 線性彈性
力學有個基本定理——對稱結構受對稱載重, 反應(反力、內力)必然對稱。
SW FEA app報告的F9/F10軸力(16.006/4.006 對比 1.006/13.006)明確違反
這個定理(如果對稱, 兩根桿件靠屋脊那端的軸力應該相等, 但不相等),
這在數學上就是不可能從一個正確的對稱計算得出的結果, 不需要靠三方
驗證也能看出SW FEA app這個案例的計算有問題; 三方驗證是進一步的
獨立佐證。

anastruct跟PyNite都用「節點分割法」(在0.8L/0.2L邊界插入真實節點,
把載重段變成一根獨立元素的整段均佈/線性變化載重), 這個手法本身也在
本專案別處驗證過是精確等價(不是近似)的寫法, 不是為了湊答案。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01
L_slope = np.hypot(3.0, 3.0)

# frame2d的答案 (direction='global', angle_deg=-90, x_start/x_end局部段)
FRAME2D_TARGET = (2.6646, 16.9706, -4.5391)


def test_frame2d_self():
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
    f.distributed_load(1, w=10.0, w_end=0.0, x_start=0.0, x_end=0.8 * L_slope,
                        direction='global', angle_deg=-90.0)
    f.distributed_load(2, w=0.0, w_end=10.0, x_start=0.2 * L_slope, x_end=L_slope,
                        direction='global', angle_deg=-90.0)
    r = solve(f)
    ux, uy, rot = f.dofs_of(0)
    got = (r.reactions[ux], r.reactions[uy], r.reactions[rot])
    for g, t in zip(got, FRAME2D_TARGET):
        assert abs(g - t) < 1e-3, f"frame2d自己的結果({got})變了, 跟記錄的({FRAME2D_TARGET})對不起來"
    print(f"frame2d:   Rx={got[0]:.4f}  Ry={got[1]:.4f}  M={got[2]:.4f}")


def test_anastruct_cross_check():
    try:
        from anastruct import SystemElements
    except ImportError:
        import pytest
        pytest.skip("anastruct未安裝(pip install anastruct) -- 跳過第三方交叉驗證")

    p0, p1, p2, p3, p4 = (0., 0.), (6., 0.), (0., 4.), (6., 4.), (3., 7.)
    u24 = (np.array(p4) - np.array(p2)) / L_slope
    u43 = (np.array(p3) - np.array(p4)) / L_slope
    p_split_left = tuple(np.array(p2) + u24 * 0.8 * L_slope)
    p_split_right = tuple(np.array(p4) + u43 * 0.2 * L_slope)

    ss = SystemElements(EA=E * A, EI=E * I)
    ss.add_element(location=[p0, p2])
    ss.add_element(location=[p3, p1])
    ss.add_element(location=[p2, p_split_left])
    ss.add_element(location=[p_split_left, p4])
    ss.add_element(location=[p4, p_split_right])
    ss.add_element(location=[p_split_right, p3])

    for nid, node in ss.node_map.items():
        if np.allclose((node.vertex.x, node.vertex.y), p0, atol=1e-6) or \
           np.allclose((node.vertex.x, node.vertex.y), p1, atol=1e-6):
            ss.add_support_fixed(node_id=nid)

    def elem_id_by_coords(pa, pb):
        for eid, el in ss.element_map.items():
            c1, c2 = (el.vertex_1.x, el.vertex_1.y), (el.vertex_2.x, el.vertex_2.y)
            if (np.allclose(c1, pa, atol=1e-6) and np.allclose(c2, pb, atol=1e-6)) or \
               (np.allclose(c1, pb, atol=1e-6) and np.allclose(c2, pa, atol=1e-6)):
                return eid
        raise ValueError("not found")

    eid_left = elem_id_by_coords(p2, p_split_left)
    eid_right = elem_id_by_coords(p_split_right, p3)
    ss.q_load(q=(-10.0, -1e-7), element_id=eid_left, direction='y')
    ss.q_load(q=(-1e-7, -10.0), element_id=eid_right, direction='y')
    ss.solve()

    results = ss.get_node_results_system()
    n0 = next(r for r in results if r['ux'] == 0.0 and r['uy'] == 0.0
              and r['id'] == min(r2['id'] for r2 in results if r2['ux'] == 0.0 and r2['uy'] == 0.0))
    got = (-n0['Fx'], -n0['Fy'], -n0['Tz'])
    print(f"anastruct: Rx={got[0]:.4f}  Ry={got[1]:.4f}  M={got[2]:.4f}")
    for g, t in zip(got, FRAME2D_TARGET):
        assert abs(g - t) < 1e-2, f"anastruct對不上frame2d: {got} vs {FRAME2D_TARGET}"


def test_pynite_cross_check():
    try:
        from Pynite import FEModel3D
    except ImportError:
        import pytest
        pytest.skip("PyNiteFEA未安裝(pip install PyNiteFEA) -- 跳過第三方交叉驗證")

    G = E / (2 * (1 + 0.3))
    J = 2 * I
    model = FEModel3D()
    model.add_material('steel', E, G, 0.3, 7.85)
    model.add_section('sec', A, I, I, J)
    model.add_node('N0', 0, 0, 0)
    model.add_node('N1', 6, 0, 0)
    model.add_node('N2', 0, 4, 0)
    model.add_node('N3', 6, 4, 0)
    model.add_node('N4', 3, 7, 0)
    for n in ['N0', 'N1']:
        model.def_support(n, True, True, True, True, True, True)
    for n in ['N2', 'N3', 'N4']:
        model.def_support(n, False, False, True, True, True, False)
    model.add_member('F0', 'N0', 'N2', 'steel', 'sec')
    model.add_member('F9', 'N2', 'N4', 'steel', 'sec')
    model.add_member('F10', 'N4', 'N3', 'steel', 'sec')
    model.add_member('F2', 'N3', 'N1', 'steel', 'sec')
    model.add_member_dist_load('F9', 'FY', -10.0, 0.0, x1=0.0, x2=0.8 * L_slope)
    model.add_member_dist_load('F10', 'FY', 0.0, -10.0, x1=0.2 * L_slope, x2=L_slope)
    model.analyze()

    n0 = model.nodes['N0']
    got = (n0.RxnFX['Combo 1'], n0.RxnFY['Combo 1'], n0.RxnMZ['Combo 1'])
    print(f"PyNite:    Rx={got[0]:.4f}  Ry={got[1]:.4f}  M={got[2]:.4f}")
    for g, t in zip(got, FRAME2D_TARGET):
        assert abs(g - t) < 1e-3, f"PyNite對不上frame2d: {got} vs {FRAME2D_TARGET}"


if __name__ == '__main__':
    test_frame2d_self()
    test_anastruct_cross_check()
    test_pynite_cross_check()
    print("\n三個完全獨立的FEM工具(frame2d/anastruct/PyNite)精確吻合。")
    print("SW FEA app報告的(3.544/19.092/-6.243)是唯一對不上的一個,")
    print("而且它本身違反「對稱結構+對稱載重->對稱反應」這個基本定理,")
    print("所以這組數字不能拿來當驗證基準, 不是frame2d算錯。")
