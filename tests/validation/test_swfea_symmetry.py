"""
SW FEA validation matrix -- 分類4: 對稱性檢查(獨立於SW FEA的健檢工具)

線性彈性力學的基本定理: 對稱結構受對稱載重, 反應(反力、內力)必然
鏡射對稱。這條定理不需要任何外部工具驗證, 純粹是frame2d自己求解出
來的結果必須滿足的內部一致性檢查——如果frame2d對一個明確對稱的案例
算出不對稱的答案, 那一定是frame2d自己的bug, 不需要等SW FEA或其他
外部工具才能發現。

這個檢查在調查SW FEA的「局部段」問題時扮演了關鍵角色: 它是第一個
獨立於任何外部工具、純粹用力學定理就能斷定「SW FEA那組數字不可能是
正確計算結果」的證據(SW FEA報告的F9/F10軸力16.006/4.006 對比
1.006/13.006, 兩根鏡射桿件靠屋脊那端的值不相等, 直接違反這條定理)。

這裡把它寫成正式測試, 確保frame2d自己在對稱案例上永遠通過這個健檢。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01
L_slope = np.hypot(3.0, 3.0)


def _build_symmetric_gable_roof():
    """左右鏡射對稱的三角形屋架: 結構對稱(x=3為對稱軸), 載重也對稱
    (左屋頂簷口重、屋脊輕, 右屋頂鏡射一樣的形狀)。"""
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
    return f


def test_reactions_are_mirror_symmetric():
    """對稱結構+對稱載重 -> 兩側Ry必須相等, Rx/M必須正負相反(機器精度等級)"""
    f = _build_symmetric_gable_roof()
    r = solve(f)
    ux0, uy0, rot0 = f.dofs_of(0)
    ux1, uy1, rot1 = f.dofs_of(1)
    assert abs(r.reactions[uy0] - r.reactions[uy1]) < 1e-9
    assert abs(r.reactions[ux0] + r.reactions[ux1]) < 1e-9
    assert abs(r.reactions[rot0] + r.reactions[rot1]) < 1e-9


def test_member_internal_forces_are_mirror_symmetric():
    """兩根鏡射桿件(F9: 簷口->屋脊, F10: 屋脊->簷口), 對齊鏡射位置後
    軸力必須逐點相等(機器精度等級, 1e-12量級, 不是近似吻合)"""
    f = _build_symmetric_gable_roof()
    r = solve(f)
    x9, N9, V9, M9 = member_internal_forces(f, r, 1, n=50)
    x10, N10, V10, M10 = member_internal_forces(f, r, 2, n=50)
    L9 = x9[-1]
    for frac in np.linspace(0, 1, 11):
        n9_at = -np.interp(frac * L9, x9, N9)
        n10_mirror = -np.interp((1 - frac) * L9, x10, N10)
        assert abs(n9_at - n10_mirror) < 1e-9, (
            f"F9@{frac:.1f}L={n9_at} 跟 F10鏡射位置={n10_mirror} 不相等, "
            f"違反對稱定理 -- frame2d自己的bug, 優先修"
        )


if __name__ == '__main__':
    test_reactions_are_mirror_symmetric()
    test_member_internal_forces_are_mirror_symmetric()
    print("✅ frame2d對稱案例的反力/內力鏡射對稱, 精確到機器精度等級。")
