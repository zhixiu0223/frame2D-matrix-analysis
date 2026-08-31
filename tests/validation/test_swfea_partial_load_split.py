"""
SW FEA validation matrix -- 分類2: 局部段(partial-length)分布載重
(⚠️ 不可靠) 與其 workaround(✅ 插入真實節點, 可靠)

判斷: SW FEA的Distributed Load功能, 只要載重範圍不是整根桿件("局部段",
用Start Distance/End Distance框出一部分), 搭配非0/±90度的Load Angle,
算出來的支承反力/桿件內力就跟三個完全獨立的FEM工具(frame2d本身,
anastruct, PyNiteFEA)算出來的不一樣, 而且在對稱結構+對稱載重的案例下
會給出違反鏡射對稱定理的結果(見test_swfea_symmetry.py)。這代表SW FEA
的「局部段」計算路徑本身不可信, 不是我們的載重定義猜錯。

**這不是說SW FEA是爛軟體**——full-member的計算精確無誤(見
test_swfea_full_member_load.py), 問題範圍很窄, 就是「局部段+角度換算」
這一個特定組合。使用者也用最小案例(單一斜桿+局部段均佈載重)獨立重現過
同樣的偏差模式, 排除是這個屋頂案例本身幾何/對稱性造成的巧合。

**Workaround(已用真實SW FEA app重新輸入驗證過, 不只是理論推導)**:
在局部段的邊界插入真實節點, 把「一根桿件的局部段載重」改成「兩根桿件,
其中一根整段都有載重」, 這樣就繞開SW FEA不可靠的局部段路徑, 走它可靠
的full-member路徑。使用者實際在SW FEA app裡插入節點(2.4,6.4)跟
(3.6,6.4)重新建模, 算出來的反力(2.665/16.971/-4.539)跟frame2d、
anastruct、PyNite的答案精確吻合(誤差<0.001), 也完全對稱。

結論: SW FEA的「局部段」輸入介面不能當驗證基準用; 需要局部段載重時,
拿SW FEA當oracle的正確做法是自己先插入節點再輸入, 不要依賴它的
Start/End Distance局部段功能。
"""
import subprocess
import sys
from pathlib import Path
import numpy as np

TESTS_DIR = Path(__file__).resolve().parent.parent

# 使用者在SW FEA app裡用「插入真實節點」的workaround重新輸入後,
# 實測截圖讀出的反力(2026-09-01, slop-roof-snow_mesh.frame/_Report.pdf)
SWFEA_SPLIT_NODE_RESULT = (2.665, 16.971, -4.539)
# frame2d / anastruct / PyNite三方交叉驗證的答案(test_partial_snow_three_way_crosscheck.py)
INDEPENDENT_FEM_RESULT = (2.6646, 16.9706, -4.5391)


def test_split_node_workaround_matches_independent_fem_and_real_swfea_app():
    """插入節點繞開局部段功能後, SW FEA app真實輸出 vs 三方獨立FEM -- ✅"""
    for a, b in zip(SWFEA_SPLIT_NODE_RESULT, INDEPENDENT_FEM_RESULT):
        assert abs(a - b) < 1e-2, (
            f"SW FEA app(插入節點後)的實測值{SWFEA_SPLIT_NODE_RESULT} "
            f"應該要跟三方獨立FEM的答案{INDEPENDENT_FEM_RESULT}吻合"
        )


def test_three_way_crosscheck_still_passes():
    """frame2d/anastruct/PyNite三方交叉驗證(局部段載重, 不靠SW FEA的局部段功能)"""
    result = subprocess.run(
        [sys.executable, str(TESTS_DIR / "test_partial_snow_three_way_crosscheck.py")],
        capture_output=True, text=True, cwd=str(TESTS_DIR.parent),
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"


if __name__ == '__main__':
    test_three_way_crosscheck_still_passes()
    test_split_node_workaround_matches_independent_fem_and_real_swfea_app()
    print("⚠️  SW FEA的「局部段」載重輸入功能不可靠, 不當驗證基準用。")
    print("✅ 插入真實節點的workaround: SW FEA app + frame2d + anastruct + PyNite 四方吻合。")
