"""
SW FEA validation matrix -- 分類1: 整根桿件的分布載重 (✅ 可靠)

判斷: SW FEA對「整根桿件都受載重」(不管均佈或線性變化, 局部方向或
全域方向)算得跟frame2d精確吻合, 可以當作驗證基準用。

這個檔案不重複計算, 只是把散落在其他測試檔案裡的關鍵結論集中引用,
當作「驗證矩陣」的索引入口——實際的詳細推導、原始SW FEA截圖數字比對,
在各自的原始測試檔案裡(這裡只斷言最終結果, 不重講一次過程)。

證據來源:
  - test_sloped_roof_global_udl.py: 全域垂直均佈載重(10kN/m, 整根桿件),
    對照SW FEA的斜屋頂案例反力(11.047/42.426/-20.157), 精確吻合。
  - test_snow_load_vs_swfea_app.py: 全域垂直線性變化載重(梯形10->0,
    整根桿件, 不是局部段), 對照SW FEA app真實截圖+.frame資料庫匯出的
    反力(3.993/21.213/-6.990)跟兩根斜梁的軸力(17.823/2.823), 精確吻合
    到4位小數。
"""
import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent.parent


def _run(test_file):
    result = subprocess.run(
        [sys.executable, str(TESTS_DIR / test_file)],
        capture_output=True, text=True, cwd=str(TESTS_DIR.parent),
    )
    assert result.returncode == 0, (
        f"{test_file} 失敗 (應該是✅可靠分類, 如果這裡開始失敗代表"
        f"full-member load的驗證基準本身被改壞了, 要優先處理):\n"
        f"{result.stdout}\n{result.stderr}"
    )
    return result.stdout


def test_full_member_uniform_load_matches_swfea():
    """全域垂直均佈載重(整根桿件) vs SW FEA -- ✅ 可靠"""
    _run("test_sloped_roof_global_udl.py")


def test_full_member_trapezoidal_load_matches_swfea_app():
    """全域垂直線性變化載重(梯形, 整根桿件) vs SW FEA app真實輸出 -- ✅ 可靠"""
    _run("test_snow_load_vs_swfea_app.py")


if __name__ == '__main__':
    test_full_member_uniform_load_matches_swfea()
    test_full_member_trapezoidal_load_matches_swfea_app()
    print("✅ 整根桿件的分布載重(均佈+線性變化): SW FEA可以當驗證基準。")
