"""
其餘測試檔案(test_cantilever.py, test_truss.py, ...)是寫成「直接執行的
驗證腳本」風格(print+assert, 執行時會印出詳細的比對數值方便人眼檢查),
不是pytest的def test_*()函式風格。這支檔案讓pytest也能發現、執行它們
(用subprocess跑, assert exit code==0), 這樣單純執行`pytest`(不加任何
參數)就能涵蓋全部驗證, 不會漏掉。

這是2026-08-29跟使用者對話時, 對方引用的chatGPT分析指出「zip裡有
output_plots/tests污染」, 查證後發現那個具體問題不存在, 但過程中額外
抓到這個真的問題(pytest --collect-only只收集到11個測試, 其餘16個檔案
完全沒被發現)才補上的。
"""
import subprocess
import sys
from pathlib import Path

import pytest

_TESTS_DIR = Path(__file__).parent
_SCRIPT_STYLE_TESTS = sorted(
    p.name for p in _TESTS_DIR.glob("test_*.py")
    if p.name != Path(__file__).name
    and "def test_" not in p.read_text(encoding="utf-8")
)


@pytest.mark.parametrize("script_name", _SCRIPT_STYLE_TESTS)
def test_script_style_validation(script_name):
    """跑一支腳本風格的驗證檔案, 用exit code判斷成功/失敗(腳本內部本身
    已經用assert在檢查數值, exit code非0代表某個assert失敗或整體丟了例外)。"""
    script_path = _TESTS_DIR / script_name
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=_TESTS_DIR.parent,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if result.returncode != 0:
        pytest.fail(f"{script_name} 失敗 (exit code {result.returncode}):\n"
                    f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}")
