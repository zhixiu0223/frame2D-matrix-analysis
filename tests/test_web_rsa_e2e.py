"""
網頁端到端測試(選用): 反應譜分析 (D3b)

啟動真正的 stdlib 後端(webapi_stdlib.server.Handler, 隨機埠), 用 jsdom 載入它提供的
index.html, 走「匯入JSON」載入範例模型 examples/portal_rsa_demo.json, 選「反應譜」分析、
規範反應譜(SDS/SD1/TL)、Solve, 檢查反應譜曲線圖、結果表(模態/節點位移/桿件內力)、單位切換、
自訂反應譜、錯誤訊息與重置, 並把網頁拿到的基底剪力/週期/累積質量比跟 frame2d 核心
spectrum_analysis() 直接算的結果逐一比對。

需要 node 與 jsdom(npm install jsdom); 沒有就印 SKIPPED 並正常結束(pytest 會回報成 skip)。
實際的操作與斷言在 tests/web/rsa_e2e.js。
"""
import json
import shutil
import subprocess
import sys
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "examples" / "portal_rsa_demo.json"


def main():
    node = shutil.which("node")
    if not node:
        print("SKIPPED: 沒有安裝 node, 略過網頁端到端測試")
        return 0
    probe = subprocess.run([node, "-e", "require('jsdom')"], capture_output=True, text=True, cwd=str(ROOT))
    if probe.returncode != 0:
        print("SKIPPED: 沒有安裝 jsdom (npm install jsdom), 略過網頁端到端測試")
        return 0

    sys.path.insert(0, str(ROOT))
    from webapi_stdlib.server import Handler, _rsa_payload

    payload = json.loads(MODEL.read_text(encoding="utf-8"))
    payload.update({"rsa_direction": "x", "rsa_combine": "CQC", "rsa_n_modes": 4, "rsa_mass_kind": "consistent",
                    "rsa_spectrum_type": "code", "rsa_code_sds": 0.6, "rsa_code_sd1": 0.35, "rsa_code_tl": 6.0})
    out = _rsa_payload(payload)
    expected = {
        "baseShear": out["base_shear"], "cumRatioTotal": out["cum_ratio_total"],
        "periods": [m["period"] for m in out["modes"]], "nModes": out["n_modes"],
        "gammas": [m["gamma"] for m in out["modes"]],
    }

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        res = subprocess.run([node, str(ROOT / "tests" / "web" / "rsa_e2e.js"), base, str(MODEL), json.dumps(expected)],
                             capture_output=True, text=True, cwd=str(ROOT), timeout=120)
    finally:
        server.shutdown()
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        print("失敗: 網頁端到端測試沒有通過")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
