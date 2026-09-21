"""
網頁端到端測試(選用): 反覆載重(遲滯)分析 (D7b)

啟動真正的 stdlib 後端(webapi_stdlib.server.Handler, 隨機埠), 用 jsdom 載入它提供的 index.html,
走「匯入JSON」載入範例模型 examples/portal_cyclic_demo.json, 選「循環」分析、填幅值/圈數/步長、Solve,
檢查遲滯迴圈圖、塑鉸 M-θp 圖、結果表、單位切換、錯誤訊息與重置, 並把網頁拿到的耗能/最大底剪力/塑鉸
降伏順序跟 frame2d 核心直接算的結果逐一比對。

需要 node 與 jsdom(npm install jsdom); 沒有就印 SKIPPED 並正常結束(pytest 會回報成 skip)。
實際的操作與斷言在 tests/web/cyclic_e2e.js。
"""
import json
import shutil
import subprocess
import sys
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "examples" / "portal_cyclic_demo.json"


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
    from webapi_stdlib.server import Handler, _cyclic_payload

    # 核心直接算的預期值(同一個 JSON 模型 + 同樣的載重歷程參數, SI 單位)
    payload = json.loads(MODEL.read_text(encoding="utf-8"))
    payload.update({"cyclic_control_nodes": [1], "cyclic_direction": "x", "cyclic_amplitudes": [0.015, 0.03, 0.045, 0.06],
                    "cyclic_n_cycles": 2, "cyclic_step": 0.002})
    out = _cyclic_payload(payload)
    expected = {
        "loops": [l["energy"] for l in out["loops"]],
        "hingeLabels": [h["label"] for h in out["hinges"]],
        "maxF": max(abs(v) for v in out["F"]),
        "nPoints": len(out["u"]),
        "nYield": sum(1 for e in out["events"] if e["kind"] == "yield"),
    }

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        res = subprocess.run([node, str(ROOT / "tests" / "web" / "cyclic_e2e.js"), base, str(MODEL), json.dumps(expected)],
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
