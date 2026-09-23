"""
網頁端到端測試(選用): 非線性地震反應分析 (D6b)

啟動真正的 stdlib 後端(webapi_stdlib.server.Handler, 隨機埠), 用 jsdom 載入它提供的
index.html, 走「匯入JSON」載入範例模型 examples/portal_seismic_demo.json, 選「非線性地震」
分析、衰減脈衝地震歷程、Solve, 檢查地震動畫(播放/暫停/滑桿)、時程圖、遲滯迴圈圖、塑鉸M-θp圖、
結果表、單位切換、自訂地震歷程、錯誤訊息與重置, 並把網頁拿到的步數/尖峰位移/動畫影格數/能量
平衡跟 frame2d 核心 nonlinear_seismic_web_analysis() 直接算的結果逐一比對。

需要 node 與 jsdom(npm install jsdom); 沒有就印 SKIPPED 並正常結束(pytest 會回報成 skip)。
實際的操作與斷言在 tests/web/seismic_e2e.js。
"""
import json
import shutil
import subprocess
import sys
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "examples" / "portal_seismic_demo.json"


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
    from webapi_stdlib.server import Handler, _nonlinear_seismic_payload

    payload = json.loads(MODEL.read_text(encoding="utf-8"))
    payload.update({"seismic_control_node": 1, "seismic_direction": "x", "seismic_zeta": 0.05,
                    "seismic_damping_modes": [0, 2], "seismic_mass_kind": "lumped",
                    "seismic_ground_motion_type": "pulse", "seismic_pulse_amplitude_g": 0.5,
                    "seismic_pulse_freq_hz": 1.5, "seismic_pulse_decay": 0.25})
    out = _nonlinear_seismic_payload(payload)
    energy_final = out["energy"]["KE"][-1] + out["energy"]["Wdamp"][-1] + out["energy"]["Wint"][-1]
    expected = {"nSteps": out["n_steps"], "peakDisp": out["peak_displacement"], "nFrames": len(out["frames"]),
               "hingeLabels": [h["label"] for h in out["hinges"]], "energyFinal": energy_final}

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        res = subprocess.run([node, str(ROOT / "tests" / "web" / "seismic_e2e.js"), base, str(MODEL), json.dumps(expected)],
                             capture_output=True, text=True, cwd=str(ROOT), timeout=150)
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
