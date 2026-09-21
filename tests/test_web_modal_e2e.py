"""
網頁端到端測試(選用): 模態分析 (D2b)

啟動真正的 stdlib 後端(webapi_stdlib.server.Handler, 隨機埠), 用 jsdom 載入它提供的
index.html, 像使用者一樣操作介面(節點質量輸入 -> 選「模態」-> Solve -> 檢查結果表、振型、
單位切換、錯誤訊息、重置), 並把網頁拿到的週期跟 frame2d 核心直接算的結果逐一比對。

需要 node 與 jsdom(npm install jsdom); 沒有就印 SKIPPED 並正常結束(pytest 會回報成 skip)。
實際的操作與斷言在 tests/web/modal_e2e.js。
"""
import json
import os
import shutil
import subprocess
import sys
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


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
    from frame2d import Frame2D
    from frame2d.modal import eigen
    from webapi_stdlib.server import Handler

    # 核心直接算的預期週期(跟網頁模型相同: SI、ρ=7850、兩個樓層節點各 5 t = 5000 kg)
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('s', E=200e9, I=8e-5, A=1e-2, rho=7850.0)
    f.add_member(0, 0, 1, 's').add_member(1, 1, 2, 's').add_member(2, 3, 2, 's')
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    md = eigen(f, n_modes=3, mass='consistent')
    expected = {"periods": [float(t) for t in md.period], "nodeMassKg": 5000.0}

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        res = subprocess.run([node, str(ROOT / "tests" / "web" / "modal_e2e.js"), base, json.dumps(expected)],
                             capture_output=True, text=True, cwd=str(ROOT), timeout=90)
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
