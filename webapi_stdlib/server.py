"""
frame2d 的最小 Web API 層 -- 純 Python 標準函式庫版本(不依賴 FastAPI/pydantic)。

會有這個版本是因為 Termux(Android/arm64 + Python 3.14 這種冷門組合)
在 PyPI 上抓不到 pydantic-core 的預編譯 wheel,pip 會 fallback 成在手機上
編譯 Rust,非常慢甚至可能失敗。這裡改用 http.server + json,零額外依賴,
邏輯跟 webapi/(FastAPI 版)完全一致,只是拿掉 pydantic 的自動驗證改成
手動檢查。FastAPI 版留著給 Colab / GCP 用(那邊是標準 x86_64 Linux,
pip 抓得到現成 wheel,不會有這個問題)。

執行方式: python3 -m webapi_stdlib.server
"""
import json
import math
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

from .diagrams import build_diagrams_and_deformed
from .storage import LocalFileStorage, InvalidNameError, NotFoundError
from .pdf_export import build_pdf_report
from .query_point import query_point

STATIC_DIR = Path(__file__).parent / "static"
# saved_models/ 放在 repo 根目錄, 跟 webapi/(FastAPI版)共用同一份存檔,
# 不管在 Termux 上啟動哪一個後端, 存檔清單都是同一份。
storage = LocalFileStorage(Path(__file__).parent.parent / "saved_models")


def _build_frame(payload: dict) -> Frame2D:
    """把 JSON payload 轉成 frame2d.Frame2D。跟 webapi/main.py 的
    _build_frame() 邏輯完全一致,只是輸入是普通 dict 不是 pydantic model。"""
    f = Frame2D()
    for n in payload.get("nodes", []):
        f.add_node(n["id"], n["x"], n["y"])
    for s in payload.get("sections", []):
        f.add_section(s["name"], E=s["E"], I=s["I"], A=s.get("A", 1e8))
    for m in payload.get("members", []):
        f.add_member(m["id"], node_i=m["node_i"], node_j=m["node_j"], section=m["section"],
                     member_type=m.get("member_type", "frame"),
                     release_i=m.get("release_i", False), release_j=m.get("release_j", False))
    for sp in payload.get("supports", []):
        f.support(sp["node"], ux=sp.get("ux"), uy=sp.get("uy"), rot=sp.get("rot"))
    for pl in payload.get("point_loads", []):
        f.point_load(pl["node"], fx=pl.get("fx", 0.0), fy=pl.get("fy", 0.0), m=pl.get("m", 0.0))
    for dl in payload.get("distributed_loads", []):
        f.distributed_load(dl["member"], w=dl["w_start"], w_end=dl.get("w_end"),
                            x_start=dl.get("x_start"), x_end=dl.get("x_end"),
                            direction=dl.get("direction", "local"), angle_deg=dl.get("angle_deg"))
    for mpl in payload.get("member_point_loads", []):
        f.member_point_load(mpl["member"], a=mpl["a"], fx=mpl.get("fx", 0.0),
                             fy=mpl.get("fy", 0.0), m=mpl.get("m", 0.0),
                             direction=mpl.get("direction", "local"),
                             F=mpl.get("F"), angle_deg=mpl.get("angle_deg"))
    return f


def _build_and_solve(payload: dict):
    """_build_frame()+solve() 的共用包裝, 統一處理「模型內有殘留參照」
    這類錯誤(例如分割/刪除桿件後, 還留著指向舊桿件編號的均佈載重),
    這種情況下 frame2d 底層會丟出 KeyError(9) 這種只印一個數字的
    錯誤, 前端顯示出來完全看不懂在講什麼, 這裡轉成 ValueError 帶
    看得懂的訊息(KeyError 的 str() 只會印裸的 key 值, ValueError
    才會把完整句子印出來)。"""
    f = _build_frame(payload)
    try:
        result = solve(f)
    except KeyError as e:
        raise ValueError(
            f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
            f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
            f"還留著指向舊桿件編號), 請檢查並移除"
        )
    return f, result


def _solve_payload(payload: dict) -> dict:
    f, result = _build_and_solve(payload)

    node_out = []
    for n in payload.get("nodes", []):
        nid = n["id"]
        ux_i, uy_i, rot_i = f.dofs_of(nid)
        node_out.append({
            "node": nid,
            "ux": float(result.displacements[ux_i]),
            "uy": float(result.displacements[uy_i]),
            "rot": float(result.displacements[rot_i]),
            "Rx": float(result.reactions[ux_i]),
            "Ry": float(result.reactions[uy_i]),
            "M": float(result.reactions[rot_i]),
        })

    member_out = []
    for mid, mr in result.member_results.items():
        x, N, V, M = member_internal_forces(f, result, mid, n=2)
        member_out.append({
            "member_id": mid,
            "L": float(mr.L),
            "angle_deg": float(math.degrees(mr.angle)),
            "N1": float(N[0]), "V1": float(V[0]), "M1": float(M[0]),
            "N2": float(N[-1]), "V2": float(V[-1]), "M2": float(M[-1]),
            "slack": bool(mr.slack),
        })

    diagrams, deformed, deform_scale = build_diagrams_and_deformed(f, result)

    return {
        "nodes": node_out,
        "members": member_out,
        "diagrams": diagrams,
        "deformed": deformed,
        "deform_scale": deform_scale,
    }


class Handler(BaseHTTPRequestHandler):
    def _send_json(self, obj, status=200):
        body = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_bytes(self, data: bytes, content_type: str, status=200, extra_headers=None):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        for k, v in (extra_headers or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(data)

    def _read_json_body(self):
        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length)
        return json.loads(raw)

    def do_GET(self):
        if self.path == "/" or self.path == "/index.html":
            html = (STATIC_DIR / "index.html").read_bytes()
            self._send_bytes(html, "text/html; charset=utf-8")
            return
        if self.path == "/models":
            self._send_json({"names": storage.list()})
            return
        if self.path.startswith("/models/"):
            name = urllib.parse.unquote(self.path[len("/models/"):])
            try:
                self._send_json(storage.load(name))
            except InvalidNameError as e:
                self._send_json({"error": str(e)}, status=400)
            except NotFoundError as e:
                self._send_json({"error": str(e)}, status=404)
            return
        self.send_error(404)

    def do_POST(self):
        if self.path == "/solve":
            try:
                payload = self._read_json_body()
                self._send_json(_solve_payload(payload))
            except Exception as e:
                self._send_json({"error": str(e)}, status=400)
            return
        if self.path == "/export/pdf":
            try:
                payload = self._read_json_body()
                f = _build_frame(payload)
                try:
                    pdf_bytes = build_pdf_report(f)
                except KeyError as e:
                    raise ValueError(
                        f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                        f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                        f"還留著指向舊桿件編號), 請檢查並移除"
                    )
                self._send_bytes(pdf_bytes, "application/pdf",
                                  extra_headers={"Content-Disposition": 'attachment; filename="frame2d_report.pdf"'})
            except Exception as e:
                self._send_json({"error": str(e)}, status=400)
            return
        if self.path == "/query_point":
            try:
                payload = self._read_json_body()
                f, result = _build_and_solve(payload)
                out = query_point(f, result, payload["member"], payload["mode"], payload["value"])
                self._send_json(out)
            except Exception as e:
                self._send_json({"error": str(e)}, status=400)
            return
        if self.path.startswith("/models/"):
            name = urllib.parse.unquote(self.path[len("/models/"):])
            try:
                payload = self._read_json_body()
                storage.save(name, payload)
                self._send_json({"saved": name})
            except InvalidNameError as e:
                self._send_json({"error": str(e)}, status=400)
            return
        self.send_error(404)

    def do_DELETE(self):
        if self.path.startswith("/models/"):
            name = urllib.parse.unquote(self.path[len("/models/"):])
            try:
                storage.delete(name)
                self._send_json({"deleted": name})
            except InvalidNameError as e:
                self._send_json({"error": str(e)}, status=400)
            except NotFoundError as e:
                self._send_json({"error": str(e)}, status=404)
            return
        self.send_error(404)

    def log_message(self, fmt, *args):
        # 精簡一點, 只印方法+路徑+狀態碼
        pass


def main(host="0.0.0.0", port=8000):
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"frame2d API (stdlib) 啟動: http://{host}:{port}  (Ctrl+C 結束)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.shutdown()


if __name__ == "__main__":
    main()
