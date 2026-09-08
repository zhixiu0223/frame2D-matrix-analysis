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
from frame2d.dofmanager import solve_pdelta, initial_hinge_states
from frame2d.pushover import run_pushover, apply_gravity
from frame2d.postprocess import member_internal_forces

from .diagrams import build_diagrams_and_deformed, build_deformed_with_scale
from .storage import LocalFileStorage, InvalidNameError, NotFoundError
from .pdf_export import build_pdf_report, build_fbd_previews, build_fbd_images_archive, build_pushover_pdf_report
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
                     release_i=m.get("release_i", False), release_j=m.get("release_j", False),
                     Mp_i=m.get("Mp_i"), Mp_j=m.get("Mp_j"),
                     R_post_yield_i=m.get("R_post_yield_i"), R_post_yield_j=m.get("R_post_yield_j"))
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


def _build_and_solve(payload: dict, analysis_type: str = None):
    """_build_frame()+solve() 的共用包裝, 統一處理「模型內有殘留參照」
    這類錯誤(例如分割/刪除桿件後, 還留著指向舊桿件編號的均佈載重),
    這種情況下 frame2d 底層會丟出 KeyError(9) 這種只印一個數字的
    錯誤, 前端顯示出來完全看不懂在講什麼, 這裡轉成 ValueError 帶
    看得懂的訊息(KeyError 的 str() 只會印裸的 key 值, ValueError
    才會把完整句子印出來)。solve_pdelta()疊代不收斂丟出的RuntimeError
    原樣往上冒出去, do_POST()那層統一當成400處理。

    analysis_type: None時用payload.get('analysis_type', 'linear');
    呼叫端可以明確覆寫(見_solve_payload裡pdelta額外算一次線性基準時
    的用法)。"""
    if analysis_type is None:
        analysis_type = payload.get("analysis_type", "linear")
    f = _build_frame(payload)
    try:
        if analysis_type == "pdelta":
            result = solve_pdelta(f)
        else:
            result = solve(f)
    except KeyError as e:
        raise ValueError(
            f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
            f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
            f"還留著指向舊桿件編號), 請檢查並移除"
        )
    return f, result


def _build_hinge_summary(hs_final):
    """跟webapi/main.py同名函式邏輯一致(把{member_id: HingeState}轉成
    JSON安全的摘要dict), /solve(pushover)跟/export/pdf(pushover)共用。"""
    return {
        str(mid): {
            "Mp": [None if math.isinf(hs.Mp[0]) else float(hs.Mp[0]),
                   None if math.isinf(hs.Mp[1]) else float(hs.Mp[1])],
            "yielded": [bool(hs.yielded[0]), bool(hs.yielded[1])],
            "theta_p": [float(hs.theta_p[0]), float(hs.theta_p[1])],
            "performance_level": [hs.performance_level(0), hs.performance_level(1)],
        }
        for mid, hs in hs_final.items()
    }


def _build_event_log_out(event_log):
    """跟webapi/main.py同名函式邏輯一致。"""
    return [
        {"u": float(ev["u"]), "F": float(ev["F"]),
         "yielded": [[mid, end_idx] for mid, end_idx in ev["yielded"]]}
        for ev in event_log
    ]


def _prepare_pushover_run(payload: dict):
    """跟webapi/main.py的_prepare_pushover_run()邏輯一致, 輸入是dict
    (這個檔案不用Pydantic的FrameIn)。回傳(f, run_kwargs), run_kwargs
    可以直接**展開餵給run_pushover()。"""
    control_nodes = payload.get("pushover_control_nodes")
    if control_nodes is None:
        control_node = payload.get("pushover_control_node")
        if control_node is None:
            raise ValueError("pushover需要指定pushover_control_node(單點)或pushover_control_nodes(多點)其中一個")
        control_nodes = [control_node]
    weights = payload.get("pushover_weights")
    if weights is None:
        weights = [1.0] * len(control_nodes)
    elif len(weights) != len(control_nodes):
        raise ValueError(f"pushover_weights長度({len(weights)})必須跟pushover_control_nodes長度({len(control_nodes)})一樣")

    target = payload.get("pushover_target")
    step = payload.get("pushover_step")
    if target is None or step is None:
        raise ValueError("pushover需要指定pushover_target、pushover_step兩個欄位")

    f = _build_frame(payload)
    hinge_states = initial_hinge_states(f)
    if not hinge_states:
        raise ValueError(
            "模型裡沒有任何桿件設定Mp_i(塑鉸容量), 無法做pushover -- "
            "至少要有一根frame桿件設定Mp_i/Mp_j/R_post_yield_i/R_post_yield_j"
        )

    direction_key = payload.get("pushover_direction", "x")
    local_idx = {"x": 0, "y": 1}[direction_key]
    control_dofs = [f.dofs_of(n)[local_idx] for n in control_nodes]
    base_reaction_dofs = list(dict.fromkeys(
        f.dofs_of(s.node)[local_idx] for s in f.supports))

    initial_cum_forces = None
    has_gravity_loads = bool(f.point_loads or f.distributed_loads or f.member_point_loads)
    if has_gravity_loads:
        initial_cum_forces, _ = apply_gravity(f, hinge_states)

    return f, dict(
        frame=f, hinge_states=hinge_states, prescribed_dofs=control_dofs, direction=weights,
        target_total=target, d_nominal=step,
        base_reaction_dofs=base_reaction_dofs, initial_cum_forces=initial_cum_forces,
        use_pdelta=payload.get("pushover_use_pdelta", False),
        mechanism_ratio_limit=payload.get("pushover_mechanism_ratio_limit", 1e-8),
        control_mode=payload.get("pushover_control_mode", "displacement"),
        geometry_update=payload.get("pushover_geometry_update", False),
    )


def _solve_pushover_payload(payload: dict) -> dict:
    """analysis_type='pushover'的獨立處理路徑, 邏輯跟webapi/main.py的
    _solve_pushover()完全一致, 只是輸入是普通dict。"""
    f, run_kwargs = _prepare_pushover_run(payload)
    history_u, history_F, event_log, hs_final, mechanism, snapshots = run_pushover(
        include_snapshots=True, **run_kwargs)

    hinge_summary = _build_hinge_summary(hs_final)
    event_log_out = _build_event_log_out(event_log)
    snapshots_out = [
        {
            "member_forces": {str(mid): f for mid, f in snap["member_forces"].items()},
            "hinge_states": {str(mid): hs for mid, hs in snap["hinge_states"].items()},
        }
        for snap in snapshots
    ]

    return {
        "analysis_type": "pushover",
        "history_u": [float(v) for v in history_u],
        "history_F": [float(v) for v in history_F],
        "event_log": event_log_out,
        "mechanism_reached": bool(mechanism),
        "hinge_summary": hinge_summary,
        "history_snapshots": snapshots_out,
    }


def _export_pushover_pdf(payload: dict) -> bytes:
    """跟webapi/main.py的_export_pushover_pdf()邏輯一致, 輸入是dict。"""
    f, run_kwargs = _prepare_pushover_run(payload)
    (history_u, history_F, event_log, hs_final, mechanism,
     u_full_cum, snapshots, cum_reaction, max_rotation) = run_pushover(
        include_final_displacement=True, include_snapshots=True,
        include_final_reactions=True, include_max_rotation=True, **run_kwargs)

    pushover_result = {
        "history_u": history_u, "history_F": history_F,
        "event_log": _build_event_log_out(event_log),
        "hinge_summary": _build_hinge_summary(hs_final),
        "mechanism_reached": bool(mechanism),
        "cum_forces": snapshots[-1]["member_forces"],
        "u_full_cum": u_full_cum,
        "cum_reaction": cum_reaction,
    }
    return build_pushover_pdf_report(f, pushover_result, units=payload.get("units"))


def _solve_payload(payload: dict) -> dict:
    if payload.get("analysis_type", "linear") == "pushover":
        return _solve_pushover_payload(payload)

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
    analysis_type = payload.get("analysis_type", "linear")

    out = {
        "nodes": node_out,
        "members": member_out,
        "diagrams": diagrams,
        "deformed": deformed,
        "deform_scale": deform_scale,
        "analysis_type": analysis_type,
    }

    if analysis_type == "pdelta":
        # 額外算一次純線性基準(P=0.0)做對照, 理由跟webapi/main.py同一段
        # 註解一致: 線性理論上比pdelta更不容易發散, 這裡不特別catch,
        # 讓例外照do_POST()既有的統一400處理走。
        _, result_linear = _build_and_solve(payload, analysis_type="linear")
        max_disp_linear = result_linear.max_displacement()
        max_disp_pdelta = result.max_displacement()
        v_linear = max_disp_linear.value if max_disp_linear is not None else 0.0
        v_pdelta = max_disp_pdelta.value if max_disp_pdelta is not None else 0.0
        out["pdelta_comparison"] = {
            "max_disp_linear": float(v_linear),
            "max_disp_pdelta": float(v_pdelta),
            "amplification": float(v_pdelta / v_linear) if v_linear > 1e-12 else None,
            "iterations": result.pdelta_iterations,
        }
        out["deformed_linear"] = build_deformed_with_scale(f, result_linear, deform_scale)

    return out


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
                if payload.get("analysis_type") == "pushover":
                    try:
                        pdf_bytes = _export_pushover_pdf(payload)
                    except KeyError as e:
                        raise ValueError(
                            f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照, 請檢查並移除"
                        )
                    self._send_bytes(pdf_bytes, "application/pdf",
                                      extra_headers={"Content-Disposition": 'attachment; filename="frame2d_pushover_report.pdf"'})
                    return
                f = _build_frame(payload)
                try:
                    pdf_bytes = build_pdf_report(f, units=payload.get("units"), member_ids=payload.get("member_ids"),
                                                  include_member_diagrams=not payload.get("fbd_only", False))
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
        if self.path == "/export/fbd_images":
            try:
                payload = self._read_json_body()
                f = _build_frame(payload)
                member_ids = payload.get("member_ids")
                if not member_ids:
                    raise ValueError("沒有指定要匯出的桿件(member_ids)")
                try:
                    zip_bytes = build_fbd_images_archive(f, member_ids, units=payload.get("units"))
                except KeyError as e:
                    raise ValueError(
                        f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                        f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                        f"還留著指向舊桿件編號), 請檢查並移除"
                    )
                if len(member_ids) == 1:
                    import io as _io
                    import zipfile as _zipfile
                    with _zipfile.ZipFile(_io.BytesIO(zip_bytes)) as zf:
                        names = zf.namelist()
                        if names:
                            png_bytes = zf.read(names[0])
                            self._send_bytes(png_bytes, "image/png",
                                              extra_headers={"Content-Disposition": f'attachment; filename="member_{member_ids[0]}_fbd.png"'})
                            return
                self._send_bytes(zip_bytes, "application/zip",
                                  extra_headers={"Content-Disposition": 'attachment; filename="frame2d_fbd_images.zip"'})
            except Exception as e:
                self._send_json({"error": str(e)}, status=400)
            return
        if self.path == "/preview/fbd":
            try:
                payload = self._read_json_body()
                f = _build_frame(payload)
                member_ids = payload.get("member_ids")
                if not member_ids:
                    raise ValueError("沒有指定要預覽的桿件(member_ids)")
                try:
                    previews = build_fbd_previews(f, member_ids, units=payload.get("units"))
                except KeyError as e:
                    raise ValueError(
                        f"找不到 ID 為 {e} 的節點或桿件, 模型內有殘留的參照"
                        f"(常見情況: 分割或刪除桿件後, 均佈載重/桿件集中力"
                        f"還留著指向舊桿件編號), 請檢查並移除"
                    )
                self._send_json({"previews": previews})
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
