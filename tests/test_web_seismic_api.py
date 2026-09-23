"""
驗證案例: 網頁後端 /nonlinear_seismic 端點 (動力分析 D6b)

 (1) stdlib 後端: /nonlinear_seismic 的結果跟核心直接呼叫 nonlinear_seismic_web_analysis()+
     seismic_to_dict() 逐項相同(規範反應譜與自訂地震歷程都測); 各種錯誤情況都是清楚的錯誤
     訊息。
 (2) FastAPI 後端(有安裝 fastapi/httpx 才測, 否則略過): 同一個 payload 回傳與 stdlib 完全
     相同; 錯誤情況回 400 與相同訊息; /solve 不受新欄位影響。
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from frame2d import Frame2D
from frame2d.seismic import nonlinear_seismic_web_analysis, seismic_to_dict

# 注意: 不能在模組層級 sys.exit()——pytest 收集測試會 import 這個檔案, SystemExit 會讓整個 pytest
# INTERNALERROR。所有實際工作放在 main() 裡, 只有直接執行(或被 test_zz 的 subprocess 執行)才會跑。

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0


def base_payload(**over):
    def member(mid, i, j, L):
        R = 0.1 * E * I / L
        return {"id": mid, "node_i": i, "node_j": j, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                "R_post_yield_i": R, "R_post_yield_j": R}
    p = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": H, "mx": 5000, "my": 5000},
                  {"id": 2, "x": W, "y": H, "mx": 5000, "my": 5000}, {"id": 3, "x": W, "y": 0}],
        "sections": [{"name": "H", "E": E, "I": I, "A": A}],
        "members": [member(0, 0, 1, H), member(1, 1, 2, W), member(2, 3, 2, H)],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "seismic_control_node": 1, "seismic_direction": "x", "seismic_zeta": 0.05,
        "seismic_damping_modes": [0, 2], "seismic_mass_kind": "lumped",
        "seismic_ground_motion_type": "pulse", "seismic_pulse_amplitude_g": 0.5,
        "seismic_pulse_freq_hz": 1.5, "seismic_pulse_decay": 0.25,
    }
    p.update(over)
    return p


def core_expected(payload):
    f = Frame2D()
    for n in payload["nodes"]:
        f.add_node(n["id"], n["x"], n["y"])
        if n.get("mx") or n.get("my") or n.get("Iz"):
            f.add_mass(n["id"], mx=n.get("mx") or 0.0, my=n.get("my") or 0.0, Iz=n.get("Iz") or 0.0)
    for s in payload["sections"]:
        f.add_section(s["name"], E=s["E"], I=s["I"], A=s["A"])
    for m in payload["members"]:
        f.add_member(m["id"], node_i=m["node_i"], node_j=m["node_j"], section=m["section"],
                     Mp_i=m.get("Mp_i"), Mp_j=m.get("Mp_j"),
                     R_post_yield_i=m.get("R_post_yield_i"), R_post_yield_j=m.get("R_post_yield_j"))
    for sp in payload["supports"]:
        f.support(sp["node"], ux=sp["ux"], uy=sp["uy"], rot=sp["rot"])
    pkg = nonlinear_seismic_web_analysis(
        f, payload["seismic_control_node"], direction=payload["seismic_direction"],
        dt=payload.get("seismic_dt"), n_steps=payload.get("seismic_n_steps"), zeta=payload["seismic_zeta"],
        damping_modes=tuple(payload["seismic_damping_modes"]), mass_kind=payload["seismic_mass_kind"],
        apply_gravity_loads=payload.get("seismic_apply_gravity_loads", True),
        ground_motion_type=payload["seismic_ground_motion_type"],
        pulse_amplitude_g=payload.get("seismic_pulse_amplitude_g"), pulse_freq_hz=payload.get("seismic_pulse_freq_hz"),
        pulse_decay=payload.get("seismic_pulse_decay", 0.0), custom_points_g=payload.get("seismic_custom_points_g"))
    return seismic_to_dict(pkg)


def expect_error(label, fn, contains):
    try:
        fn()
    except Exception as e:
        assert contains in str(e), f"{label}: 錯誤訊息應包含「{contains}」, 實際: {e}"
        print(f"  {label}: 錯誤訊息 OK ({str(e)[:44]}...)")
        return
    raise AssertionError(f"{label}: 應該要報錯")


def main():
    print("=== (1) stdlib 後端 ===")
    try:
        from webapi_stdlib.server import _nonlinear_seismic_payload
    except Exception as e:
        print(f"SKIPPED: 無法匯入 webapi_stdlib.server ({type(e).__name__}: {e})")
        return

    payload = base_payload()
    out = _nonlinear_seismic_payload(payload)
    want = core_expected(payload)
    assert json.loads(json.dumps(out, allow_nan=False)) == json.loads(json.dumps(want)), \
        "/nonlinear_seismic 結果應與核心直接算的完全相同"
    print(f"  規範脈衝地震歷程: /nonlinear_seismic 與核心逐項相同, n_steps={out['n_steps']}, "
          f"peak_disp={out['peak_displacement'] * 1000:.2f} mm, 塑鉸 {[h['label'] for h in out['hinges']]}")

    custom_payload = base_payload(seismic_ground_motion_type="custom", seismic_pulse_amplitude_g=None,
                                  seismic_pulse_freq_hz=None,
                                  seismic_custom_points_g=[[0, 0], [0.3, 0.5], [1.0, -0.3], [3.0, 0.0]],
                                  seismic_dt=0.005, seismic_n_steps=800)
    out_c = _nonlinear_seismic_payload(custom_payload)
    want_c = core_expected(custom_payload)
    assert json.loads(json.dumps(out_c, allow_nan=False)) == json.loads(json.dumps(want_c))
    print(f"  自訂地震歷程: /nonlinear_seismic 與核心逐項相同, n_steps={out_c['n_steps']}")

    expect_error("缺控制節點", lambda: _nonlinear_seismic_payload(base_payload(seismic_control_node=None)),
                "seismic_control_node")
    expect_error("pulse缺振幅", lambda: _nonlinear_seismic_payload(base_payload(seismic_pulse_amplitude_g=None)),
                "pulse_amplitude_g")
    expect_error("custom缺點", lambda: _nonlinear_seismic_payload(base_payload(
        seismic_ground_motion_type="custom", seismic_pulse_amplitude_g=None, seismic_pulse_freq_hz=None)),
                "custom_points_g")
    no_hinge = base_payload(members=[{k: v for k, v in m.items() if k in ("id", "node_i", "node_j", "section")}
                                     for m in base_payload()["members"]])
    expect_error("沒有塑鉸容量", lambda: _nonlinear_seismic_payload(no_hinge), "塑鉸容量")
    expect_error("超過步數上限", lambda: _nonlinear_seismic_payload(
        base_payload(seismic_dt=1e-5, seismic_n_steps=5000)), "上限")
    expect_error("找不到控制節點", lambda: _nonlinear_seismic_payload(base_payload(seismic_control_node=99)),
                "找不到控制節點")

    print("=== (2) FastAPI 後端 ===")
    try:
        from fastapi.testclient import TestClient
        from webapi.main import app
    except Exception as e:
        print(f"  SKIPPED(這一段): 沒有 fastapi/httpx ({type(e).__name__}), 略過 FastAPI 後端比對")
        print("\n全部通過(stdlib 後端; FastAPI 段略過)")
        return

    client = TestClient(app)
    r = client.post("/nonlinear_seismic", json=payload)
    assert r.status_code == 200, r.text
    assert r.json() == json.loads(json.dumps(out)), "FastAPI 與 stdlib 後端的 /nonlinear_seismic 回傳必須完全相同"
    print("  FastAPI 與 stdlib 後端 /nonlinear_seismic 回傳完全相同")
    for label, bad, contains in (
        ("缺控制節點", base_payload(seismic_control_node=None), "seismic_control_node"),
        ("沒有塑鉸容量", no_hinge, "塑鉸容量"),
        ("超過步數上限", base_payload(seismic_dt=1e-5, seismic_n_steps=5000), "上限"),
    ):
        rr = client.post("/nonlinear_seismic", json=bad)
        assert rr.status_code == 400 and contains in rr.json()["detail"], f"{label}: 應回400含「{contains}」: {rr.text[:80]}"
        print(f"  {label}: 400 OK")
    solve_payload = {k: v for k, v in payload.items() if not k.startswith("seismic")}
    solve_payload["point_loads"] = [{"node": 1, "fx": 10000, "fy": 0, "m": 0}]
    rs = client.post("/solve", json=solve_payload)
    assert rs.status_code == 200 and abs(rs.json()["nodes"][1]["ux"]) > 0, "/solve 不應受到新欄位影響"
    print("  /solve 照常求解 OK")
    print("\n全部通過: stdlib 與 FastAPI 兩個後端的 /nonlinear_seismic 行為一致、錯誤訊息清楚、舊呼叫端不受影響")


if __name__ == "__main__":
    main()
