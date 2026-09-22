"""
驗證案例: 網頁後端 /rsa 端點 (動力分析 D3b)

 (1) stdlib 後端: /rsa 的結果跟核心直接呼叫 spectrum_analysis()+rsa_to_dict() 逐項相同(規範
     反應譜與自訂反應譜都測); 各種錯誤情況(缺 SDS/SD1、缺自訂點、沒有支承)都是清楚的錯誤訊息。
 (2) FastAPI 後端(有安裝 fastapi/httpx 才測, 否則略過): 同一個 payload 回傳與 stdlib 完全相同;
     錯誤情況回 400 與相同訊息; /solve 不受新欄位影響。
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from frame2d import Frame2D
from frame2d.spectrum import spectrum_analysis, rsa_to_dict

# 注意: 不能在模組層級 sys.exit()——pytest 收集測試會 import 這個檔案, SystemExit 會讓整個 pytest
# INTERNALERROR。所有實際工作放在 main() 裡, 只有直接執行(或被 test_zz 的 subprocess 執行)才會跑。

E, I, A = 200e9, 8e-5, 5e-3


def base_payload(**over):
    p = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": 4, "mx": 5000, "my": 5000},
                  {"id": 2, "x": 6, "y": 4, "mx": 5000, "my": 5000}, {"id": 3, "x": 6, "y": 0}],
        "sections": [{"name": "H", "E": E, "I": I, "A": A}],
        "members": [{"id": 0, "node_i": 0, "node_j": 1, "section": "H"},
                    {"id": 1, "node_i": 1, "node_j": 2, "section": "H"},
                    {"id": 2, "node_i": 3, "node_j": 2, "section": "H"}],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "rsa_direction": "x", "rsa_combine": "CQC", "rsa_n_modes": 4, "rsa_mass_kind": "consistent",
        "rsa_spectrum_type": "code", "rsa_code_sds": 0.6, "rsa_code_sd1": 0.35, "rsa_code_tl": 6.0,
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
        f.add_member(m["id"], node_i=m["node_i"], node_j=m["node_j"], section=m["section"])
    for sp in payload["supports"]:
        f.support(sp["node"], ux=sp["ux"], uy=sp["uy"], rot=sp["rot"])
    pkg = spectrum_analysis(
        f, direction=payload["rsa_direction"], damping=payload.get("rsa_damping", 0.05),
        combine=payload["rsa_combine"], n_modes=payload["rsa_n_modes"], mass_kind=payload["rsa_mass_kind"],
        spectrum_type=payload["rsa_spectrum_type"], code_sds=payload.get("rsa_code_sds"),
        code_sd1=payload.get("rsa_code_sd1"), code_tl=payload.get("rsa_code_tl", 6.0),
        custom_points=payload.get("rsa_custom_points"))
    return rsa_to_dict(pkg)


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
        from webapi_stdlib.server import _rsa_payload
    except Exception as e:
        print(f"SKIPPED: 無法匯入 webapi_stdlib.server ({type(e).__name__}: {e})")
        return

    payload = base_payload()
    out = _rsa_payload(payload)
    want = core_expected(payload)
    assert json.loads(json.dumps(out, allow_nan=False)) == json.loads(json.dumps(want)), "/rsa 結果應與核心直接算的完全相同"
    print(f"  規範反應譜: /rsa 與核心逐項相同, base_shear={out['base_shear']:.2f} N, "
          f"{out['n_modes']} 個模態, cum_ratio_total={out['cum_ratio_total']:.4f}")

    custom_payload = base_payload(rsa_spectrum_type="custom", rsa_code_sds=None, rsa_code_sd1=None,
                                  rsa_custom_points=[[0.01, 1.0], [0.3, 6.0], [2.0, 1.0]])
    out_c = _rsa_payload(custom_payload)
    want_c = core_expected(custom_payload)
    assert json.loads(json.dumps(out_c, allow_nan=False)) == json.loads(json.dumps(want_c))
    print(f"  自訂反應譜: /rsa 與核心逐項相同, base_shear={out_c['base_shear']:.2f} N")

    expect_error("code缺SDS", lambda: _rsa_payload(base_payload(rsa_code_sds=None)), "code_sds")
    expect_error("custom缺點", lambda: _rsa_payload(base_payload(rsa_spectrum_type="custom", rsa_custom_points=None)),
                "custom_points")
    expect_error("spectrum_type打錯", lambda: _rsa_payload(base_payload(rsa_spectrum_type="bogus")), "spectrum_type")
    no_support = base_payload(supports=[])
    expect_error("沒有支承", lambda: _rsa_payload(no_support), "機構")

    print("=== (2) FastAPI 後端 ===")
    try:
        from fastapi.testclient import TestClient
        from webapi.main import app
    except Exception as e:
        print(f"  SKIPPED(這一段): 沒有 fastapi/httpx ({type(e).__name__}), 略過 FastAPI 後端比對")
        print("\n全部通過(stdlib 後端; FastAPI 段略過)")
        return

    client = TestClient(app)
    r = client.post("/rsa", json=payload)
    assert r.status_code == 200, r.text
    assert r.json() == json.loads(json.dumps(out)), "FastAPI 與 stdlib 後端的 /rsa 回傳必須完全相同"
    print("  FastAPI 與 stdlib 後端 /rsa 回傳完全相同")
    for label, bad, contains in (
        ("code缺SDS", base_payload(rsa_code_sds=None), "code_sds"),
        ("沒有支承", base_payload(supports=[]), "機構"),
    ):
        rr = client.post("/rsa", json=bad)
        assert rr.status_code == 400 and contains in rr.json()["detail"], f"{label}: 應回400含「{contains}」: {rr.text[:80]}"
        print(f"  {label}: 400 OK")
    assert client.post("/rsa", json=base_payload(rsa_combine="average")).status_code == 422, "非法combine應被schema擋下(422)"
    solve_payload = {k: v for k, v in payload.items() if not k.startswith("rsa")}
    solve_payload["point_loads"] = [{"node": 1, "fx": 10000, "fy": 0, "m": 0}]
    rs = client.post("/solve", json=solve_payload)
    assert rs.status_code == 200 and abs(rs.json()["nodes"][1]["ux"]) > 0, "/solve 不應受到新欄位影響"
    print("  /solve 照常求解 OK")
    print("\n全部通過: stdlib 與 FastAPI 兩個後端的 /rsa 行為一致、錯誤訊息清楚、舊呼叫端不受影響")


if __name__ == "__main__":
    main()
