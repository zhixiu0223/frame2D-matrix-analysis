"""
驗證案例: 網頁後端 /modal 端點 (動力分析 D2b)

 (1) stdlib 後端: 節點質量 mx/my/Iz 與斷面 rho 經 JSON 抵達 frame2d, /modal 的結果跟核心直接算
     eigen()+modal_to_dict() 逐項相同; 錯誤情況(無支承、沒有質量、cable、equal_dof、非零支座
     沉陷、非法質量種類)都是清楚的錯誤訊息。
 (2) FastAPI 後端(有安裝 fastapi/httpx 才測, 否則略過這一段): 同一個 payload 回傳與 stdlib
     後端完全相同; 錯誤情況回 400 與相同訊息; /solve 不受新欄位影響。
兩個後端各自有一份 _build_frame, 這支測試守住它們對質量欄位的處理不漂移。
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from frame2d import Frame2D
from frame2d.modal import eigen, modal_to_dict

# 注意: 不能在模組層級 sys.exit()——pytest 收集測試會 import 這個檔案, SystemExit 會讓整個 pytest
# INTERNALERROR。所有實際工作放在 main() 裡, 只有直接執行(或被 test_zz 的 subprocess 執行)才會跑。


def base_payload(**over):
    p = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": 4, "mx": 5000, "my": 5000},
                  {"id": 2, "x": 6, "y": 4, "mx": 5000, "my": 5000, "Iz": 100}, {"id": 3, "x": 6, "y": 0}],
        "sections": [{"name": "s", "E": 200e9, "I": 8e-5, "A": 1e-2, "rho": 7850}],
        "members": [{"id": 0, "node_i": 0, "node_j": 1, "section": "s"},
                    {"id": 1, "node_i": 1, "node_j": 2, "section": "s"},
                    {"id": 2, "node_i": 3, "node_j": 2, "section": "s"}],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "modal_n_modes": 4, "modal_mass_kind": "consistent",
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
        f.add_section(s["name"], E=s["E"], I=s["I"], A=s["A"], rho=s.get("rho"))
    for m in payload["members"]:
        f.add_member(m["id"], node_i=m["node_i"], node_j=m["node_j"], section=m["section"])
    for sp in payload["supports"]:
        f.support(sp["node"], ux=sp["ux"], uy=sp["uy"], rot=sp["rot"])
    return modal_to_dict(eigen(f, n_modes=payload["modal_n_modes"], mass=payload["modal_mass_kind"]))


def expect_error(label, fn, contains):
    try:
        fn()
    except Exception as e:      # stdlib: ValueError; FastAPI: HTTPException 由呼叫端轉換
        assert contains in str(e), f"{label}: 錯誤訊息應包含「{contains}」, 實際: {e}"
        print(f"  {label}: 錯誤訊息 OK ({str(e)[:44]}...)")
        return
    raise AssertionError(f"{label}: 應該要報錯")



def main():
    # ------------------------------------------------------------------
    print("=== (1) stdlib 後端 ===")
    try:
        from webapi_stdlib.server import _modal_payload
    except Exception as e:
        print(f"SKIPPED: 無法匯入 webapi_stdlib.server ({type(e).__name__}: {e})")
        return

    payload = base_payload()
    out = _modal_payload(payload)
    want = core_expected(payload)
    assert json.loads(json.dumps(out, allow_nan=False)) == json.loads(json.dumps(want)), "/modal 結果應與核心直接算的完全相同"
    print(f"  /modal 與核心直接算逐項相同: {out['n_modes']} 個模態, T = {[round(m['period'], 5) for m in out['modes']]}")
    assert abs(out['mass_total']['x'] - (5000 + 5000 + 7850 * 1e-2 * (4 + 6 + 4))) < 1e-6, "總質量 = 節點質量 + 桿件質量"
    # 沒有 Iz 欄位、Iz=None 的舊格式不受影響
    old_style = base_payload(nodes=[{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": 4, "mx": 1000},
                                     {"id": 2, "x": 6, "y": 4, "my": 1000, "Iz": None}, {"id": 3, "x": 6, "y": 0}])
    _modal_payload(old_style)
    print("  Iz 缺欄位或為 None 的節點 OK")

    bad = base_payload(supports=[])
    expect_error("無支承", lambda: _modal_payload(bad), "機構")
    bad = base_payload(sections=[{"name": "s", "E": 200e9, "I": 8e-5, "A": 1e-2}],
                       nodes=[{"id": i, "x": x, "y": y} for i, (x, y) in enumerate([(0, 0), (0, 4), (6, 4), (6, 0)])])
    expect_error("完全沒有質量", lambda: _modal_payload(bad), "沒有任何有質量")
    bad = base_payload(modal_mass_kind="diagonal")
    expect_error("非法質量種類", lambda: _modal_payload(bad), "lumped")
    bad = base_payload(supports=[{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": -0.01, "rot": 0}])
    expect_error("非零支座沉陷", lambda: _modal_payload(bad), "非零指定位移")
    bad = base_payload(equal_dofs=[{"master_node": 1, "slave_node": 2, "ux": True, "uy": False, "rot": False}])
    expect_error("equalDOF", lambda: _modal_payload(bad), "equal_dof")
    bad = base_payload()
    bad["sections"].append({"name": "c", "E": 200e9, "I": 8e-5, "A": 1e-3, "rho": 7850})
    bad["members"].append({"id": 9, "node_i": 0, "node_j": 2, "section": "c", "member_type": "cable"})
    expect_error("cable", lambda: _modal_payload(bad), "cable")

    # ------------------------------------------------------------------
    print("=== (2) FastAPI 後端 ===")
    try:
        from fastapi.testclient import TestClient
        from webapi.main import app
    except Exception as e:
        print(f"  SKIPPED(這一段): 沒有 fastapi/httpx ({type(e).__name__}), 略過 FastAPI 後端比對")
        print("\n全部通過(stdlib 後端; FastAPI 段略過)")
        return

    client = TestClient(app)
    r = client.post("/modal", json=payload)
    assert r.status_code == 200, r.text
    assert r.json() == json.loads(json.dumps(out)), "FastAPI 與 stdlib 後端的 /modal 回傳必須完全相同"
    print("  FastAPI 與 stdlib 後端 /modal 回傳完全相同")
    for label, bad_payload, contains in (
        ("無支承", base_payload(supports=[]), "機構"),
        ("非零支座沉陷", base_payload(supports=[{"node": 0, "ux": 0, "uy": 0, "rot": 0},
                                              {"node": 3, "ux": 0, "uy": -0.01, "rot": 0}]), "非零指定位移"),
    ):
        rr = client.post("/modal", json=bad_payload)
        assert rr.status_code == 400 and contains in rr.json()["detail"], f"{label}: 應回400含「{contains}」: {rr.text[:80]}"
        print(f"  {label}: 400 OK")
    assert client.post("/modal", json=base_payload(modal_mass_kind="diagonal")).status_code == 422, "非法質量種類應被schema擋下(422)"
    solve_payload = {k: v for k, v in payload.items() if not k.startswith("modal")}
    solve_payload["point_loads"] = [{"node": 1, "fx": 10000, "fy": 0, "m": 0}]
    rs = client.post("/solve", json=solve_payload)
    assert rs.status_code == 200 and abs(rs.json()["nodes"][1]["ux"]) > 0, "/solve 不應受到質量欄位影響"
    print("  /solve 帶質量欄位照常求解 OK")
    print("\n全部通過: stdlib 與 FastAPI 兩個後端的 /modal 行為一致、錯誤訊息清楚、舊呼叫端不受影響")


if __name__ == "__main__":
    main()
