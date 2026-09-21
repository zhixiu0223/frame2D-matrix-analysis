"""
驗證案例: 網頁後端 /cyclic 端點 (動力分析 D7b)

 (1) stdlib 後端: 桿件的塑鉸容量(Mp_i/Mp_j/R_post_yield_*)經 JSON 抵達 frame2d, /cyclic 的結果跟核心直接
     呼叫 cyclic_analysis()+cyclic_to_dict() 逐項相同; 各種錯誤情況(缺欄位、沒有塑鉸容量、預載超過
     Mp、步數過多)都是清楚的錯誤訊息。
 (2) FastAPI 後端(有安裝 fastapi/httpx 才測, 否則略過這一段): 同一個 payload 回傳與 stdlib 完全相同;
     錯誤情況回 400 與相同訊息; /solve 不受新欄位影響。
兩個後端各自有一份 _build_frame 與端點, 這支測試守住它們不漂移。
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from frame2d import Frame2D
from frame2d.cyclic import cyclic_analysis, cyclic_to_dict

# 注意: 不能在模組層級 sys.exit()——pytest 收集測試會 import 這個檔案, SystemExit 會讓整個 pytest
# INTERNALERROR。所有實際工作放在 main() 裡, 只有直接執行(或被 test_zz 的 subprocess 執行)才會跑。

E, I, A = 200e9, 8e-5, 5e-3


def base_payload(**over):
    def member(mid, i, j, L):
        R = 0.1 * E * I / L
        return {"id": mid, "node_i": i, "node_j": j, "section": "H", "Mp_i": 100e3, "Mp_j": 100e3,
                "R_post_yield_i": R, "R_post_yield_j": R}
    p = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": 4}, {"id": 2, "x": 6, "y": 4}, {"id": 3, "x": 6, "y": 0}],
        "sections": [{"name": "H", "E": E, "I": I, "A": A}],
        "members": [member(0, 0, 1, 4.0), member(1, 1, 2, 6.0), member(2, 3, 2, 4.0)],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "cyclic_control_nodes": [1], "cyclic_direction": "x",
        "cyclic_amplitudes": [0.03, 0.06], "cyclic_n_cycles": 2, "cyclic_step": 0.004,
    }
    p.update(over)
    return p


def core_expected(payload):
    f = Frame2D()
    for n in payload["nodes"]:
        f.add_node(n["id"], n["x"], n["y"])
    for s in payload["sections"]:
        f.add_section(s["name"], E=s["E"], I=s["I"], A=s["A"])
    for m in payload["members"]:
        f.add_member(m["id"], node_i=m["node_i"], node_j=m["node_j"], section=m["section"],
                     Mp_i=m.get("Mp_i"), Mp_j=m.get("Mp_j"),
                     R_post_yield_i=m.get("R_post_yield_i"), R_post_yield_j=m.get("R_post_yield_j"))
    for sp in payload["supports"]:
        f.support(sp["node"], ux=sp["ux"], uy=sp["uy"], rot=sp["rot"])
    for pl in payload.get("point_loads", []):
        f.point_load(pl["node"], fx=pl["fx"], fy=pl["fy"], m=pl["m"])
    res = cyclic_analysis(f, payload["cyclic_control_nodes"], payload.get("cyclic_weights"),
                          payload["cyclic_direction"], payload["cyclic_amplitudes"],
                          payload["cyclic_n_cycles"], payload["cyclic_step"])
    return cyclic_to_dict(res, payload["cyclic_control_nodes"], payload["cyclic_direction"], payload["cyclic_n_cycles"])


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
        from webapi_stdlib.server import _cyclic_payload
    except Exception as e:
        print(f"SKIPPED: 無法匯入 webapi_stdlib.server ({type(e).__name__}: {e})")
        return

    payload = base_payload()
    out = _cyclic_payload(payload)
    want = core_expected(payload)
    assert json.loads(json.dumps(out, allow_nan=False)) == json.loads(json.dumps(want)), "/cyclic 結果應與核心直接算的完全相同"
    print(f"  /cyclic 與核心直接算逐項相同: {len(out['u'])} 個點, 塑鉸 {[h['label'] for h in out['hinges']]}, "
          f"迴圈 {[round(l['energy'], 1) if l['energy'] is not None else None for l in out['loops']]}")
    assert len(out['hinges']) >= 2 and out['loops'][-1]['closed'] and out['loops'][-1]['energy'] > 0

    # 有載重(當重力預載)
    pl = base_payload(point_loads=[{"node": 1, "fx": 5000, "fy": 0, "m": 0}])
    out_pl = _cyclic_payload(pl)
    assert out_pl['events'][0]['F'] != out['events'][0]['F'], "重力預載應該改變首次降伏時的底剪力"
    assert json.loads(json.dumps(out_pl)) == json.loads(json.dumps(core_expected(pl)))
    print(f"  重力預載: 首次降伏底剪力 {out['events'][0]['F']:.1f} → {out_pl['events'][0]['F']:.1f} N OK")

    expect_error("缺幅值", lambda: _cyclic_payload(base_payload(cyclic_amplitudes=None)), "cyclic_amplitudes")
    expect_error("缺控制節點", lambda: _cyclic_payload(base_payload(cyclic_control_nodes=[])), "cyclic_control_nodes")
    expect_error("缺步長", lambda: _cyclic_payload(base_payload(cyclic_step=None)), "cyclic_step")
    no_hinge = base_payload(members=[{k: v for k, v in m.items() if k in ("id", "node_i", "node_j", "section")}
                                     for m in base_payload()["members"]])
    expect_error("沒有塑鉸容量", lambda: _cyclic_payload(no_hinge), "塑鉸容量")
    expect_error("幅值非正數", lambda: _cyclic_payload(base_payload(cyclic_amplitudes=[0.03, -0.06])), "正數")
    expect_error("不存在的控制節點", lambda: _cyclic_payload(base_payload(cyclic_control_nodes=[99])), "找不到控制節點")
    expect_error("預載超過Mp", lambda: _cyclic_payload(base_payload(point_loads=[{"node": 1, "fx": 300000, "fy": 0, "m": 0}])), "已超出塑性彎矩")
    expect_error("步數過多", lambda: _cyclic_payload(base_payload(cyclic_step=1e-9)), "預估要走")

    print("=== (2) FastAPI 後端 ===")
    try:
        from fastapi.testclient import TestClient
        from webapi.main import app
    except Exception as e:
        print(f"  SKIPPED(這一段): 沒有 fastapi/httpx ({type(e).__name__}), 略過 FastAPI 後端比對")
        print("\n全部通過(stdlib 後端; FastAPI 段略過)")
        return

    client = TestClient(app)
    r = client.post("/cyclic", json=payload)
    assert r.status_code == 200, r.text
    assert r.json() == json.loads(json.dumps(out)), "FastAPI 與 stdlib 後端的 /cyclic 回傳必須完全相同"
    print("  FastAPI 與 stdlib 後端 /cyclic 回傳完全相同")
    for label, bad, contains in (
        ("缺幅值", base_payload(cyclic_amplitudes=None), "cyclic_amplitudes"),
        ("沒有塑鉸容量", no_hinge, "塑鉸容量"),
        ("預載超過Mp", base_payload(point_loads=[{"node": 1, "fx": 300000, "fy": 0, "m": 0}]), "已超出塑性彎矩"),
    ):
        rr = client.post("/cyclic", json=bad)
        assert rr.status_code == 400 and contains in rr.json()["detail"], f"{label}: 應回400含「{contains}」: {rr.text[:80]}"
        print(f"  {label}: 400 OK")
    assert client.post("/cyclic", json=base_payload(cyclic_direction="z")).status_code == 422, "非法方向應被schema擋下(422)"
    solve_payload = {k: v for k, v in payload.items() if not k.startswith("cyclic")}
    solve_payload["point_loads"] = [{"node": 1, "fx": 10000, "fy": 0, "m": 0}]
    rs = client.post("/solve", json=solve_payload)
    assert rs.status_code == 200 and abs(rs.json()["nodes"][1]["ux"]) > 0, "/solve 不應受到新欄位影響"
    print("  /solve 照常求解 OK")
    print("\n全部通過: stdlib 與 FastAPI 兩個後端的 /cyclic 行為一致、錯誤訊息清楚、舊呼叫端不受影響")


if __name__ == "__main__":
    main()
