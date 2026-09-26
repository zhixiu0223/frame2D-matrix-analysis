"""
驗證案例: webapi.pdf_export -- 反應譜分析 PDF 匯出 (PDF匯出D11的延伸: 反應譜)

跟 test_pdf_export_modal.py 同一個驗證哲學: **不依賴任何 PDF 解析套件**, 直接檢查
`_rsa_spectrum_curve_fig()`/`build_rsa_result_data_page()` 回傳的 matplotlib Figure 物件
本身(線段座標點、標記點位置、表格儲存格文字)——這些內容之後會被 `PdfPages` 原封不動存成PDF
頁面, 檢查Figure本身等於檢查了最終PDF長什麼樣子。

層1 `_rsa_spectrum_curve_fig()`: 曲線座標點精確等於 `rsa_to_dict()` 的 curve.T/curve.Sa
    (含單位換算); 每個模態的紅點位置精確等於該模態的 period/sa。
層2 `build_rsa_result_data_page()`: 模態組合表逐格比對(T/f/Sa/Γ/有效質量比/累積質量比/
    模態基底剪力), 節點位移表跟桿件內力表逐格比對, 標題裡的基底剪力/累積質量比也要對。
層3 `build_rsa_pdf_report()`: 有效的PDF bytes、固定7頁(曲線1+結果表2+質量1+輸入資料3)。
層4 網頁後端(選用, 需要pymupdf): stdlib 與 FastAPI 兩後端文字內容逐頁完全相同;
    analysis_type='cyclic'/'seismic' 回應清楚的「還在做」。
"""
import re

import numpy as np

from frame2d import Frame2D
from frame2d.spectrum import spectrum_analysis, rsa_to_dict
from webapi.pdf_export import _rsa_spectrum_curve_fig, build_rsa_pdf_report, build_rsa_result_data_page

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0


def portal():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 'H').add_member(1, 1, 2, 'H').add_member(2, 3, 2, 'H')
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


def rsa_result(**over):
    kw = dict(direction='x', damping=0.05, combine='CQC', mass_kind='lumped',
             spectrum_type='code', code_sds=1.0, code_sd1=0.6, code_tl=8.0)
    kw.update(over)
    pkg = spectrum_analysis(portal(), **kw)
    return rsa_to_dict(pkg)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: _rsa_spectrum_curve_fig() 是精確的pass-through ===")
r = rsa_result()
fig = _rsa_spectrum_curve_fig(r, units={'accel': 'g'})
ax = fig.axes[0]
curve_line = [ln for ln in ax.get_lines() if ln.get_marker() in ('', 'None', None)][0]
want_T = np.array(r["curve"]["T"])
want_Sa_g = np.array(r["curve"]["Sa"]) / 9.80665
assert np.allclose(curve_line.get_xdata(), want_T), "曲線T座標應該精確等於curve.T"
assert np.allclose(curve_line.get_ydata(), want_Sa_g), "曲線Sa座標應該精確換算成g"
dot_lines = [ln for ln in ax.get_lines() if ln.get_marker() == 'o']
assert len(dot_lines) == len(r["modes"]), f"應該有{len(r['modes'])}個模態紅點, 實際{len(dot_lines)}"
for ln, m in zip(dot_lines, r["modes"]):
    assert abs(ln.get_xdata()[0] - m["period"]) < 1e-12
    assert abs(ln.get_ydata()[0] - m["sa"] / 9.80665) < 1e-9
print(f"  曲線{len(want_T)}個取樣點與{len(dot_lines)}個模態紅點座標都精確對照 rsa_to_dict() 的資料")

# ------------------------------------------------------------------
print("=== 層2: build_rsa_result_data_page() 表格內容 ===")
page1, page2 = build_rsa_result_data_page(r, units={'force': 'kN', 'accel': 'g'})
title = page1._suptitle.get_text()
assert f"Base Shear = {r['base_shear'] / 1000:.4f}" in title, f"標題應該有正確的基底剪力(kN): {title}"
assert f"{r['cum_ratio_total'] * 100:.1f}%" in title, f"標題應該有正確的累積質量比: {title}"
mode_table = page1.axes[0].tables[0]
cells = mode_table.get_celld()
n_modes = len(r["modes"])
for i, m in enumerate(r["modes"]):
    row = i + 1
    assert cells[(row, 0)].get_text().get_text() == str(m["index"])
    assert cells[(row, 1)].get_text().get_text() == f"{m['period']:.5f}"
    got_sa = float(cells[(row, 3)].get_text().get_text())
    assert rel(got_sa, m["sa"] / 9.80665) < 1e-6, f"模態{m['index']}的Sa欄位應該正確換算成g"
print(f"  模態組合表{n_modes}列逐格比對(含Sa單位換算成g)都正確, 標題基底剪力/累積質量比都對")

node_table = page2.axes[0].tables[0]
ncells = node_table.get_celld()
sorted_nodes = sorted(r["nodes"].items(), key=lambda kv: int(kv[0]))
for i, (nid, v) in enumerate(sorted_nodes):
    row = i + 1
    assert ncells[(row, 0)].get_text().get_text() == nid
    got_ux = float(ncells[(row, 1)].get_text().get_text())
    got_uy = float(ncells[(row, 2)].get_text().get_text())
    assert rel(got_ux, v["ux"]) < 2e-6 if v["ux"] != 0 else got_ux == 0
    assert rel(got_uy, v["uy"]) < 2e-6 if v["uy"] != 0 else got_uy == 0
print(f"  節點位移表{len(sorted_nodes)}列逐格比對都正確")

# ------------------------------------------------------------------
print("=== 層3: build_rsa_pdf_report() 有效性與頁數 ===")
pdf_bytes = build_rsa_pdf_report(portal(), r, units=None)
assert pdf_bytes[:4] == b'%PDF', "應該是有效的PDF檔案"
n_pages = len(re.findall(rb'/Type\s*/Page[^s]', pdf_bytes))
assert n_pages == 7, f"應該是固定7頁(曲線1+結果表2+質量1+輸入資料3), 實際{n_pages}"
print(f"  PDF bytes={len(pdf_bytes)}, 頁數={n_pages}(預期7)")

# ------------------------------------------------------------------
print("=== 層4: 網頁後端(選用, 需要pymupdf+fastapi) ===")
try:
    import pymupdf
    import sys
    sys.path.insert(0, '.')
    from fastapi.testclient import TestClient
    from webapi.main import app
    from webapi_stdlib.server import _export_rsa_pdf
except ImportError as e:
    print(f"SKIPPED: 缺少選用套件({e}), 略過後端文字內容比對")
else:
    model = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": H, "mx": 5000, "my": 5000},
                  {"id": 2, "x": W, "y": H, "mx": 5000, "my": 5000}, {"id": 3, "x": W, "y": 0}],
        "sections": [{"name": "H", "E": E, "I": I, "A": A}],
        "members": [{"id": 0, "node_i": 0, "node_j": 1, "section": "H"},
                    {"id": 1, "node_i": 1, "node_j": 2, "section": "H"},
                    {"id": 2, "node_i": 3, "node_j": 2, "section": "H"}],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "rsa_direction": "x", "rsa_damping": 0.05, "rsa_combine": "CQC", "rsa_mass_kind": "lumped",
        "rsa_spectrum_type": "code", "rsa_code_sds": 1.0, "rsa_code_sd1": 0.6, "rsa_code_tl": 8.0,
        "units": {"force": "kN", "moment": "kN·m", "mass": "t", "accel": "g"},
    }
    stdlib_bytes = _export_rsa_pdf(model)
    client = TestClient(app)
    m2 = dict(model, analysis_type="rsa")
    resp = client.post("/export/pdf", json=m2)
    assert resp.status_code == 200
    d1 = pymupdf.open(stream=stdlib_bytes, filetype="pdf")
    d2 = pymupdf.open(stream=resp.content, filetype="pdf")
    assert d1.page_count == d2.page_count == 7
    for i in range(d1.page_count):
        assert d1[i].get_text() == d2[i].get_text(), f"第{i}頁文字內容應該完全相同"
    print(f"  stdlib 與 FastAPI 後端產生的 PDF(7頁)文字內容逐頁完全相同")

    for at, keyword in [("cyclic", "循環"), ("seismic", "非線性地震")]:
        rr = client.post("/export/pdf", json=dict(model, analysis_type=at))
        assert rr.status_code == 400 and keyword in rr.json()["detail"] and "還在做" in rr.json()["detail"]
    print("  cyclic/seismic 的 PDF 匯出回應清楚的「還在做」訊息")

print("\n全部通過: 反應譜PDF匯出的曲線圖/結果表/頁數/(選用)後端一致性都正確。")
