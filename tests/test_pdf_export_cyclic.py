"""
驗證案例: webapi.pdf_export -- 循環(遲滯)分析 PDF 匯出 (PDF匯出D11的延伸: 循環)

跟 test_pdf_export_modal.py/test_pdf_export_rsa.py 同一個驗證哲學: 不依賴任何 PDF 解析套件,
直接檢查 `_cyclic_loop_fig()`/`_cyclic_hinge_fig()`/`build_cyclic_result_data_page()` 回傳的
matplotlib Figure 物件本身。

層1 `_hinge_label_en()`: `cyclic_to_dict()`的塑鉸label是中文(含「端」字, matplotlib預設字型
    不支援, 會變缺字框), 這裡驗證PDF報告改用的英文label(Member N, end i/j)格式正確, 且用
    `warnings.simplefilter('error')`確認產生圖表過程**不會**觸發字型缺字警告(這是實際發生
    過的真實bug, 見下面的說明)。
層2 `_cyclic_loop_fig()`: 迴圈曲線座標點精確等於 `cyclic_to_dict()` 的 u/F(含單位換算);
    降伏事件的三角形標記位置精確等於 events 裡對應的 u/F。
層3 `build_cyclic_result_data_page()`: 耗能表跟降伏順序表逐格比對。
層4 `build_cyclic_pdf_report()`: 有效的PDF bytes、分頁公式(曲線1+結果表1+ceil(塑鉸數/6)+
    質量1+輸入資料3)。
層5 網頁後端(選用, 需要pymupdf): stdlib 與 FastAPI 兩後端文字內容逐頁完全相同;
    analysis_type='seismic' 回應清楚的「還在做」。
"""
import re
import warnings

import numpy as np

from frame2d import Frame2D
from frame2d.cyclic import cyclic_analysis, cyclic_to_dict
from webapi.pdf_export import (_cyclic_hinge_fig, _cyclic_loop_fig, _hinge_label_en,
                               build_cyclic_pdf_report, build_cyclic_result_data_page)

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0


def portal(mp=60e3):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    for mid, (i, j) in enumerate([(0, 1), (1, 2), (3, 2)]):
        f.add_member(mid, node_i=i, node_j=j, section='H', Mp_i=mp, Mp_j=mp, R_post_yield_i=1000.0, R_post_yield_j=1000.0)
    f.fix(0).fix(3)
    return f


def cyclic_result(**over):
    kw = dict(control_nodes=[1], weights=None, direction='x',
             amplitudes=[0.005, 0.01, 0.015], n_cycles=2, step=0.001)
    kw.update(over)
    res = cyclic_analysis(portal(), kw['control_nodes'], kw['weights'], kw['direction'],
                          kw['amplitudes'], kw['n_cycles'], kw['step'])
    return cyclic_to_dict(res, kw['control_nodes'], kw['direction'], kw['n_cycles'])


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: 中文塑鉸label改用英文, 不觸發字型缺字警告 ===")
r = cyclic_result()
assert len(r["hinges"]) >= 2, "這個案例應該至少有2個塑鉸降伏才測得到重點"
for h in r["hinges"]:
    en = _hinge_label_en(h)
    assert en == f"Member {h['member']}, end {'i' if h['end'] == 0 else 'j'}", f"英文label格式錯誤: {en}"
    assert "端" not in en and all(ord(c) < 128 for c in en), f"英文label不應該含中文字元: {en}"
with warnings.catch_warnings():
    warnings.simplefilter("error")
    fig1 = _cyclic_loop_fig(r, units=None)
    fig2 = _cyclic_hinge_fig(r, r["hinges"], units=None)
    fig3 = build_cyclic_result_data_page(r, units=None)
print(f"  {len(r['hinges'])}個塑鉸的英文label格式正確, 產圖過程沒有觸發任何警告(含字型缺字)")

# ------------------------------------------------------------------
print("=== 層2: _cyclic_loop_fig() 是精確的pass-through ===")
fig = _cyclic_loop_fig(r, units={'force': 'kN'})
ax = fig.axes[0]
loop_line = [ln for ln in ax.get_lines() if ln.get_marker() in ('', 'None', None) and len(ln.get_xdata()) > 5][0]
want_u = np.array(r["u"])
want_F = np.array(r["F"]) / 1000.0
assert np.allclose(loop_line.get_xdata(), want_u), "迴圈u座標應該精確等於cyclic_to_dict()的u"
assert np.allclose(loop_line.get_ydata(), want_F), "迴圈F座標應該精確換算成kN"
scatter = [c for c in ax.collections]
assert len(scatter) == 1, "應該有1組降伏事件散佈點"
yield_events = [e for e in r["events"] if e["kind"] == "yield"]
offsets = scatter[0].get_offsets()
assert len(offsets) == len(yield_events), f"降伏三角形數量應該等於降伏事件數{len(yield_events)}, 實際{len(offsets)}"
for (xu, yf), e in zip(offsets, yield_events):
    assert abs(xu - e["u"]) < 1e-9 and abs(yf - e["F"] / 1000.0) < 1e-6
print(f"  迴圈曲線{len(want_u)}點與{len(yield_events)}個降伏事件三角形座標都精確對照 cyclic_to_dict() 的資料")

# ------------------------------------------------------------------
print("=== 層3: build_cyclic_result_data_page() 表格內容 ===")
page = build_cyclic_result_data_page(r, units={'force': 'kN', 'moment': 'kN·m'})
loop_table = page.axes[0].tables[0]
lcells = loop_table.get_celld()
for i, loop in enumerate(r["loops"]):
    row = i + 1
    got_amp = float(lcells[(row, 0)].get_text().get_text())
    assert rel(got_amp, loop["amplitude"]) < 2e-6
    if loop["energy"] is not None:
        got_e = float(lcells[(row, 2)].get_text().get_text())
        assert rel(got_e, loop["energy"] / 1000.0) < 1e-4, f"耗能欄位應該正確換算成kN·m: {got_e} vs {loop['energy']/1000}"
hinge_table = page.axes[1].tables[0]
hcells = hinge_table.get_celld()
for i, h in enumerate(r["hinges"]):
    row = i + 1
    assert hcells[(row, 1)].get_text().get_text() == _hinge_label_en(h)
    assert hcells[(row, 4)].get_text().get_text() == str(h["n_yield"])
print(f"  耗能表{len(r['loops'])}列與降伏順序表{len(r['hinges'])}列逐格比對都正確")

# ------------------------------------------------------------------
print("=== 層4: build_cyclic_pdf_report() 有效性與頁數 ===")
pdf_bytes = build_cyclic_pdf_report(portal(), r, units=None)
assert pdf_bytes[:4] == b'%PDF'
n_pages = len(re.findall(rb'/Type\s*/Page[^s]', pdf_bytes))
want_pages = 1 + 1 + -(-len(r["hinges"]) // 6) + 1 + 3
print(f"  PDF bytes={len(pdf_bytes)}, 頁數={n_pages}, 預期={want_pages}")
assert n_pages == want_pages

# build_cyclic_pdf_report() 內部呼叫 _cyclic_hinge_fig() 用的分頁區間: 這個案例只有2個塑鉸,
# 不管分頁大小是5還是6都只會產生1頁, 頁數本身測不出來——用monkeypatch監聽實際傳入的塑鉸清單
# 才抓得到(跟test_pdf_export_modal.py的分頁驗證同一個手法)。
import webapi.pdf_export as pdf_export_mod
captured_chunks = []
orig_fn = pdf_export_mod._cyclic_hinge_fig
def _spy(cyclic_result_, hinges_subset_, units=None, ncols=2):
    captured_chunks.append([h["member"] for h in hinges_subset_])
    return orig_fn(cyclic_result_, hinges_subset_, units=units, ncols=ncols)
pdf_export_mod._cyclic_hinge_fig = _spy
try:
    build_cyclic_pdf_report(portal(), r, units=None)
finally:
    pdf_export_mod._cyclic_hinge_fig = orig_fn
assert captured_chunks == [[h["member"] for h in r["hinges"]]], \
    f"應該一頁裝下全部{len(r['hinges'])}個塑鉸(分頁大小6), 實際{captured_chunks}"
print(f"  build_cyclic_pdf_report() 內部分頁精確是 {captured_chunks}")

# 這個模型最多只有6個塑鉸(3根桿件×2端), 分頁大小5跟6在這裡剛好看不出差異(切片超出長度時
# 兩者結果一樣)——用合成的假塑鉸清單(複製既有塑鉸的資料結構, 只換member編號)湊出7個, 才真的
# 測得到「分頁大小是不是精確6」這件事。
synthetic_hinges = [dict(r["hinges"][0], member=100 + k) for k in range(7)]
r_big = dict(r, hinges=synthetic_hinges)
captured_chunks2 = []
def _spy2(cyclic_result_, hinges_subset_, units=None, ncols=2):
    captured_chunks2.append([h["member"] for h in hinges_subset_])
    return orig_fn(cyclic_result_, hinges_subset_, units=units, ncols=ncols)
pdf_export_mod._cyclic_hinge_fig = _spy2
try:
    build_cyclic_pdf_report(portal(), r_big, units=None)
finally:
    pdf_export_mod._cyclic_hinge_fig = orig_fn
assert captured_chunks2 == [[100, 101, 102, 103, 104, 105], [106]], \
    f"7個(合成)塑鉸應該精確分成前6個、後1個兩頁, 實際{captured_chunks2}"
print(f"  7個塑鉸(合成資料, 純測分頁邏輯): 分頁精確是 {captured_chunks2}")

# ------------------------------------------------------------------
print("=== 層5: 網頁後端(選用, 需要pymupdf+fastapi) ===")
try:
    import pymupdf
    import sys
    sys.path.insert(0, '.')
    from fastapi.testclient import TestClient
    from webapi.main import app
    from webapi_stdlib.server import _export_cyclic_pdf
except ImportError as e:
    print(f"SKIPPED: 缺少選用套件({e}), 略過後端文字內容比對")
else:
    model = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": H}, {"id": 2, "x": W, "y": H}, {"id": 3, "x": W, "y": 0}],
        "sections": [{"name": "H", "E": E, "I": I, "A": A}],
        "members": [{"id": 0, "node_i": 0, "node_j": 1, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                    "R_post_yield_i": 1000.0, "R_post_yield_j": 1000.0},
                    {"id": 1, "node_i": 1, "node_j": 2, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                    "R_post_yield_i": 1000.0, "R_post_yield_j": 1000.0},
                    {"id": 2, "node_i": 3, "node_j": 2, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                    "R_post_yield_i": 1000.0, "R_post_yield_j": 1000.0}],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "cyclic_control_nodes": [1], "cyclic_weights": None, "cyclic_direction": "x",
        "cyclic_amplitudes": [0.005, 0.01, 0.015], "cyclic_n_cycles": 2, "cyclic_step": 0.001,
        "units": {"force": "kN", "moment": "kN·m", "mass": "t"},
    }
    stdlib_bytes = _export_cyclic_pdf(model)
    client = TestClient(app)
    m2 = dict(model, analysis_type="cyclic")
    resp = client.post("/export/pdf", json=m2)
    assert resp.status_code == 200
    d1 = pymupdf.open(stream=stdlib_bytes, filetype="pdf")
    d2 = pymupdf.open(stream=resp.content, filetype="pdf")
    assert d1.page_count == d2.page_count
    for i in range(d1.page_count):
        assert d1[i].get_text() == d2[i].get_text(), f"第{i}頁文字內容應該完全相同"
    print(f"  stdlib 與 FastAPI 後端產生的 PDF({d1.page_count}頁)文字內容逐頁完全相同")

    rr = client.post("/export/pdf", json=dict(model, analysis_type="seismic"))
    assert rr.status_code == 400 and "非線性地震" in rr.json()["detail"] and "還在做" in rr.json()["detail"]
    print("  seismic 的 PDF 匯出回應清楚的「還在做」訊息")

print("\n全部通過: 循環PDF匯出的迴圈圖/塑鉸圖/結果表/頁數/(選用)後端一致性都正確。")
