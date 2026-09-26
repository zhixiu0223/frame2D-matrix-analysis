"""
驗證案例: webapi.pdf_export -- 模態分析 PDF 匯出 (網頁 D10 的延伸: PDF 匯出, 模態先做)

這裡**不依賴任何 PDF 解析套件**(pymupdf 等不是這個專案既有的依賴, 不想額外引入)——驗證手法
是直接檢查 `_mode_shape_fig()`/`build_mass_data_page()` 回傳的 matplotlib Figure 物件本身
(座標軸標題、線段資料點、表格儲存格文字), 這些內容之後會被 `PdfPages` 原封不動存成 PDF 頁面,
檢查 Figure 本身等於檢查了最終 PDF 會長什麼樣子, 不需要真的解析 PDF 檔案。

層1 `_mode_shape_fig()`: 每個子圖的標題數字(週期/頻率/Γ)、變形曲線的座標點, 都要精確等於
    `modal_to_dict()` 算出來的資料(這一層應該是完全pass-through, 不重算任何東西)。
    分頁邏輯(每頁最多6個模態)用 7 個模態的案例驗證。
層2 `build_mass_data_page()`: 表格儲存格文字對照節點質量/斷面密度的實際數值(含單位換算)。
層3 `build_modal_pdf_report()`: 產生有效的 PDF bytes(以 %PDF 開頭), 頁數符合公式
    1(表格)+ceil(n_modes/6)(振型)+1(質量)+3(輸入資料); 模型完全沒有節點質量/斷面密度時
    也不出錯(表格是空的, 頁數照常)。
層4 網頁後端(`/export/pdf`, analysis_type='modal'): stdlib 與 FastAPI 兩後端萃取出來的
    純文字內容逐頁完全相同(用來對照「同一份資料兩邊算出同一份報告」, 不比較原始bytes——PDF
    裡嵌了產生時間戳記, 兩次呼叫的bytes本來就不會相同, 這不是bug); analysis_type='seismic'
    回應清楚的「還在做」訊息, 不是壞掉的PDF或不明錯誤(rsa/cyclic已經分別在
    test_pdf_export_rsa.py/test_pdf_export_cyclic.py另外測過, 這裡不重複測)。這一層需要
    pymupdf, 沒裝就印 SKIPPED 並正常結束。
"""
import numpy as np

from frame2d import Frame2D
from frame2d.modal import eigen, modal_to_dict
from webapi.pdf_export import _mode_shape_fig, build_mass_data_page, build_modal_pdf_report

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0


def portal(n_extra_nodes=0):
    """門型剛架, 可選擇性地在橫梁上多切幾段(用來湊出更多模態, 測分頁邏輯)。"""
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, H)
    ids = [1]
    n_seg = n_extra_nodes + 1
    for i in range(1, n_seg):
        nid = 100 + i
        f.add_node(nid, W * i / n_seg, H)
        ids.append(nid)
    f.add_node(2, W, H)
    ids.append(2)
    f.add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 'H')
    for i in range(len(ids) - 1):
        f.add_member(10 + i, ids[i], ids[i + 1], 'H')
    f.add_member(2, 3, 2, 'H')
    f.fix(0).fix(3)
    for nid in ids:
        f.add_mass(nid, mx=3000.0, my=3000.0)
    return f


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: _mode_shape_fig() 是精確的pass-through ===")
f = portal()
md = eigen(f, n_modes=4, mass='lumped')
d = modal_to_dict(md)
fig = _mode_shape_fig(f, d, [m["index"] for m in d["modes"]])
assert len(fig.axes) == 4, f"4個模態應該有4個子圖, 實際{len(fig.axes)}"
for ax, m in zip(fig.axes, d["modes"]):
    title = ax.get_title()
    assert f"T={m['period']:.4f}s" in title, f"標題應該包含週期: {title}"
    assert f"f={m['frequency']:.4f}Hz" in title, f"標題應該包含頻率: {title}"
    assert f"Γx={m['gamma_x']:.3f}" in title, f"標題應該包含正確的Γx: {title}"
    assert f"Γy={m['gamma_y']:.3f}" in title, f"標題應該包含正確的Γy: {title}"
    # 找出這個子圖裡「非虛線」的那條線(變形曲線), 逐點比對modal_to_dict()算出來的座標
    deformed_lines = [ln for ln in ax.get_lines() if ln.get_linestyle() != '--']
    all_x = np.concatenate([ln.get_xdata() for ln in deformed_lines])
    all_y = np.concatenate([ln.get_ydata() for ln in deformed_lines])
    want_x = np.concatenate([np.array(c["X"]) for c in m["curves"].values()])
    want_y = np.concatenate([np.array(c["Y"]) for c in m["curves"].values()])
    assert len(all_x) == len(want_x)
    assert np.allclose(sorted(all_x), sorted(want_x)) and np.allclose(sorted(all_y), sorted(want_y)), \
        f"模態{m['index']}的變形曲線座標應該精確等於modal_to_dict()算出來的資料"
print(f"  4個子圖標題與座標點都精確對照 modal_to_dict() 的資料")

print("  分頁: 7個模態應該分成 ceil(7/6)=2 頁 ===")
f7 = portal(n_extra_nodes=4)
md7 = eigen(f7, n_modes=7, mass='lumped')
d7 = modal_to_dict(md7)
all_idx = [m["index"] for m in d7["modes"]]
page1 = _mode_shape_fig(f7, d7, all_idx[:6])
page2 = _mode_shape_fig(f7, d7, all_idx[6:])
assert len(page1.axes) == 6 and len(page2.axes) == 2
print(f"  第1頁6個模態、第2頁1個模態 OK")

# build_modal_pdf_report() 內部呼叫 _mode_shape_fig() 用的分頁區間是不是真的「每6個一頁」:
# 監聽(monkeypatch)實際呼叫時傳進去的 mode_indices, 而不是隔著PDF bytes猜頁數(頁數不變
# 不代表每頁的模態數量真的對——例如分頁大小從6改成5, 對7個模態一樣是2頁, 但內容不對)。
import webapi.pdf_export as pdf_export_mod
captured_chunks = []
orig_fn = pdf_export_mod._mode_shape_fig
def _spy(f_, modal_result_, mode_indices_, ncols=2):
    captured_chunks.append(list(mode_indices_))
    return orig_fn(f_, modal_result_, mode_indices_, ncols=ncols)
pdf_export_mod._mode_shape_fig = _spy
try:
    build_modal_pdf_report(f7, d7, units=None)
finally:
    pdf_export_mod._mode_shape_fig = orig_fn
assert captured_chunks == [[1, 2, 3, 4, 5, 6], [7]], f"7個模態應該精確分成[1..6]、[7]兩頁, 實際{captured_chunks}"
print(f"  build_modal_pdf_report() 內部分頁精確是 {captured_chunks}")

# ------------------------------------------------------------------
print("=== 層2: build_mass_data_page() 表格內容 ===")
f2 = portal()
f2.sections['H'].rho = 7850.0
page = build_mass_data_page(f2, units={'mass': 't', 'density': 't/m³'})
cells = page.axes[0].tables[1] if False else None
# 直接抓兩張表的所有儲存格文字
tables = [ax.tables[0] for ax in page.axes if ax.tables]
assert len(tables) == 2, f"應該有2張表(密度+節點質量), 實際{len(tables)}"
rho_texts = [c.get_text().get_text() for c in tables[0].get_celld().values()]
assert 'H' in rho_texts and '7.85' in rho_texts, f"密度表應該有斷面名稱H跟7850kg/m³換算成7.85t/m³: {rho_texts}"
mass_texts = [c.get_text().get_text() for c in tables[1].get_celld().values()]
assert '1' in mass_texts, f"節點質量表應該有節點1: {mass_texts}"
# 逐格對照, 不是只看「表格裡某處有沒有出現3」——mx/my欄位各自要對, 才抓得到「my漏了單位換算」
# 這種只有一欄壞掉、另一欄剛好還是對的的錯誤
mass_cells = tables[1].get_celld()
n_data_rows = len(f2.node_masses)
for r in range(1, n_data_rows + 1):   # row 0是表頭
    node_id_text = mass_cells[(r, 0)].get_text().get_text()
    mx_text = mass_cells[(r, 1)].get_text().get_text()
    my_text = mass_cells[(r, 2)].get_text().get_text()
    nm = next(x for x in f2.node_masses if str(x.node) == node_id_text)
    assert mx_text == '3', f"節點{node_id_text}的mx欄位應該是3(t): 實際{mx_text}"
    assert my_text == '3', f"節點{node_id_text}的my欄位應該是3(t): 實際{my_text}"
print(f"  密度表(斷面H, ρ=7.85 t/m³)跟節點質量表(節點1, mx=my=3t)都正確")

# ------------------------------------------------------------------
print("=== 層3: build_modal_pdf_report() 有效性與頁數 ===")
for n_modes, extra in [(4, 0), (7, 4)]:
    f3 = portal(n_extra_nodes=extra)
    md3 = eigen(f3, n_modes=n_modes, mass='lumped')
    d3 = modal_to_dict(md3)
    pdf_bytes = build_modal_pdf_report(f3, d3, units=None)
    assert pdf_bytes[:4] == b'%PDF', "應該是有效的PDF檔案"
    import re
    n_pages = len(re.findall(rb'/Type\s*/Page[^s]', pdf_bytes))
    want_pages = 1 + -(-n_modes // 6) + 1 + 3   # 表格 + ceil(n/6)振型頁 + 質量 + 3頁輸入資料
    print(f"  n_modes={n_modes}: PDF bytes={len(pdf_bytes)}, 頁數(regex估計)={n_pages}, 預期={want_pages}")
    assert n_pages == want_pages, f"頁數應該是{want_pages}, 用regex估計是{n_pages}"

f_nomass_rho = Frame2D()
f_nomass_rho.add_node(0, 0, 0); f_nomass_rho.add_node(1, 0, H)
f_nomass_rho.add_section('H', E=E, I=I, A=A)
f_nomass_rho.add_member(0, 0, 1, 'H')
f_nomass_rho.fix(0)
f_nomass_rho.add_mass(1, mx=100.0)   # 至少要有質量eigen()才不會報錯, 但完全沒有斷面密度rho
md_nr = eigen(f_nomass_rho, n_modes=1, mass='lumped')
pdf_nr = build_modal_pdf_report(f_nomass_rho, modal_to_dict(md_nr), units=None)
assert pdf_nr[:4] == b'%PDF', "沒有斷面密度時也應該正常產生PDF(密度表是空的)"
print("  沒有斷面密度(只有節點質量)時也能正常產生PDF")

# ------------------------------------------------------------------
print("=== 層4: 網頁後端(選用, 需要pymupdf+fastapi) ===")
try:
    import pymupdf
    import sys
    sys.path.insert(0, '.')
    from fastapi.testclient import TestClient
    from webapi.main import app
    from webapi_stdlib.server import _export_modal_pdf
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
        "modal_n_modes": 4, "modal_mass_kind": "lumped", "units": {"force": "kN", "moment": "kN·m", "mass": "t"},
    }
    stdlib_bytes = _export_modal_pdf(model)
    client = TestClient(app)
    m2 = dict(model, analysis_type="modal")
    r = client.post("/export/pdf", json=m2)
    assert r.status_code == 200
    d1 = pymupdf.open(stream=stdlib_bytes, filetype="pdf")
    d2 = pymupdf.open(stream=r.content, filetype="pdf")
    assert d1.page_count == d2.page_count
    for i in range(d1.page_count):
        assert d1[i].get_text() == d2[i].get_text(), f"第{i}頁文字內容應該完全相同"
    print(f"  stdlib 與 FastAPI 後端產生的 PDF({d1.page_count}頁)文字內容逐頁完全相同")

    for at, keyword in [("seismic", "非線性地震")]:
        rr = client.post("/export/pdf", json=dict(model, analysis_type=at))
        assert rr.status_code == 400 and keyword in rr.json()["detail"] and "還在做" in rr.json()["detail"], \
            f"{at}: 應該回400且說明還在做: {rr.status_code} {rr.text[:100]}"
    print("  seismic 的 PDF 匯出回應清楚的「還在做」訊息(不是壞掉的PDF)")

print("\n全部通過: 模態PDF匯出的振型圖/質量表/頁數/(選用)後端一致性都正確。")
