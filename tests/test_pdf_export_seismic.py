"""
驗證案例: webapi.pdf_export -- 非線性地震反應分析 PDF 匯出(四種分析類型PDF匯出的最後一種)

跟其他三個 test_pdf_export_*.py 同一個驗證哲學: 不依賴任何 PDF 解析套件, 直接檢查
`_seismic_history_fig()`/`_seismic_hysteresis_fig()`/`_seismic_energy_fig()`/
`build_seismic_result_data_page()` 回傳的 matplotlib Figure 物件本身。

層1 三張時程/遲滯/能量圖: 座標點精確等於 `seismic_to_dict()` 的 t/ground_motion/
    control_disp/equivalent_force/energy(含單位換算); 用
    `warnings.simplefilter('error')` 確認不會觸發字型缺字這類警告(D11c循環PDF踩過這個坑,
    這裡直接照抄同一個防護)。
層2 `build_seismic_result_data_page()`: 摘要表逐格比對(含能量平衡兩個數字要相等這件事本身
    也驗證一次——不是只信任frame2d核心, PDF報告呈現的兩個數字也要真的相等); 降伏塑鉸表
    重用跟 test_pdf_export_cyclic.py 一樣的 `_hinge_label_en()` 格式。
層3 `build_seismic_pdf_report()`: 有效PDF bytes、分頁公式(時程1+遲滯1+能量1+摘要1+
    ceil(塑鉸數/6)+質量1+輸入資料3)、塑鉸圖表重用`_cyclic_hinge_fig()`(用monkeypatch確認
    真的呼叫到它, 不是自己重寫一份)。
層4 網頁後端(選用, 需要pymupdf): stdlib 與 FastAPI 兩後端文字內容逐頁完全相同——四種分析
    類型的PDF匯出到這裡全部做完, 這一層也確認不再有任何「還在做」的分析類型剩下。
"""
import re
import warnings

import numpy as np

from frame2d import Frame2D
from frame2d.seismic import nonlinear_seismic_web_analysis, seismic_to_dict
from webapi.pdf_export import (_seismic_energy_fig, _seismic_history_fig, _seismic_hysteresis_fig,
                               build_seismic_pdf_report, build_seismic_result_data_page)

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0
G = 9.80665


def portal(mp=60e3):
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    for mid, (i, j, L) in enumerate([(0, 1, H), (1, 2, W), (3, 2, H)]):
        R = 0.1 * E * I / L
        f.add_member(mid, node_i=i, node_j=j, section='H', Mp_i=mp, Mp_j=mp, R_post_yield_i=R, R_post_yield_j=R)
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


def seismic_result(**over):
    kw = dict(direction='x', zeta=0.05, ground_motion_type='pulse',
             pulse_amplitude_g=0.5, pulse_freq_hz=1.5, pulse_decay=0.25)
    kw.update(over)
    pkg = nonlinear_seismic_web_analysis(portal(), 1, **kw)
    return seismic_to_dict(pkg)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: 三張圖是精確的pass-through, 不觸發任何警告 ===")
r = seismic_result()
assert len(r["hinges"]) >= 2, "這個案例應該至少有2個塑鉸降伏才測得到重點"
with warnings.catch_warnings():
    warnings.simplefilter("error")
    fig_h = _seismic_history_fig(r, units={'disp': 'mm'})
    fig_y = _seismic_hysteresis_fig(r, units={'force': 'kN'})
    fig_e = _seismic_energy_fig(r, units={'moment': 'kN·m'})

ax1, ax2 = fig_h.axes
line1 = ax1.get_lines()[0]
assert np.allclose(line1.get_xdata(), r["t"])
assert np.allclose(line1.get_ydata(), np.array(r["ground_motion"]) / G)
line2 = ax2.get_lines()[0]
assert np.allclose(line2.get_ydata(), np.array(r["control_disp"]) * 1000), "位移應該正確換算成mm"
print("  地震時程圖(地面加速度+位移)座標點精確對照 seismic_to_dict() 的資料")

axh = fig_y.axes[0]
lineh = [ln for ln in axh.get_lines() if len(ln.get_xdata()) > 5][0]
assert np.allclose(lineh.get_xdata(), r["control_disp"])
assert np.allclose(lineh.get_ydata(), np.array(r["equivalent_force"]) / 1000.0), "力應該正確換算成kN"
print("  遲滯圖座標點精確對照 control_disp/equivalent_force(含kN換算)")

axe = fig_e.axes[0]
want_wext = np.array(r["energy"]["Wext"]) / 1000.0
line_wext = [ln for ln in axe.get_lines() if ln.get_linestyle() == '--'][0]
assert np.allclose(line_wext.get_ydata(), want_wext), "外力作功虛線應該正確換算成kN·m"
print("  能量平衡圖的外力作功虛線精確對照 energy.Wext(含kN·m換算)")

# ------------------------------------------------------------------
print("=== 層2: build_seismic_result_data_page() 表格內容 ===")
page = build_seismic_result_data_page(r, units={'force': 'kN', 'moment': 'kN·m', 'disp': 'mm'})
summary_table = page.axes[0].tables[0]
scells = summary_table.get_celld()
n_rows = 9
summary = {scells[(i, 0)].get_text().get_text(): scells[(i, 1)].get_text().get_text() for i in range(1, n_rows + 1)}
e_final = r["energy"]["Wext"][-1] / 1000.0
e_check = (r["energy"]["KE"][-1] + r["energy"]["Wdamp"][-1] + r["energy"]["Wint"][-1]) / 1000.0
got_ext = float(summary["Energy balance: external work (kN·m)"])
got_check = float(summary["Energy balance: KE+damping+internal (kN·m, should match above exactly)"])
assert rel(got_ext, e_final) < 1e-5 and rel(got_check, e_check) < 1e-5
assert abs(got_ext - got_check) < 1e-6 * max(abs(got_ext), 1), \
    f"PDF報告裡的能量平衡兩個數字本身應該相等: {got_ext} vs {got_check}"
n_yield = sum(1 for e in r["events"] if e["kind"] == "yield")
got_yield = summary["Yield / unload events"].split(" / ")[0]
assert got_yield == str(n_yield)
print(f"  摘要表逐項比對正確, 能量平衡兩個數字在報告裡本身相等({got_ext} vs {got_check})")

# 用一個「刻意讓Wext跟KE+Wdamp+Wint不相等」的假資料測試, 這樣才能真的分清楚表格第二列是不是
# 真的來自KE+Wdamp+Wint、不是不小心又抄了一次Wext——正常案例下這兩個值本來就該相等(能量平衡
# 成立), 用真實資料測不出「抄錯來源但剛好數字一樣」這種錯誤。
import copy
r_fake = copy.deepcopy(r)
r_fake["energy"]["Wext"][-1] = 999.0
r_fake["energy"]["KE"][-1] = 1.0
r_fake["energy"]["Wdamp"][-1] = 2.0
r_fake["energy"]["Wint"][-1] = 3.0   # KE+Wdamp+Wint = 6.0, 刻意跟Wext(999.0)不同
page_fake = build_seismic_result_data_page(r_fake, units=None)
scells_fake = page_fake.axes[0].tables[0].get_celld()
summary_fake = {scells_fake[(i, 0)].get_text().get_text(): scells_fake[(i, 1)].get_text().get_text() for i in range(1, n_rows + 1)}
got_ext_fake = float(summary_fake["Energy balance: external work (N·m)"])
got_check_fake = float(summary_fake["Energy balance: KE+damping+internal (N·m, should match above exactly)"])
assert abs(got_ext_fake - 999.0) < 1e-6, f"外力作功欄位應該是999.0, 實際{got_ext_fake}"
assert abs(got_check_fake - 6.0) < 1e-6, f"KE+damping+internal欄位應該是1+2+3=6.0(不是999.0), 實際{got_check_fake}"
print(f"  用刻意不相等的假資料驗證兩個欄位分別來自正確的來源(999.0 跟 6.0, 不是兩個都999.0)")

hinge_table = page.axes[1].tables[0]
hcells = hinge_table.get_celld()
for i, h in enumerate(r["hinges"]):
    row = i + 1
    assert hcells[(row, 0)].get_text().get_text() == f"Member {h['member']}, end {'i' if h['end'] == 0 else 'j'}"
    assert hcells[(row, 1)].get_text().get_text() == str(h["n_yield"])
print(f"  降伏塑鉸表{len(r['hinges'])}列逐格比對都正確")

# ------------------------------------------------------------------
print("=== 層3: build_seismic_pdf_report() 有效性、頁數、重用_cyclic_hinge_fig() ===")
import webapi.pdf_export as pdf_export_mod
captured = []
orig_fn = pdf_export_mod._cyclic_hinge_fig
def _spy(cyclic_result_, hinges_subset_, units=None, ncols=2):
    captured.append([h["member"] for h in hinges_subset_])
    return orig_fn(cyclic_result_, hinges_subset_, units=units, ncols=ncols)
pdf_export_mod._cyclic_hinge_fig = _spy
try:
    pdf_bytes = build_seismic_pdf_report(portal(), r, units=None)
finally:
    pdf_export_mod._cyclic_hinge_fig = orig_fn
assert captured == [[h["member"] for h in r["hinges"]]], \
    f"塑鉸圖應該重用_cyclic_hinge_fig(), 監聽到的呼叫參數應該是全部塑鉸: {captured}"
print(f"  塑鉸M-θp圖確實重用了 _cyclic_hinge_fig()(監聽到呼叫參數 {captured})")

assert pdf_bytes[:4] == b'%PDF'
n_pages = len(re.findall(rb'/Type\s*/Page[^s]', pdf_bytes))
want_pages = 4 + -(-len(r["hinges"]) // 6) + 1 + 3   # 時程+遲滯+能量+摘要=4, 塑鉸, 質量, 輸入資料3
print(f"  PDF bytes={len(pdf_bytes)}, 頁數={n_pages}, 預期={want_pages}")
assert n_pages == want_pages

# ------------------------------------------------------------------
print("=== 層4: 網頁後端(選用, 需要pymupdf+fastapi) ===")
try:
    import pymupdf
    import sys
    sys.path.insert(0, '.')
    from fastapi.testclient import TestClient
    from webapi.main import app
    from webapi_stdlib.server import _export_seismic_pdf
except ImportError as e:
    print(f"SKIPPED: 缺少選用套件({e}), 略過後端文字內容比對")
else:
    model = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 0, "y": H, "mx": 5000, "my": 5000},
                  {"id": 2, "x": W, "y": H, "mx": 5000, "my": 5000}, {"id": 3, "x": W, "y": 0}],
        "sections": [{"name": "H", "E": E, "I": I, "A": A}],
        "members": [{"id": 0, "node_i": 0, "node_j": 1, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                    "R_post_yield_i": 0.1 * E * I / H, "R_post_yield_j": 0.1 * E * I / H},
                    {"id": 1, "node_i": 1, "node_j": 2, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                    "R_post_yield_i": 0.1 * E * I / W, "R_post_yield_j": 0.1 * E * I / W},
                    {"id": 2, "node_i": 3, "node_j": 2, "section": "H", "Mp_i": 60000.0, "Mp_j": 60000.0,
                    "R_post_yield_i": 0.1 * E * I / H, "R_post_yield_j": 0.1 * E * I / H}],
        "supports": [{"node": 0, "ux": 0, "uy": 0, "rot": 0}, {"node": 3, "ux": 0, "uy": 0, "rot": 0}],
        "seismic_control_node": 1, "seismic_direction": "x", "seismic_zeta": 0.05,
        "seismic_damping_modes": [0, 2], "seismic_mass_kind": "lumped",
        "seismic_ground_motion_type": "pulse", "seismic_pulse_amplitude_g": 0.5,
        "seismic_pulse_freq_hz": 1.5, "seismic_pulse_decay": 0.25,
        "units": {"force": "kN", "moment": "kN·m", "mass": "t"},
    }
    stdlib_bytes = _export_seismic_pdf(model)
    client = TestClient(app)
    m2 = dict(model, analysis_type="seismic")
    resp = client.post("/export/pdf", json=m2)
    assert resp.status_code == 200
    d1 = pymupdf.open(stream=stdlib_bytes, filetype="pdf")
    d2 = pymupdf.open(stream=resp.content, filetype="pdf")
    assert d1.page_count == d2.page_count
    for i in range(d1.page_count):
        assert d1[i].get_text() == d2[i].get_text(), f"第{i}頁文字內容應該完全相同"
    print(f"  stdlib 與 FastAPI 後端產生的 PDF({d1.page_count}頁)文字內容逐頁完全相同")

    for at in ("modal", "rsa", "cyclic", "seismic"):
        extra = {"modal_n_modes": 4, "modal_mass_kind": "lumped",
                "rsa_direction": "x", "rsa_code_sds": 1.0, "rsa_code_sd1": 0.6,
                "cyclic_control_nodes": [1], "cyclic_amplitudes": [0.005], "cyclic_step": 0.001}
        rr = client.post("/export/pdf", json=dict(model, analysis_type=at, **extra))
        assert rr.status_code == 200, f"{at}: 四種分析類型現在都應該能匯出PDF, 實際{rr.status_code} {rr.text[:150]}"
    print("  確認 modal/rsa/cyclic/seismic 四種分析類型現在都能成功匯出PDF, 沒有任何一種還卡在「還在做」")

print("\n全部通過: 非線性地震PDF匯出的時程圖/遲滯圖/能量圖/摘要表/塑鉸圖(重用)/頁數/(選用)後端一致性都正確。")
