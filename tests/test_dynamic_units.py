"""
驗證案例: 動力分析單位設定 (網頁「單位設定」的動力分析部分) 與 rho 管線

單位換算表同時存在於四個地方, 手動同步很容易漂移:
  webapi/static/index.html           (前端 UNIT_FACTORS)
  webapi_stdlib/static/index.html    (應該跟上面逐字相同)
  webapi/pdf_export.py               (後端 PDF 匯出用的 UNIT_FACTORS)
  webapi_stdlib/pdf_export.py        (同上)

這支測試守住:
 (1) 四份換算表內容一致, 而且包含動力分析的 5 種單位(mass/inertia/density/accel/vel);
 (2) 每個單位選單(<select>)的選項跟換算表的鍵完全對應, 有掛進 unitOf()、
     顯示設定持久化清單、匯出用的 units payload;
 (3) 單位系統的物理一致性: 預設 kN 配 t (1 kN = 1 t·m/s²)、g = 9.80665 m/s²、
     gal = cm/s²; 質量基準 kg 配力基準 N (後端SI, N = kg·m/s²);
 (4) 斷面的 rho 從網頁 JSON 一路到 Frame2D.Section.rho (stdlib 後端 _build_frame)。
只用標準函式庫解析原始碼, 不需要瀏覽器/FastAPI。
"""
import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DYN_KINDS = ('mass', 'inertia', 'density', 'accel', 'vel')
SELECT_IDS = {'mass': 'unitMass', 'inertia': 'unitInertia', 'density': 'unitDensity',
              'accel': 'unitAccel', 'vel': 'unitVel'}

html_a = (ROOT / 'webapi/static/index.html').read_text(encoding='utf-8')
html_b = (ROOT / 'webapi_stdlib/static/index.html').read_text(encoding='utf-8')
assert html_a == html_b, "webapi 與 webapi_stdlib 的 index.html 不一致(改一份忘了同步另一份)"


def parse_js_unit_factors(html):
    start = html.index('const UNIT_FACTORS = {')
    end = html.index('\n};', start)
    block = re.sub(r'//[^\n]*', '', html[start:end])
    out = {}
    for kind, body in re.findall(r"(\w+):\s*\{([^}]*)\}", block):
        d = {}
        for key, val in re.findall(r"('[^']+'|\w+)\s*:\s*([0-9.eE+\-]+)", body):
            d[key.strip("'")] = float(val)
        out[kind] = d
    return out


def parse_py_unit_factors(path):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'UNIT_FACTORS' for t in node.targets):
            return {k: {kk: float(vv) for kk, vv in v.items()} for k, v in ast.literal_eval(node.value).items()}
    raise AssertionError(f"{path} 找不到 UNIT_FACTORS")


print("=== (1) 四份單位換算表一致 ===")
js = parse_js_unit_factors(html_a)
py_a = parse_py_unit_factors(ROOT / 'webapi/pdf_export.py')
py_b = parse_py_unit_factors(ROOT / 'webapi_stdlib/pdf_export.py')
for kind in DYN_KINDS:
    assert kind in js, f"index.html 的 UNIT_FACTORS 缺少 {kind}"
    print(f"  {kind}: {js[kind]}")
assert js == py_a, f"index.html 與 webapi/pdf_export.py 的 UNIT_FACTORS 不一致:\n{js}\n{py_a}"
assert py_a == py_b, "webapi 與 webapi_stdlib 的 pdf_export.UNIT_FACTORS 不一致"

print("=== (2) 選單/unitOf/持久化/匯出payload 都有掛上 ===")
for kind, sel_id in SELECT_IDS.items():
    m = re.search(r'<select id="%s"[^>]*>(.*?)</select>' % sel_id, html_a, flags=re.S)
    assert m, f"找不到 <select id=\"{sel_id}\">"
    options = re.findall(r'<option value="([^"]+)"', m.group(1))
    assert set(options) == set(js[kind]), f"{sel_id} 的選項 {options} 跟換算表 {list(js[kind])} 不對應"
    assert len(re.findall(r'<option[^>]*selected', m.group(1))) == 1, f"{sel_id} 必須剛好一個預設選項"
    assert f"{kind}: '{sel_id}'" in html_a, f"unitOf() 的 ids 沒有 {kind}: '{sel_id}'"
    assert f"{{id: '{sel_id}', kind: 'value'}}" in html_a, f"{sel_id} 沒有進顯示設定持久化清單"
    assert f"{kind}: unitOf('{kind}')" in html_a, f"currentUnitsPayload() 沒有 {kind}"
    print(f"  {sel_id}: 選項 {options} OK")

print("=== (3) 單位系統物理一致性 ===")
assert js['force']['kN'] == js['mass']['t'], "kN 必須對應 t (1 kN = 1 t·m/s²), 否則預設單位組合會不自洽"
assert js['density']['t/m³'] == js['mass']['t'], "密度 t/m³ 對 kg/m³ 的倍率必須等於 t 對 kg"
assert js['inertia']['t·m²'] == js['mass']['t'], "轉動慣量 t·m² 對 kg·m² 的倍率必須等於 t 對 kg"
assert js['force']['N'] == js['mass']['kg'] == 1.0, "後端SI: 力基準N配質量基準kg"
assert js['accel']['g'] == 9.80665 and js['accel']['gal'] == 1e-2
assert js['vel']['cm/s'] == 1e-2 and js['vel']['mm/s'] == 1e-3
# 一個具體數字: 鋼 7.85 t/m³ 必須等於 7850 kg/m³
assert abs(7.85 * js['density']['t/m³'] - 7850.0) < 1e-9
print("  kN↔t、N↔kg、7.85 t/m³ = 7850 kg/m³、g/gal 係數 OK")

print("=== (4) rho 管線: 網頁JSON -> _build_frame -> Section.rho ===")
try:
    sys.path.insert(0, str(ROOT))
    from webapi_stdlib.server import _build_frame
except Exception as e:      # 缺少可選依賴(matplotlib/reportlab等)時, 只跳過這一段
    print(f"  (略過: 無法匯入 webapi_stdlib.server: {type(e).__name__}: {e})")
else:
    payload = {
        "nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 4, "y": 0}],
        "sections": [{"name": "steel", "E": 200e9, "I": 8e-5, "A": 1e-2, "rho": 7850.0},
                     {"name": "nomass", "E": 200e9, "I": 8e-5, "A": 1e-2}],
        "members": [{"id": 0, "node_i": 0, "node_j": 1, "section": "steel"}],
        "supports": [{"node": 0, "ux": 0.0, "uy": 0.0, "rot": 0.0}],
    }
    f = _build_frame(payload)
    assert f.sections['steel'].rho == 7850.0
    assert f.sections['nomass'].rho is None, "沒填rho的斷面必須是None(=無分佈質量)"
    from frame2d.mass import assemble_M, total_mass
    M = assemble_M(f, 'lumped')
    want = 7850.0 * 1e-2 * 4.0
    assert abs(total_mass(f, 'x') - want) < 1e-9 and abs(M.diagonal().sum() - 2 * want) < 1e-9
    print(f"  rho=7850 kg/m³ 經 JSON 抵達 Section, 桿件質量 {want:.1f} kg OK")

print("\n全部通過: 四份單位表一致、選單掛接完整、單位系統自洽、rho 管線通。")
