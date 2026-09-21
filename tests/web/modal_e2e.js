/*
 * 網頁端到端測試(選用): 模態分析 (D2b)
 *
 * 用 jsdom 載入真正由後端提供的 index.html, 把 window.fetch 接到真的 Python 後端(stdlib),
 * 然後「像使用者一樣」操作介面: 節點屬性面板輸入質量 -> 選「模態」分析 -> 按 Solve ->
 * 檢查結果表/振型/單位切換/錯誤訊息。
 *
 * 用法(由 tests/test_web_modal_e2e.py 呼叫, 一般不需要手動執行):
 *   node tests/web/modal_e2e.js <後端網址> <預期結果JSON>
 * 需要 jsdom:  npm install jsdom
 */
const { JSDOM, VirtualConsole } = require('jsdom');

const base = process.argv[2];
const expected = JSON.parse(process.argv[3]);   // {periods:[...], nodeMassKg:5000}

function assert(cond, msg) { if (!cond) { throw new Error('ASSERT: ' + msg); } }
const sleep = ms => new Promise(r => setTimeout(r, ms));
async function waitFor(fn, what, ms = 8000) {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) { if (fn()) return; await sleep(20); }
  throw new Error('逾時等待: ' + what);
}

(async () => {
  const errs = [];
  const vc = new VirtualConsole();
  vc.on('jsdomError', e => errs.push((e.detail && e.detail.message) || e.message));
  const dom = await JSDOM.fromURL(base + '/', {
    runScripts: 'dangerously', pretendToBeVisual: true, resources: 'usable', virtualConsole: vc,
    beforeParse(window) {
      window.fetch = (u, o) => fetch(new URL(u, base), o);      // 接到真的後端
    },
  });
  const w = dom.window, d = w.document;
  await sleep(300);
  const $ = id => d.getElementById(id);
  const ev = code => w.eval(code);

  // ---- 建模: 門型鋼架(SI: Pa, m, kg), 斷面帶密度 ρ ----
  ev(`
    model.nodes = [{id:0,x:0,y:0},{id:1,x:0,y:4},{id:2,x:6,y:4},{id:3,x:6,y:0}];
    model.sections = [{name:'s',E:200e9,I:8e-5,A:1e-2,rho:7850}];
    model.members = [
      {id:0,node_i:0,node_j:1,section:'s',member_type:'frame',release_i:false,release_j:false},
      {id:1,node_i:1,node_j:2,section:'s',member_type:'frame',release_i:false,release_j:false},
      {id:2,node_i:3,node_j:2,section:'s',member_type:'frame',release_i:false,release_j:false}];
    model.supports = [{node:0,ux:0.0,uy:0.0,rot:0.0},{node:3,ux:0.0,uy:0.0,rot:0.0}];
    model.point_loads = []; model.distributed_loads = [];
    render();`);

  // ---- 節點屬性面板輸入質量(預設質量單位 t) ----
  for (const nid of [1, 2]) {
    ev(`selected = {type:'node', id:${nid}}; renderPropPanel();`);
    assert($('p_mx') && $('p_my') && $('p_Iz'), '節點面板應有質量輸入欄');
    $('p_mx').value = '5'; $('p_my').value = '5'; $('p_Iz').value = '0';
    $('p_apply').click();
  }
  const massOf = nid => ev(`(() => { const n = nodeById(${nid}); return [n.mx, n.my, n.Iz]; })()`);
  const m1 = Array.from(massOf(1));
  assert(Math.abs(m1[0] - expected.nodeMassKg) < 1e-6 && Math.abs(m1[1] - expected.nodeMassKg) < 1e-6,
         `節點質量應該存成 SI 的 ${expected.nodeMassKg} kg, 實際 ${m1}`);
  assert(ev(`nodeById(0).mx === undefined`), '沒填質量的節點不應該有 mx 欄位');
  console.log('節點質量輸入: 5 t -> 存成', m1[0], 'kg OK');

  // 負數質量必須被拒絕、不修改模型
  ev(`selected = {type:'node', id:1}; renderPropPanel();`);
  $('p_mx').value = '-3'; $('p_apply').click();
  assert(ev(`nodeById(1).mx`) === expected.nodeMassKg, '負質量應該被拒絕且不改動原本的質量');
  assert($('status').textContent.includes('不能是負數'), '負質量應該顯示錯誤訊息');
  console.log('負質量被拒絕 OK');
  ev(`clearSelection();`);

  // ---- 選「模態」分析 ----
  const sel = $('analysisType');
  sel.value = 'modal'; sel.dispatchEvent(new w.Event('change'));
  assert($('modalConfigBar').style.display === 'block', '選模態後應顯示模態設定列');
  assert($('pushoverConfigBar').style.display === 'none', 'Pushover設定列應隱藏');
  $('md_nmodes').value = String(expected.periods.length);
  $('md_kind').value = 'consistent';
  $('btnSolve').click();
  await waitFor(() => ev(`modalResult !== null`), '模態分析結果');

  const r = ev(`JSON.parse(JSON.stringify(modalResult))`);
  assert(r.n_modes === expected.periods.length, `模態數 ${r.n_modes}`);
  for (let i = 0; i < r.n_modes; i++) {
    assert(Math.abs(r.modes[i].period / expected.periods[i] - 1) < 1e-9,
           `第${i + 1}模態週期 ${r.modes[i].period} 應等於核心算出的 ${expected.periods[i]}`);
  }
  console.log('週期與核心一致:', r.modes.map(m => m.period.toFixed(5)).join(', '));

  assert(ev(`viewMode`) === 'modal', 'Solve後應自動切到模態形狀視圖');
  assert($('btnViewModal').style.display !== 'none' && $('btnViewModal').classList.contains('active'), '模態形狀按鈕應顯示且啟用');
  assert($('modalPanel').style.display === 'block', '結果面板應顯示');
  let rows = d.querySelectorAll('#modalTable tbody tr');
  assert(rows.length === r.n_modes, `結果表應有 ${r.n_modes} 列, 實際 ${rows.length}`);
  assert(rows[0].classList.contains('active'), '第1列預設選取');
  let lines = d.querySelectorAll('#content polyline.mode-line');
  assert(lines.length === 3, `振型曲線應有 3 條(3根桿件), 實際 ${lines.length}`);
  assert($('overlay').textContent.includes('第 1 模態'), '標題應顯示第 1 模態');
  assert($('status').textContent.includes('模態分析完成'), '狀態列應顯示完成');
  const cellTexts = [...d.querySelectorAll('#modalTable tbody td')].map(td => td.textContent);
  assert(!cellTexts.some(t => /e-(1\d|[2-9]\d)/.test(t)), '結果表不應出現 e-10 以下的數值雜訊(應顯示成 0): ' + cellTexts.filter(t => /e-\d\d/.test(t)).join(','));
  console.log('結果表 3 列、振型曲線 3 條、標題與狀態列、無數值雜訊 OK');

  // ---- 「◀ ▶」切換鈕(不用捲到結果表) ----
  assert($('modalNav').style.display !== 'none', '模態視圖應顯示 ◀ ▶ 切換鈕');
  assert($('modalNavLabel').textContent === `模態 1/${r.n_modes}`, '標籤應顯示 模態 1/N: ' + $('modalNavLabel').textContent);
  assert($('btnModePrev').disabled && !$('btnModeNext').disabled, '第1個模態時「◀」停用、「▶」可用');
  $('btnModeNext').click(); await sleep(30);
  assert(ev(`modalModeIndex`) === 1 && $('overlay').textContent.includes('第 2 模態'), '按 ▶ 應切到第 2 模態');
  $('btnModeNext').click(); await sleep(30);
  assert(ev(`modalModeIndex`) === r.n_modes - 1 && $('overlay').textContent.includes(`第 ${r.n_modes} 模態`), '再按 ▶ 應切到最後一個模態');
  assert($('btnModeNext').disabled, '最後一個模態時「▶」停用');
  $('btnModeNext').click(); await sleep(30);
  assert(ev(`modalModeIndex`) === r.n_modes - 1, '最後一個時再按 ▶ 不應越界');
  $('btnModePrev').click(); await sleep(30);
  assert(ev(`modalModeIndex`) === r.n_modes - 2 && $('modalNavLabel').textContent === `模態 ${r.n_modes - 1}/${r.n_modes}`, '按 ◀ 應回上一個, 標籤同步');
  $('btnModePrev').click(); await sleep(30);
  assert(ev(`modalModeIndex`) === 0 && $('btnModePrev').disabled, '回到第 1 個, ◀ 停用');
  console.log('◀ ▶ 切換鈕(邊界停用、標籤同步、振型跟著換) OK');

  // ---- 點第2列切換振型 ----
  rows[1].click();
  await sleep(50);
  assert(ev(`modalModeIndex`) === 1, '點第2列後 modalModeIndex 應為 1');
  assert($('overlay').textContent.includes('第 2 模態'), '標題應改成第 2 模態');
  rows = d.querySelectorAll('#modalTable tbody tr');
  assert(rows[1].classList.contains('active') && !rows[0].classList.contains('active'), '選取列應跟著換');
  console.log('點選切換模態 OK');

  // ---- 切換質量單位: 表頭與數值要跟著換 ----
  const headBefore = d.querySelector('#modalTable thead').textContent;
  const effBefore = parseFloat(d.querySelectorAll('#modalTable tbody tr')[0].children[4].textContent);
  assert(headBefore.includes('(t)'), '預設質量單位 t');
  $('unitMass').value = 'kg'; $('unitMass').dispatchEvent(new w.Event('change'));
  const headAfter = d.querySelector('#modalTable thead').textContent;
  const effAfter = parseFloat(d.querySelectorAll('#modalTable tbody tr')[0].children[4].textContent);
  assert(headAfter.includes('(kg)'), '切換單位後表頭應顯示 kg');
  assert(Math.abs(effAfter / effBefore - 1000) < 1e-2, `有效質量數值應放大1000倍: ${effBefore} -> ${effAfter}`);
  console.log(`質量單位切換: 有效質量 ${effBefore} t -> ${effAfter} kg OK`);

  // ---- 錯誤路徑: 拿掉支承 -> 機構, 顯示清楚訊息, 舊結果不被破壞 ----
  const supBackup = ev(`JSON.stringify(model.supports)`);
  ev(`model.supports = [];`);
  $('btnSolve').click();
  await waitFor(() => $('status').textContent.includes('模態分析失敗'), '錯誤訊息');
  assert($('status').textContent.includes('機構') || $('status').textContent.includes('支承'), '錯誤訊息應說明機構/支承不足: ' + $('status').textContent);
  assert(ev(`modalResult !== null && modalResult.n_modes`) === r.n_modes, '失敗時不應清掉上一次的結果');
  console.log('機構/支承不足: 顯示清楚訊息 OK ->', $('status').textContent.slice(0, 60));
  ev(`model.supports = ${supBackup};`);

  // ---- 模態數輸入驗證 ----
  $('md_nmodes').value = '0'; $('btnSolve').click();
  await sleep(50);
  assert($('status').textContent.includes('不小於 1'), '模態數 < 1 應被前端擋下');
  console.log('模態數驗證 OK');

  // ---- 新建/載入時重置 ----
  ev(`resetViewToStructure();`);
  assert(ev(`modalResult`) === null && $('btnViewModal').style.display === 'none' && $('modalConfigBar').style.display === 'none',
         '新建/載入檔案時應清掉模態結果與按鈕');
  assert(ev(`viewMode`) === 'structure' && $('analysisType').value === 'linear', '應回到結構視圖與線性分析');
  assert($('modalNav').style.display === 'none', '離開模態視圖後 ◀ ▶ 應隱藏');
  console.log('新建/載入重置 OK');

  // ---- 節點質量標記 ----
  ev(`render();`);
  const marks = Array.from(d.querySelectorAll('#content text')).filter(t => t.textContent === 'm');
  assert(marks.length === 2, `結構視圖應在 2 個有質量的節點旁標 m, 實際 ${marks.length}`);
  console.log('節點質量標記 m: 2 個 OK');

  assert(errs.length === 0, '頁面有 JS 錯誤: ' + errs.slice(0, 3).join(' | '));
  console.log('\n全部通過: 網頁模態分析端到端(質量輸入、Solve、結果表、振型、單位切換、錯誤訊息、重置)');
  process.exit(0);
})().catch(e => { console.error('FAILED:', e.message); process.exit(1); });
