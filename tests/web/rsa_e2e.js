/*
 * 網頁端到端測試(選用): 反應譜分析 (D3b)
 *
 * 用 jsdom 載入真正由後端提供的 index.html, window.fetch 接到真的 Python 後端, 然後「像使用者
 * 一樣」: 匯入範例模型 → 選「反應譜」→ 規範反應譜(SDS/SD1/TL)Solve → 檢查反應譜曲線圖、結果表
 * (模態、節點位移、桿件內力)、單位切換、自訂反應譜、錯誤訊息、重置。網頁拿到的基底剪力/週期/
 * 累積質量比/參與係數跟核心 spectrum_analysis() 直接算的逐項比對(包含渲染出來的表格文字,
 * 不只是底層 JS 資料物件, 這樣才抓得到「資料對但畫錯欄位」這種只發生在渲染層的錯誤)。
 *
 * 用法(由 tests/test_web_rsa_e2e.py 呼叫):
 *   node tests/web/rsa_e2e.js <後端網址> <範例JSON路徑> <預期結果JSON>
 * 需要 jsdom:  npm install jsdom
 */
const { JSDOM, VirtualConsole } = require('jsdom');
const fs = require('fs');

const base = process.argv[2];
const modelFile = process.argv[3];
const expected = JSON.parse(process.argv[4]);   // {baseShear, periods:[...], cumRatioTotal, nModes, gammas:[...]}

function assert(cond, msg) { if (!cond) { throw new Error('ASSERT: ' + msg); } }
const sleep = ms => new Promise(r => setTimeout(r, ms));
async function waitFor(fn, what, ms = 15000) {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) { if (fn()) return; await sleep(25); }
  throw new Error('逾時等待: ' + what);
}

(async () => {
  const errs = [];
  const vc = new VirtualConsole();
  vc.on('jsdomError', e => errs.push((e.detail && e.detail.message) || e.message));
  const dom = await JSDOM.fromURL(base + '/', {
    runScripts: 'dangerously', pretendToBeVisual: true, resources: 'usable', virtualConsole: vc,
    beforeParse(window) { window.fetch = (u, o) => fetch(new URL(u, base), o); },
  });
  const w = dom.window, d = w.document;
  await sleep(300);
  const $ = id => d.getElementById(id);
  const ev = code => w.eval(code);

  // ---- 匯入範例模型 ----
  const input = $('importJsonFile');
  const file = new w.File([fs.readFileSync(modelFile, 'utf8')], 'portal_rsa_demo.json', {type: 'application/json'});
  Object.defineProperty(input, 'files', {value: [file], configurable: true});
  input.dispatchEvent(new w.Event('change'));
  await sleep(300);
  assert($('status').textContent.includes('已匯入'), '匯入應成功: ' + $('status').textContent);
  console.log('匯入範例模型 OK');

  // ---- 選「反應譜」, 預設是規範來源 ----
  const sel = $('analysisType');
  sel.value = 'rsa'; sel.dispatchEvent(new w.Event('change'));
  assert($('rsaConfigBar').style.display === 'block', '選反應譜後應顯示設定列');
  assert($('cyclicConfigBar').style.display === 'none' && $('modalConfigBar').style.display === 'none', '其他設定列應隱藏');
  assert($('rsa_spec_type').value === 'code' && $('rsa_code_panel').style.display !== 'none', '預設應該是規範來源');
  assert($('rsa_custom_panel').style.display === 'none', '自訂面板預設隱藏');

  // ---- 輸入驗證 ----
  $('rsa_sds').value = '-1';
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('S_DS 必須是正數'), '負SDS應被擋下: ' + $('status').textContent);
  $('rsa_sds').value = '0.6'; $('rsa_sd1').value = '0.35'; $('rsa_tl').value = '6.0';
  $('rsa_damping').value = '1.5';
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('阻尼比'), '阻尼比超出範圍應被擋下: ' + $('status').textContent);
  $('rsa_damping').value = '0.05';
  assert(ev(`rsaResult`) === null, '驗證失敗時不應有結果');
  console.log('前端輸入驗證(負SDS、阻尼比超出範圍) OK');

  // ---- 正式求解(規範反應譜) ----
  $('rsa_nmodes').value = '4'; $('rsa_kind').value = 'consistent'; $('rsa_combine').value = 'CQC';
  $('btnSolve').click();
  await waitFor(() => ev(`rsaResult !== null`), '反應譜分析結果');
  const r = ev(`JSON.parse(JSON.stringify(rsaResult))`);
  assert(r.n_modes === expected.nModes, `模態數 ${r.n_modes} 應等於核心的 ${expected.nModes}`);
  assert(Math.abs(r.base_shear / expected.baseShear - 1) < 1e-9, `基底剪力 ${r.base_shear} 應等於核心的 ${expected.baseShear}`);
  assert(Math.abs(r.cum_ratio_total - expected.cumRatioTotal) < 1e-9, `累積質量比 ${r.cum_ratio_total} 應等於核心的 ${expected.cumRatioTotal}`);
  for (let i = 0; i < expected.periods.length; i++) {
    assert(Math.abs(r.modes[i].period / expected.periods[i] - 1) < 1e-9, `第${i + 1}模態週期不符`);
    assert(Math.abs(r.modes[i].gamma - expected.gammas[i]) < 1e-9 * Math.max(1, Math.abs(expected.gammas[i])),
           `第${i + 1}模態 Γ ${r.modes[i].gamma} 應等於核心的 ${expected.gammas[i]}`);
  }
  console.log(`結果與核心逐項一致(含Γ): base_shear=${r.base_shear.toFixed(2)} N, cum_ratio_total=${r.cum_ratio_total.toFixed(4)}, `
    + `週期 ${r.modes.map(m => m.period.toFixed(5)).join(', ')}`);

  // ---- 反應譜曲線視圖 ----
  assert(ev(`viewMode`) === 'rsa', 'Solve後應自動切到反應譜視圖');
  assert($('btnViewRsa').style.display !== 'none' && $('btnViewRsa').classList.contains('active'), '反應譜按鈕應顯示且啟用');
  assert($('rsaPanel').style.display === 'block', '結果面板應顯示');
  assert(d.querySelectorAll('#overlay polyline.rsaCurveLine').length === 1, '應有 1 條反應譜曲線');
  const curvePts = d.querySelector('#overlay polyline.rsaCurveLine').getAttribute('points').trim().split(/\s+/);
  assert(curvePts.length === 300, `曲線取樣點數應為300, 實際 ${curvePts.length}`);
  assert(!curvePts.some(p => p.includes('NaN')), '曲線座標不應有 NaN');
  assert($('overlay').textContent.includes('反應譜曲線'), '應顯示標題');
  console.log('反應譜曲線圖: 300 個取樣點、無 NaN、標題正確 OK');

  // ---- 結果表(渲染出來的文字, 不是底層資料, 這樣才抓得到只發生在渲染層的錯誤) ----
  let mrows = d.querySelectorAll('#rsaModalTable tbody tr');
  assert(mrows.length === 4, `模態表應有 4 列, 實際 ${mrows.length}`);
  for (let i = 0; i < expected.gammas.length; i++) {
    const cellText = mrows[i].children[4].textContent;   // 欄位順序: 模態,T,f,S_a,Γ,...
    assert(Math.abs(parseFloat(cellText) - expected.gammas[i]) < 1e-3 * Math.max(1, Math.abs(expected.gammas[i])),
           `模態表第${i + 1}列的Γ欄應顯示 ${expected.gammas[i].toFixed(4)}, 實際 ${cellText}`);
  }
  console.log('模態表的 Γ 欄位(渲染文字)跟核心逐項一致 OK');
  let nrows = d.querySelectorAll('#rsaNodeTable tbody tr');
  assert(nrows.length === 4, `節點位移表應有 4 列, 實際 ${nrows.length}`);
  let mbrows = d.querySelectorAll('#rsaMemberTable tbody tr');
  assert(mbrows.length === 3, `桿件內力表應有 3 列, 實際 ${mbrows.length}`);
  assert($('rsa_summary').textContent.includes('CQC') && $('rsa_summary').textContent.includes('一致質量'),
         '摘要應顯示組合方法與質量矩陣: ' + $('rsa_summary').textContent);
  assert($('status').textContent.includes('反應譜分析完成'), '狀態列: ' + $('status').textContent);
  console.log('模態表 4 列、節點位移表 4 列、桿件內力表 3 列、摘要正確 OK');

  // ---- 切換單位: 表頭與數值跟著換 ----
  const fmax0 = parseFloat(mrows[0].children[7].textContent);   // 模態基底剪力欄
  $('unitForce').value = 'N'; $('unitForce').dispatchEvent(new w.Event('change'));
  assert(d.querySelector('#rsaModalTable thead').textContent.includes('(N)'), '表頭單位應改成 N');
  mrows = d.querySelectorAll('#rsaModalTable tbody tr');
  const fmax1 = parseFloat(mrows[0].children[7].textContent);
  assert(Math.abs(fmax1 / fmax0 - 1000) < 1e-2, `力單位切換: 模態基底剪力應放大1000倍: ${fmax0} kN -> ${fmax1} N`);
  $('unitForce').value = 'kN'; $('unitForce').dispatchEvent(new w.Event('change'));
  console.log(`單位切換: 模態基底剪力 ${fmax0} kN -> ${fmax1} N OK`);

  // ---- 切到自訂反應譜 ----
  $('rsa_spec_type').value = 'custom'; $('rsa_spec_type').dispatchEvent(new w.Event('change'));
  assert($('rsa_custom_panel').style.display !== 'none' && $('rsa_code_panel').style.display === 'none', '應切到自訂面板');
  let prows = d.querySelectorAll('#rsa_points_table tbody tr');
  assert(prows.length === 2, `自訂反應譜預設應有 2 個點, 實際 ${prows.length}`);
  $('rsa_pt_add').click(); await sleep(20);
  prows = d.querySelectorAll('#rsa_points_table tbody tr');
  assert(prows.length === 3, '新增一點後應有 3 列');
  prows[2].querySelector('input').value = '1.0';
  prows[2].querySelector('input').dispatchEvent(new w.Event('input'));
  assert(ev(`rsaCustomPoints[2][0]`) === 1.0, '編輯後的週期應該同步到 rsaCustomPoints');
  prows[0].querySelectorAll('button')[0].click(); await sleep(20);
  prows = d.querySelectorAll('#rsa_points_table tbody tr');
  assert(prows.length === 2, '刪除一列後應剩 2 列');
  assert(prows[0].querySelector('button').disabled, '只剩 2 點時刪除按鈕應該停用(至少要2點)');
  console.log('自訂反應譜: 新增/編輯/刪除資料點 OK(下限2點時刪除鈕停用)');

  // ---- 只有1個點時, buildRsaPayload() 必須在送出前就擋下(不依賴後端才擋) ----
  const savedPoints = ev(`JSON.stringify(rsaCustomPoints)`);
  ev(`rsaCustomPoints.length = 1;`);
  assert(ev(`buildRsaPayload()`) === null, '只有1個自訂點時 buildRsaPayload() 應該回傳 null(前端就擋下, 不用等後端回應)');
  assert($('status').textContent.includes('至少需要'), '應顯示至少需要2點的訊息: ' + $('status').textContent);
  ev(`rsaCustomPoints = ${savedPoints};`);
  console.log('只有1個自訂點: 前端在送出前就擋下 OK');

  $('btnSolve').click();
  await waitFor(() => $('status').textContent.includes('反應譜分析完成') && ev(`rsaResult.n_modes`) === 4
                && d.querySelector('#rsaModalTable tbody tr td'), '自訂反應譜結果');
  await sleep(100);
  assert(ev(`viewMode`) === 'rsa', '自訂反應譜求解後仍在反應譜視圖');
  console.log('自訂反應譜求解成功 OK');

  // ---- 後端錯誤: 拿掉支承 -> 清楚訊息, 舊結果不被破壞 ----
  const oldNModes = ev(`rsaResult.n_modes`);
  const supBackup = ev(`JSON.stringify(model.supports)`);
  ev(`model.supports = [];`);
  $('rsa_spec_type').value = 'code'; $('rsa_spec_type').dispatchEvent(new w.Event('change'));
  $('btnSolve').click();
  await waitFor(() => $('status').textContent.includes('反應譜分析失敗'), '錯誤訊息');
  assert($('status').textContent.includes('機構') || $('status').textContent.includes('支承'),
         '錯誤訊息應說明機構/支承不足: ' + $('status').textContent);
  assert(ev(`rsaResult.n_modes`) === oldNModes, '失敗時不應清掉上一次的結果');
  console.log('機構/支承不足: 顯示清楚訊息 OK ->', $('status').textContent.slice(0, 50));
  ev(`model.supports = ${supBackup};`);

  // ---- 新建/載入重置 ----
  ev(`resetViewToStructure();`);
  assert(ev(`rsaResult`) === null && $('btnViewRsa').style.display === 'none' && $('rsaConfigBar').style.display === 'none'
         && $('rsaPanel').style.display === 'none', '重置應清掉結果、按鈕、設定列與面板');
  assert(ev(`viewMode`) === 'structure' && $('analysisType').value === 'linear', '回到結構視圖與線性分析');
  console.log('新建/載入重置 OK');

  assert(errs.length === 0, '頁面有 JS 錯誤: ' + errs.slice(0, 3).join(' | '));
  console.log('\n全部通過: 網頁反應譜分析端到端(匯入模型、設定、Solve、曲線圖、結果表、單位切換、自訂反應譜、錯誤訊息、重置)');
  process.exit(0);
})().catch(e => { console.error('FAILED:', e.message); process.exit(1); });
