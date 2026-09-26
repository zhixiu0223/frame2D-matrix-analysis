/*
 * 網頁端到端測試(選用): 反覆載重(遲滯)分析 (D7b)
 *
 * 用 jsdom 載入真正由後端提供的 index.html, window.fetch 接到真的 Python 後端, 然後「像使用者一樣」:
 * 走「匯入JSON」載入範例模型(examples/portal_cyclic_demo.json)→ 選「循環」分析 → 填控制節點/幅值/
 * 圈數/步長 → Solve → 檢查遲滯迴圈圖、塑鉸 M-θp 圖、結果表、單位切換、錯誤訊息、重置。
 * 網頁拿到的耗能/最大底剪力跟核心(frame2d.cyclic)直接算的逐項比對。
 *
 * 用法(由 tests/test_web_cyclic_e2e.py 呼叫):
 *   node tests/web/cyclic_e2e.js <後端網址> <範例JSON路徑> <預期結果JSON>
 * 需要 jsdom:  npm install jsdom
 */
const { JSDOM, VirtualConsole } = require('jsdom');
const fs = require('fs');

const base = process.argv[2];
const modelFile = process.argv[3];
const expected = JSON.parse(process.argv[4]);   // {loops:[energyJ...], hingeLabels:[...], maxF, nPoints, nYield}

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

  // ---- 走「匯入JSON」載入範例模型 ----
  const input = $('importJsonFile');
  const file = new w.File([fs.readFileSync(modelFile, 'utf8')], 'portal_cyclic_demo.json', {type: 'application/json'});
  Object.defineProperty(input, 'files', {value: [file], configurable: true});
  input.dispatchEvent(new w.Event('change'));
  await sleep(300);
  assert($('status').textContent.includes('已匯入'), '匯入應成功: ' + $('status').textContent);
  assert(ev(`model.members[0].Mp_i`) === 100000, '範例模型的塑鉸容量應被載入');
  console.log('匯入範例模型 OK');

  // ---- 選「循環」分析 ----
  const sel = $('analysisType');
  sel.value = 'cyclic'; sel.dispatchEvent(new w.Event('change'));
  assert($('cyclicConfigBar').style.display === 'block', '選循環後應顯示設定列');
  assert($('pushoverConfigBar').style.display === 'none' && $('modalConfigBar').style.display === 'none', '其他設定列應隱藏');
  assert($('cy_amp_unit').textContent === 'mm' && $('cy_step_unit').textContent === 'mm', '幅值/步長單位標籤應跟位移單位一致(預設 mm)');

  // ---- 輸入驗證(不送出請求) ----
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('需要填'), '欄位全空應提示: ' + $('status').textContent);
  $('cy_nodes').value = '1'; $('cy_amps').value = '15,-30'; $('cy_step').value = '2';
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('正數'), '負幅值應被前端擋下: ' + $('status').textContent);
  $('cy_amps').value = '15,30,45,60'; $('cy_step').value = '0.00001';
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('預估要走'), '步長極小應被前端擋下(預估步數): ' + $('status').textContent);
  $('cy_nodes').value = '1,x';
  $('cy_step').value = '2'; $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('控制節點格式錯誤'), '節點格式錯誤: ' + $('status').textContent);
  assert(ev(`cyclicResult`) === null, '驗證失敗時不應有結果');
  console.log('前端輸入驗證(全空、負幅值、步數過多、節點格式) OK');

  // ---- 正式求解 ----
  $('cy_nodes').value = '1'; $('cy_ncycles').value = '2';
  $('btnSolve').click();
  await waitFor(() => ev(`cyclicResult !== null`), '反覆載重結果');
  const r = ev(`JSON.parse(JSON.stringify(cyclicResult))`);
  assert(r.u.length === expected.nPoints, `記錄點數 ${r.u.length} 應等於核心的 ${expected.nPoints}`);
  const maxF = Math.max(...r.F.map(Math.abs));
  assert(Math.abs(maxF / expected.maxF - 1) < 1e-9, `最大底剪力 ${maxF} 應等於核心的 ${expected.maxF}`);
  for (let i = 0; i < expected.loops.length; i++) {
    const got = r.loops[i].energy, want = expected.loops[i];
    assert(Math.abs(got - want) <= 1e-9 * Math.max(Math.abs(want), 1), `第${i + 1}級迴圈耗能 ${got} 應等於核心的 ${want}`);
  }
  assert(JSON.stringify(r.hinges.map(h => h.label)) === JSON.stringify(expected.hingeLabels), '塑鉸降伏順序應與核心一致: ' + r.hinges.map(h => h.label));
  console.log(`結果與核心逐項一致: ${r.u.length} 點, 耗能 ${r.loops.map(l => l.energy.toFixed(1)).join(', ')} J, 最大底剪力 ${maxF.toFixed(1)} N`);

  // ---- 遲滯迴圈圖 ----
  assert(ev(`viewMode`) === 'cyclic_loop', 'Solve後應自動切到遲滯迴圈視圖');
  assert($('btnViewCyclic').style.display !== 'none' && $('btnViewCyclicHinge').style.display !== 'none', '兩個結果按鈕應顯示');
  assert($('btnViewCyclic').classList.contains('active'), '遲滯迴圈按鈕應啟用');
  assert($('cyclicPanel').style.display === 'block', '結果面板應顯示');
  assert(d.querySelectorAll('#overlay polyline.cyclicLoopLine').length === 1, '整體迴圈應有 1 條折線');
  assert(d.querySelectorAll('#overlay circle.cyclicYieldDot').length === expected.nYield, `降伏事件紅點應有 ${expected.nYield} 個`);
  assert($('overlay').textContent.includes('整體結構遲滯迴圈'), '圖標題');
  const pts = d.querySelector('#overlay polyline.cyclicLoopLine').getAttribute('points').trim().split(/\s+/);
  assert(pts.length === r.u.length, `折線點數應等於記錄點數(${pts.length} vs ${r.u.length})`);
  assert(!pts.some(p => p.includes('NaN')), '折線座標不應有 NaN');
  let rows = d.querySelectorAll('#cyclicLoopTable tbody tr');
  assert(rows.length === 4, `迴圈耗能表應有 4 列, 實際 ${rows.length}`);
  assert(rows[0].children[2].textContent === '0', '彈性幅值的耗能應顯示 0 而不是 e-雜訊: ' + rows[0].children[2].textContent);
  assert(parseFloat(rows[3].children[3].textContent) > 0.25, '最大幅值的 ξ_eq 應該 > 0.25: ' + rows[3].children[3].textContent);
  // 彈性幅值的耗能若是浮點雜訊(例如 3e-13 J)也要顯示成 0, 不顯示科學記號雜訊
  const keepE = ev(`cyclicResult.loops[0].energy`), keepX = ev(`cyclicResult.loops[0].xi_eq`);
  ev(`cyclicResult.loops[0].energy = 3.2e-13; cyclicResult.loops[0].xi_eq = 1.1e-16; renderCyclicPanel();`);
  rows = d.querySelectorAll('#cyclicLoopTable tbody tr');
  assert(rows[0].children[2].textContent === '0' && rows[0].children[3].textContent === '0', '浮點雜訊的耗能/ξ_eq 應顯示 0: ' + rows[0].children[2].textContent + ' / ' + rows[0].children[3].textContent);
  ev(`cyclicResult.loops[0].energy = ${keepE}; cyclicResult.loops[0].xi_eq = ${keepX}; renderCyclicPanel();`);
  const hrows = d.querySelectorAll('#cyclicHingeTable tbody tr');
  assert(hrows.length === r.hinges.length && hrows[0].children[1].textContent === expected.hingeLabels[0], '降伏順序表');
  assert($('status').textContent.includes('反覆載重完成'), '狀態列: ' + $('status').textContent);
  console.log('遲滯迴圈圖(折線、降伏紅點、標題)、耗能表、降伏順序表 OK');

  // ---- 塑鉸 M-θp 圖 ----
  $('btnViewCyclicHinge').click(); await sleep(50);
  assert(ev(`viewMode`) === 'cyclic_hinge', '切到塑鉸視圖');
  assert(d.querySelectorAll('#overlay polyline.cyclicHingeLine').length === r.hinges.length, `塑鉸小圖應有 ${r.hinges.length} 張`);
  assert($('overlay').textContent.includes('1. ' + expected.hingeLabels[0]), '第一張小圖標題含降伏順序與塑鉸名稱');
  assert(d.querySelectorAll('#overlay polyline.cyclicLoopLine').length === 0, '塑鉸視圖不應再畫整體迴圈');
  assert($('cyclicPanel').style.display === 'block', '塑鉸視圖也顯示結果面板');
  console.log('塑鉸 M-θp 小圖網格 OK');

  // ---- 切換單位: 表頭與數值跟著換 ----
  const amp0 = parseFloat(d.querySelectorAll('#cyclicLoopTable tbody tr')[0].children[0].textContent);
  const fmax0 = parseFloat(d.querySelectorAll('#cyclicLoopTable tbody tr')[0].children[1].textContent);
  assert(Math.abs(amp0 - 15) < 1e-9, '預設位移單位 mm: 第一級幅值 15: ' + amp0);
  $('unitDisp').value = 'm'; $('unitDisp').dispatchEvent(new w.Event('change'));
  $('unitForce').value = 'N'; $('unitForce').dispatchEvent(new w.Event('change'));
  assert($('cy_amp_unit').textContent === 'm', '設定列的單位標籤應跟著換');
  assert(d.querySelector('#cyclicLoopTable thead').textContent.includes('幅值(m)') && d.querySelector('#cyclicLoopTable thead').textContent.includes('F_max(N)'), '表頭單位');
  const amp1 = parseFloat(d.querySelectorAll('#cyclicLoopTable tbody tr')[0].children[0].textContent);
  const fmax1 = parseFloat(d.querySelectorAll('#cyclicLoopTable tbody tr')[0].children[1].textContent);
  assert(Math.abs(amp1 - 0.015) < 1e-12 && Math.abs(fmax1 / fmax0 - 1000) < 1e-3, `單位換算: 幅值 ${amp0}→${amp1}, F_max ${fmax0}→${fmax1}`);
  assert($('overlay').textContent.includes('(N)') || d.querySelectorAll('#overlay text').length > 0, '圖的軸標題應重畫');
  $('unitDisp').value = 'mm'; $('unitDisp').dispatchEvent(new w.Event('change'));
  $('unitForce').value = 'kN'; $('unitForce').dispatchEvent(new w.Event('change'));
  console.log(`單位切換: 幅值 ${amp0} mm → ${amp1} m, F_max ${fmax0} kN → ${fmax1} N OK`);

  // ---- 後端錯誤: 拿掉塑鉸容量 -> 清楚訊息, 舊結果不被破壞 ----
  const backup = ev(`JSON.stringify(model.members)`);
  ev(`model.members = model.members.map(m => ({...m, Mp_i: null, Mp_j: null}));`);
  $('btnViewCyclic').click(); await sleep(20);
  $('btnSolve').click();
  await waitFor(() => $('status').textContent.includes('反覆載重分析失敗'), '錯誤訊息');
  assert($('status').textContent.includes('塑鉸容量'), '應說明沒有塑鉸容量: ' + $('status').textContent);
  assert(ev(`cyclicResult.u.length`) === r.u.length, '失敗時不應清掉上一次的結果');
  ev(`model.members = ${backup};`);
  console.log('沒有塑鉸容量: 顯示清楚訊息 OK ->', $('status').textContent.slice(0, 50));

  // ---- 匯出 Markdown ----
  ev(`window.__capturedBlob = null; window.__capturedName = null;
      window.downloadBlob = (blob, name) => { window.__capturedBlob = blob; window.__capturedName = name; };`);
  $('btnExportMd').click();
  await sleep(50);
  const exportName = ev(`window.__capturedName`);
  assert(exportName && exportName.endsWith('.md'), '應該產生一個 .md 檔案: ' + exportName);
  const mdText = await ev(`window.__capturedBlob.text()`);
  assert(mdText.includes('# frame2d 循環(遲滯)分析結果'), '應該是循環分析報告的標題');
  assert(mdText.includes('## 輸入資料'), '應該包含輸入資料段落(節點/斷面/桿件/塑鉸容量, 讓人可以重建模型)');
  assert(mdText.includes('## 分析設定'), '應該包含控制節點/方向/位移幅值/圈數等分析設定');
  const wantAmps = ev(`cyclicResult.amplitudes.length`);
  assert(mdText.includes('位移幅值'), '應該記錄位移幅值');
  const wantNPts = ev(`cyclicResult.u.length`);
  const rowsMatch = mdText.match(/### 完整力-位移曲線[\s\S]*/);
  assert(rowsMatch, '應該有完整力-位移曲線段落');
  const nRows = (rowsMatch[0].match(/^\| \d+ \|/gm) || []).length;
  assert(nRows === wantNPts, `完整曲線應該有 ${wantNPts} 列(逐點對照畫面上的迴圈), 實際 ${nRows}`);
  console.log(`匯出 Markdown: 標題/輸入資料/分析設定/完整力-位移曲線(${nRows}點)都正確 OK`);

  // ---- 匯出 PDF ----
  ev(`window.__capturedBlob = null; window.__capturedName = null;`);
  $('btnExportPdf').click();
  await waitFor(() => ev(`window.__capturedBlob !== null`), 'PDF匯出完成', 15000);
  const pdfName = ev(`window.__capturedName`);
  assert(pdfName && pdfName.endsWith('_cyclic.pdf'), '檔名應該以 _cyclic.pdf 結尾: ' + pdfName);
  const pdfBlob = ev(`window.__capturedBlob`);
  assert(pdfBlob.type === 'application/pdf', 'blob類型應該是application/pdf: ' + pdfBlob.type);
  const pdfBuf = await ev(`window.__capturedBlob.arrayBuffer()`);
  const pdfHead = Buffer.from(pdfBuf.slice(0, 4)).toString('ascii');
  assert(pdfHead === '%PDF', 'PDF檔案應該以%PDF開頭: ' + pdfHead);
  assert(pdfBuf.byteLength > 10000, `PDF檔案大小應該有相當內容, 實際只有${pdfBuf.byteLength}bytes`);
  console.log(`匯出 PDF: application/pdf、檔名_cyclic.pdf、有效PDF檔頭、大小${pdfBuf.byteLength}bytes OK`);

  // ---- 新建/載入重置 ----
  ev(`resetViewToStructure();`);
  assert(ev(`cyclicResult`) === null && $('btnViewCyclic').style.display === 'none' && $('btnViewCyclicHinge').style.display === 'none'
         && $('cyclicConfigBar').style.display === 'none' && $('cyclicPanel').style.display === 'none', '重置應清掉結果、按鈕、設定列與面板');
  assert(ev(`viewMode`) === 'structure' && $('analysisType').value === 'linear', '回到結構視圖與線性分析');
  console.log('新建/載入重置 OK');

  assert(errs.length === 0, '頁面有 JS 錯誤: ' + errs.slice(0, 3).join(' | '));
  console.log('\n全部通過: 網頁遲滯分析端到端(匯入模型、設定、Solve、迴圈圖、M-θp 圖、結果表、單位切換、錯誤訊息、重置)');
  process.exit(0);
})().catch(e => { console.error('FAILED:', e.message); process.exit(1); });
