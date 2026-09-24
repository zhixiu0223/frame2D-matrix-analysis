/*
 * 網頁端到端測試(選用): 非線性地震反應分析 (D6b)
 *
 * 用 jsdom 載入真正由後端提供的 index.html, window.fetch 接到真的 Python 後端, 然後「像使用者
 * 一樣」: 匯入範例模型 → 選「非線性地震」→ 衰減脈衝地震歷程 Solve → 檢查地震動畫視圖、播放/
 * 暫停/滑桿、時程圖、遲滯迴圈圖、塑鉸M-θp圖、結果表、單位切換、自訂地震歷程、錯誤訊息、重置。
 * 網頁拿到的尖峰位移/降伏次數/能量平衡跟核心 nonlinear_seismic_web_analysis() 直接算的逐項
 * 比對。
 *
 * 用法(由 tests/test_web_seismic_e2e.py 呼叫):
 *   node tests/web/seismic_e2e.js <後端網址> <範例JSON路徑> <預期結果JSON>
 * 需要 jsdom:  npm install jsdom
 */
const { JSDOM, VirtualConsole } = require('jsdom');
const fs = require('fs');

const base = process.argv[2];
const modelFile = process.argv[3];
const expected = JSON.parse(process.argv[4]);   // {nSteps, peakDisp, nFrames, hingeLabels, energyFinal}

function assert(cond, msg) { if (!cond) { throw new Error('ASSERT: ' + msg); } }
const sleep = ms => new Promise(r => setTimeout(r, ms));
async function waitFor(fn, what, ms = 20000) {
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
  const file = new w.File([fs.readFileSync(modelFile, 'utf8')], 'portal_seismic_demo.json', {type: 'application/json'});
  Object.defineProperty(input, 'files', {value: [file], configurable: true});
  input.dispatchEvent(new w.Event('change'));
  await sleep(300);
  assert($('status').textContent.includes('已匯入'), '匯入應成功: ' + $('status').textContent);
  console.log('匯入範例模型 OK');

  // ---- 選「非線性地震」, 預設是脈衝來源 ----
  const sel = $('analysisType');
  sel.value = 'seismic'; sel.dispatchEvent(new w.Event('change'));
  assert($('seismicConfigBar').style.display === 'block', '選非線性地震後應顯示設定列');
  assert($('cyclicConfigBar').style.display === 'none' && $('rsaConfigBar').style.display === 'none', '其他設定列應隱藏');
  assert($('sm_gm_type').value === 'pulse' && $('sm_pulse_panel').style.display !== 'none', '預設應該是脈衝來源');
  assert($('sm_custom_panel').style.display === 'none', '自訂面板預設隱藏');

  // ---- 輸入驗證 ----
  $('sm_node').value = '';
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('請填控制節點'), '缺控制節點應被擋下: ' + $('status').textContent);
  $('sm_node').value = '1'; $('sm_modes').value = 'a,b';
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('Rayleigh控制模態格式錯誤'), '模態格式錯誤應被擋下: ' + $('status').textContent);
  $('sm_modes').value = '1,3';
  assert(ev(`seismicResult`) === null, '驗證失敗時不應有結果');
  console.log('前端輸入驗證(缺控制節點、模態格式錯誤) OK');

  // ---- 正式求解(脈衝地震歷程) ----
  $('sm_zeta').value = '0.05'; $('sm_kind').value = 'lumped';
  $('sm_pulse_amp').value = '0.5'; $('sm_pulse_freq').value = '1.5'; $('sm_pulse_decay').value = '0.25';
  $('btnSolve').click();
  await waitFor(() => ev(`seismicResult !== null`), '非線性地震反應結果', 30000);
  const r = ev(`JSON.parse(JSON.stringify(seismicResult))`);
  assert(r.n_steps === expected.nSteps, `步數 ${r.n_steps} 應等於核心的 ${expected.nSteps}`);
  assert(Math.abs(r.peak_displacement / expected.peakDisp - 1) < 1e-9, `尖峰位移 ${r.peak_displacement} 應等於核心的 ${expected.peakDisp}`);
  assert(r.frames.length === expected.nFrames, `動畫影格數 ${r.frames.length} 應等於核心的 ${expected.nFrames}`);
  const energyFinal = r.energy.KE[r.energy.KE.length - 1] + r.energy.Wdamp[r.energy.Wdamp.length - 1] + r.energy.Wint[r.energy.Wint.length - 1];
  assert(Math.abs(energyFinal / expected.energyFinal - 1) < 1e-6, `能量平衡最終值 ${energyFinal} 應等於核心的 ${expected.energyFinal}`);
  console.log(`結果與核心逐項一致: n_steps=${r.n_steps}, peak_disp=${(r.peak_displacement * 1000).toFixed(2)}mm, `
    + `動畫影格=${r.frames.length}, 塑鉸=${r.hinges.map(h => h.label).join(',')}`);

  // ---- 地震動畫視圖 ----
  assert(ev(`viewMode`) === 'seismic_anim', 'Solve後應自動切到地震動畫視圖');
  assert($('btnViewSeisAnim').style.display !== 'none' && $('btnViewSeisAnim').classList.contains('active'), '動畫按鈕應顯示且啟用');
  assert($('seisPlaybackBar').style.display !== 'none', '播放列應顯示');
  assert($('sm_slider').max === String(r.frames.length - 1), '滑桿max應等於影格數-1');
  assert(d.querySelectorAll('#content line.deformed-line').length === ev(`model.members.length`), '應該畫出每根桿件的變形直線');
  // 獨立重算「變形放大倍率」該有的值(對照 computeSeismicDrawScale() 的公式), 這樣「忘記除以
  // maxDisp」這類讓倍率大一個數量級的錯誤才抓得到(只看畫面上有沒有線、線的端點座標有沒有
  // NaN, 對這種「倍率算得離譜但仍是個有效數字」的錯誤是抓不到的)
  let maxDispCheck = 1e-9;
  for (const fr of r.frames) for (const nid in fr.nodes) maxDispCheck = Math.max(maxDispCheck, Math.hypot(fr.nodes[nid].ux, fr.nodes[nid].uy));
  const xsCheck = ev(`model.nodes.map(n => n.x)`), ysCheck = ev(`model.nodes.map(n => n.y)`);
  const diagCheck = Math.hypot(Math.max(...xsCheck) - Math.min(...xsCheck), Math.max(...ysCheck) - Math.min(...ysCheck)) || 1;
  const expectedScale = (0.15 * diagCheck) / maxDispCheck;
  const gotScale = ev(`seismicDrawScale`);
  assert(Math.abs(gotScale / expectedScale - 1) < 1e-6, `變形放大倍率 ${gotScale} 應該等於獨立算出的 ${expectedScale}`);
  console.log(`  變形放大倍率 = ${gotScale.toFixed(3)}(獨立算出的期望值 ${expectedScale.toFixed(3)}) OK`);
  console.log('地震動畫視圖(播放列、變形直線) OK');

  // ---- 播放/暫停/滑桿 ----
  assert($('sm_play').textContent.includes('播放'), '初始應顯示「播放」');
  $('sm_play').click(); await sleep(200);
  assert($('sm_play').textContent.includes('暫停'), '按下後應顯示「暫停」: ' + $('sm_play').textContent);
  const idxAfterPlay = ev(`seismicFrameIndex`);
  assert(idxAfterPlay > 0, '播放一段時間後影格索引應該前進: ' + idxAfterPlay);
  $('sm_play').click(); await sleep(50);
  assert($('sm_play').textContent.includes('播放'), '再按一次應該回到「播放」(暫停)');
  const idxAfterPause = ev(`seismicFrameIndex`);
  await sleep(150);
  assert(ev(`seismicFrameIndex`) === idxAfterPause, '暫停後影格索引不應再變動');
  $('sm_slider').value = '5'; $('sm_slider').dispatchEvent(new w.Event('input'));
  assert(ev(`seismicFrameIndex`) === 5, '拖動滑桿應該直接跳到該影格');
  assert($('sm_frame_label').textContent.includes('5/'), '影格標籤應該更新: ' + $('sm_frame_label').textContent);
  console.log('播放/暫停/滑桿 OK');

  // ---- 時程圖 ----
  const btnHist = d.querySelector('.viewBtn[data-view="seismic_history"]');
  btnHist.click(); await sleep(50);
  assert(ev(`viewMode`) === 'seismic_history', '應切到地震時程視圖');
  assert(d.querySelectorAll('#overlay polyline.seisGroundLine').length === 1, '應有1條地面加速度曲線');
  assert(d.querySelectorAll('#overlay polyline.seisDispLine').length === 1, '應有1條位移曲線');
  assert(d.querySelectorAll('#overlay circle.seisNowDot, #overlay circle[r]').length >= 2, '應該有目前時間點的標記');
  console.log('地震時程視圖(地面加速度+位移雙圖) OK');

  // ---- 遲滯迴圈圖 ----
  const btnHyst = d.querySelector('.viewBtn[data-view="seismic_hyst"]');
  btnHyst.click(); await sleep(50);
  assert(ev(`viewMode`) === 'seismic_hyst', '應切到遲滯迴圈視圖');
  assert(d.querySelectorAll('#overlay polyline.seisHystLine').length === 1, '應有1條遲滯迴圈曲線');
  const hystPts = d.querySelector('#overlay polyline.seisHystLine').getAttribute('points').trim().split(/\s+/);
  assert(hystPts.length === r.control_disp.length, `遲滯迴圈點數應等於完整時間序列長度(${r.control_disp.length}), 實際 ${hystPts.length}`);
  assert(!hystPts.some(p => p.includes('NaN')), '遲滯迴圈座標不應有 NaN');
  console.log(`遲滯迴圈視圖(${hystPts.length} 點, 用完整時間序列不是抽稀後的動畫幀) OK`);

  // ---- 塑鉸 M-θp 圖 ----
  const btnHinge = d.querySelector('.viewBtn[data-view="seismic_hinge"]');
  btnHinge.click(); await sleep(50);
  assert(ev(`viewMode`) === 'seismic_hinge', '應切到塑鉸視圖');
  assert(d.querySelectorAll('#overlay polyline.seisHystLine').length === r.hinges.length, `應有 ${r.hinges.length} 張塑鉸小圖`);
  console.log('塑鉸 M-θp 小圖網格 OK');

  // ---- 結果表 ----
  assert($('sm_summary').textContent.includes('控制節點1'), '摘要應顯示控制節點: ' + $('sm_summary').textContent);
  const hrows = d.querySelectorAll('#seismicHingeTable tbody tr');
  assert(hrows.length === r.hinges.length, `塑鉸表應有 ${r.hinges.length} 列`);
  console.log('結果面板摘要與塑鉸表 OK');

  // ---- 單位切換 ----
  const peak0 = r.peak_displacement;
  $('unitDisp').value = 'm'; $('unitDisp').dispatchEvent(new w.Event('change'));
  assert($('sm_summary').textContent.includes('尖峰位移'), '摘要應該重新渲染');
  $('unitDisp').value = 'mm'; $('unitDisp').dispatchEvent(new w.Event('change'));
  console.log('單位切換不報錯 OK');

  // ---- 切到自訂地震歷程 ----
  $('sm_gm_type').value = 'custom'; $('sm_gm_type').dispatchEvent(new w.Event('change'));
  assert($('sm_custom_panel').style.display !== 'none' && $('sm_pulse_panel').style.display === 'none', '應切到自訂面板');
  let prows = d.querySelectorAll('#sm_points_table tbody tr');
  assert(prows.length === 2, `自訂地震歷程預設應有 2 個點, 實際 ${prows.length}`);
  $('sm_pt_add').click(); await sleep(20);
  prows = d.querySelectorAll('#sm_points_table tbody tr');
  assert(prows.length === 3, '新增一點後應有 3 列');
  prows[0].querySelectorAll('button')[0].click(); await sleep(20);
  prows = d.querySelectorAll('#sm_points_table tbody tr');
  assert(prows.length === 2, '刪除一列後應剩 2 列');
  assert(prows[0].querySelector('button').disabled, '只剩 2 點時刪除按鈕應該停用');
  console.log('自訂地震歷程: 新增/刪除資料點 OK(下限2點時刪除鈕停用)');

  // ---- PEER NGA 檔案上傳 ----
  $('sm_gm_type').value = 'peer_nga'; $('sm_gm_type').dispatchEvent(new w.Event('change'));
  assert($('sm_peer_nga_panel').style.display !== 'none' && $('sm_custom_panel').style.display === 'none', '應切到PEER NGA面板');
  assert($('sm_peer_nga_status').textContent.includes('尚未選擇'), '初始應顯示尚未選擇檔案');
  $('btnSolve').click(); await sleep(50);
  assert($('status').textContent.includes('請先選擇一個'), '沒選檔案時應該被前端擋下: ' + $('status').textContent);
  const dtGm = 0.02, nGm = 150;
  const linesGm = ['SYNTHETIC EQ FOR E2E TEST', 'STATION Y', 'ACCELERATION TIME HISTORY IN UNITS OF G', `NPTS=${nGm}, DT=${dtGm} SEC`];
  for (let i = 0; i < nGm; i += 5) {
    const row = [];
    for (let k = i; k < Math.min(i + 5, nGm); k++) {
      const t = k * dtGm;
      row.push((0.3 * Math.sin(2 * Math.PI * 1.0 * t) * Math.exp(-0.1 * t)).toFixed(6));
    }
    linesGm.push(row.join(' '));
  }
  const at2Text = linesGm.join('\n') + '\n';
  const at2File = new w.File([at2Text], 'synthetic.AT2', {type: 'text/plain'});
  Object.defineProperty($('sm_peer_nga_file'), 'files', {value: [at2File], configurable: true});
  $('sm_peer_nga_file').dispatchEvent(new w.Event('change'));
  await sleep(100);
  assert($('sm_peer_nga_status').textContent.includes('synthetic.AT2'), '應該顯示已選擇的檔名: ' + $('sm_peer_nga_status').textContent);
  assert(ev(`seismicPeerNgaText`) === at2Text, '讀到的檔案內容應該跟原始文字完全一致');
  $('sm_dt').value = '0.01'; $('sm_nsteps').value = '600';
  $('btnSolve').click();
  await waitFor(() => $('status').textContent.includes('非線性地震反應分析完成') && ev(`seismicResult.n_steps`) === 600, 'PEER NGA地震歷程結果');
  console.log('PEER NGA檔案上傳: 讀檔、送出分析、Solve成功 OK');

  // ---- 後端錯誤: 拿掉塑鉸容量 -> 清楚訊息, 舊結果不被破壞 ----
  const backup = ev(`JSON.stringify(model.members)`);
  ev(`model.members.forEach(m => { m.Mp_i = null; m.Mp_j = null; });`);
  $('sm_gm_type').value = 'pulse'; $('sm_gm_type').dispatchEvent(new w.Event('change'));
  const btnAnim = d.querySelector('.viewBtn[data-view="seismic_anim"]');
  btnAnim.click(); await sleep(20);
  $('btnSolve').click();
  await waitFor(() => $('status').textContent.includes('非線性地震反應分析失敗'), '錯誤訊息');
  assert($('status').textContent.includes('塑鉸容量'), '應說明沒有塑鉸容量: ' + $('status').textContent);
  const nStepsBeforeFailure = ev(`seismicResult.n_steps`);   // PEER NGA那次求解之後seismicResult已經更新過, 不能再拿最早的r比對
  assert(nStepsBeforeFailure === 600, '失敗時不應清掉上一次(PEER NGA)的結果');
  console.log('沒有塑鉸容量: 顯示清楚訊息 OK ->', $('status').textContent.slice(0, 50));
  ev(`model.members = ${backup};`);

  // ---- 新建/載入重置 ----
  ev(`resetViewToStructure();`);
  assert(ev(`seismicResult`) === null && $('btnViewSeisAnim').style.display === 'none'
         && $('seismicConfigBar').style.display === 'none' && $('seismicPanel').style.display === 'none'
         && $('seisPlaybackBar').style.display === 'none', '重置應清掉結果、按鈕、設定列、面板與播放列');
  assert(ev(`viewMode`) === 'structure' && $('analysisType').value === 'linear', '回到結構視圖與線性分析');
  console.log('新建/載入重置 OK');

  assert(errs.length === 0, '頁面有 JS 錯誤: ' + errs.slice(0, 3).join(' | '));
  console.log('\n全部通過: 網頁非線性地震反應端到端(匯入模型、設定、Solve、動畫播放、時程圖、遲滯迴圈圖、塑鉸圖、結果表、單位切換、自訂地震歷程、錯誤訊息、重置)');
  process.exit(0);
})().catch(e => { console.error('FAILED:', e.message); process.exit(1); });
