# frame2d 開發規劃

這份文件記錄「接下來要做什麼、為什麼」。README.md 專心回答「現在能幹嘛、
怎麼用」,兩份文件分工。

## 目前狀態 (已完成)

> 2026-09-21 更新: 同步實際程式碼現況。這份文件原本寫「非線性不在範疇內」,
> 但 P-Delta / 塑鉸 / pushover / corotational 已經在 `feature/geometric-stiffness`
> 完成並合併進 main。各求解函式的對照以
> [ANALYSIS_ARCHITECTURE.md](ANALYSIS_ARCHITECTURE.md) 為準。

版本標記: `v1.0-linear-elastic`(線彈性靜力核心)、`v1.1-nonlinear-static`
(加入 P-Delta / 塑鉸 / pushover / corotational / 溫度載重 / EqualDOF / Web GUI)。

```
frame2d
├── Model:    Node / Section / Member(frame/truss/cable, release, 塑鉸容量Mp) /
│             Support(None=自由 / 0.0=固定 / 數值=指定位移) / EqualDOF
├── Elements: 樑柱(彎曲+軸向) / 桁架(僅軸向,可拉可壓) / 纜線(僅軸向,只受拉+自動鬆弛迭代) /
│             含塑鉸樑元素(hinge.py) / co-rotational元素(corotational.py)
├── Loads:    節點集中力/力矩 / 桿件全長與局部段均佈載重(local, global_y, global任意角度,
│             線性變化) / 桿件內部集中力與力矩 / 分佈力矩 / 溫度效應(均勻+梯度)
├── Solvers:  solve()(DOFManager, 主要) / solve_condensation()(靜力凝縮, 參考實作) /
│             solve_pdelta()(線性化P-Delta) / solve_with_hinges()
├── 非線性靜力: run_pushover / run_pushover_converged / run_pushover_newton /
│             run_pushover_corotational_oneshot, 統一入口 analyze_pushover(geometry=, solver=)
├── Postprocess: N/V/M圖、變形圖、反力、塑鉸狀態、容量曲線
├── GUI:      webapi/(FastAPI) / webapi_stdlib/(純標準庫, Termux用)
└── Validation: 解析解 / sd_framework+anastruct(11案例46數值) / SW FEA(third-party) /
                桁架節點法手算 / 懸臂梁斜張橋(位移法vs力量法交叉驗證) /
                OpenSeesPy(塑鉸機制與Kg, 經portal-frame-pushover-scratch驗證後移植)
```

**目前已知限制(動力分析路線會逐一處理, 見文末):**

- 塑鉸是單向的(彈性→降伏, `HingeState.check_yield`), 沒有卸載/循環行為
- Newton / co-rotational 路徑目前不支援 release 端, 而且是 pushover 驅動
  (`prescribed_dofs`、`target_total`), 不是時間積分
- EqualDOF 用懲罰法(高勁度彈簧), 不是精確的自由度消去
- 完全沒有質量、阻尼、特徵值、時間積分, 也就是還沒有任何動力分析

這已經是一個能處理「frame + truss + cable 混合結構」的通用矩陣位移法核心,
不再是單純的驗證腳本集合。**不再繼續擴張「結構種類」**(不做第二種斜張橋案例、
不做更多bridge專屬功能), 力氣放在讓「載重/邊界條件/元素/後處理」這四個核心
抽象更完整、更穩定, 以及下面的動力分析路線。

## Load System v2 (已全部完成)

SW FEA 的載重介面對照 frame2d 現況:

| SW FEA 類型 | frame2d 現況 |
|---|---|
| Nodal Point Load | ✅ `point_load(node, fx, fy, m)` |
| Distributed(全長, 垂直桿件) | ✅ `distributed_load(member, w_start, w_end)`, 含線性變化 |
| Point Load(桿件中間任意位置) | ✅ `member_point_load(member, a, fx, fy, m)` (Phase 1) |
| Moment(桿件中間任意位置) | ✅ 同上, `m=` 參數 (Phase 1) |
| Distributed(局部段, 不用整根桿件) | ✅ `x_start` / `x_end` (Phase 2) |
| Support Displacement(強制位移/沉陷) | ✅ `Support` 改為指定值 (Phase 3) |
| 角度/局部-全域座標系選項 | ✅ `direction='local'` / `'global_y'` / `'global'` + `angle_deg` |
| 分佈力矩、溫度效應 | ✅ `distributed_moment()` / `thermal_load()` |

**實作順序(以下各階段皆已完成, 保留為實作記錄):**

### Phase 1: 桿件中間的集中力 + 集中力矩 ✅ 已完成
最常用、也是驗證「element load → equivalent nodal load」這個矩陣位移法
核心模式的最佳題目。做法照現有 `distributed_load` 的模式:算出固定端反力
→ 等效節點載重疊加進F → 回代時扣回來。**不新增節點**(不切割member,
不污染拓樸),用等效節點載重表示,這樣M(x)/V(x)公式也可以直接沿用現有的
`member_internal_forces`架構(局部段落分開處理)。

**實作記錄**: `member_point_load(member, a, fx, fy, m)`,固定端反力公式用
sympy從梁的微分方程(EI*v''=M(x)分段雙重積分)直接推導,不是憑記憶抄書
(避免重蹈`postprocess.py`的M(x)公式曾經正負號抄反的錯誤)。驗證優先用
簡支梁的「決定性反力」(純靜力學R=Pb/L等, 完全不依賴固定端反力公式本身
對不對)當第一層基準,再用懸臂梁解析解撓度公式交叉確認。四個案例(橫向
點載重反力、懸臂梁撓度、點力矩反力、軸向點載重)全數一次通過,機器精度。
`postprocess.member_internal_forces`跟`plotting.plot_loads`也都更新支援,
N/V/M圖能正確畫出跳躍不連續(剪力跳躍、彎矩斜率變化、彎矩跳躍三種情況都
驗證過)。見`tests/test_member_point_load.py`。

### Phase 2: 局部段均佈載重 ✅ 已完成
`distributed_load()` 的 `DistributedLoad` dataclass 加上了 `x_start`,
`x_end` 兩個欄位(預設None=整根桿件, 向下相容現有呼叫方式)。實作方式是
Phase 1點載重公式沿桿長積分的推廣: 不手動謄寫龐大的sympy封閉式展開式
(降低抄寫出錯風險), 改用6點高斯積分對`fixed_end_forces_point_load()`
在`[x_start,x_end]`區間積分, 對這種低次多項式被積函數是機器精度的數值
精確解。退化情況(x_start=0,x_end=L)已驗證跟既有全長UDL公式精確一致
(誤差~1e-14), 局部段案例(均佈+梯形兩種)用簡支梁決定性反力交叉驗證。
`postprocess.member_internal_forces`/`plotting.plot_loads`也都更新支援
多筆、各自範圍不同的局部段載重疊加。見`tests/test_partial_udl.py`。

### Phase 3: Support改成「指定值」而非布林值 ✅ 已完成
`Support(node, ux=True/False)` → `Support(node, ux=None/0.0/指定值)`,
`None`=自由, `0.0`=固定在原位, 非零值=強制位移(沉陷/施工誤差分析)。
這個改動統一了fixed/pin/roller/settlement, 不用另外做一個
`SupportDisplacement` API。核心公式: 邊界條件從「劃掉」變成
`K_ff u_f = F_f - K_fc u_c`(u_c為已知的支承位移向量), 改動集中在
`solve.py`, 不影響其他模組。

**實作記錄**: `fix()`/`pin()`/`roller_y()`維持原本行為(內部用0.0),
新增`support(node, ux=, uy=, rot=)`通用API可以直接設定沉陷量。驗證用
兩層: (1) 靜定結構支承沉陷應該零反力零內力(純剛體轉動, 不依賴任何公式,
是結構學最基本的性質), (2) 一次靜不定梁支承沉陷對照傾角變位法經典公式
M_A=-3EIΔ/L²(跟使用者自己的sd_framework同一種方法)。過程中因為
`Support`從布林值改成數值, 意外抓到`plotting.py`一個真的bug: 原本用
`if support.rot and ...`判斷是否為固定端, 但`0.0`在Python布林判斷裡是
假值, 導致固定端符號被誤判——已修正成`is not None`判斷。全部既有測試
(11個既有tests/*.py)重跑確認零回歸。見`tests/test_support_displacement.py`。

### Phase 4: Internal hinge / element release ✅ 已完成 (比原規劃簡單)
原本規劃「需要先做DOFManager升級才能開始」, 但實際做的時候發現有更簡單的
路線: 用靜力凝縮(static condensation)處理單一桿件的局部勁度矩陣就好
(標準4EI/L,6EI/L²改成3EI/L,3EI/L², 釋放端那列/行全為0), **完全不用動
`dofs_of()`的DOF系統**——這是結構分析裡處理桿端鉸接的標準經典做法,
之前判斷「需要先升級DOFManager」是過度保守。

**實作記錄**: `add_member(..., release_i=, release_j=)`。均佈載重的固定端
反力公式一開始正負號推導錯誤(用Gerber梁純靜力學驗證時抓到, 單獨測試
給出45/15不是預期的簡支梁30/30), 修正後兩個方向都通過獨立驗證。驗證分
兩層: 古典Gerber梁(兩跨連續梁中間放鉸接, 從一次靜不定變成靜定, 拆解成
獨立簡支梁純靜力學驗證)+ 對照SW FEA門型鋼架案例(3桿件x11點BM逐點比對)。
見`tests/test_element_release.py`、`tests/test_element_release_vs_swfea.py`。

**目前限制**(下一步如有需求再擴充): 均佈載重只支援整根桿件+真正均佈
(不支援線性變化、局部段); 桿件內部集中力還不支援; 兩端同時釋放的固定端
反力公式沒處理(等同truss但沒有內部載重公式)。**這些限制只在主要求解器
(靜力凝縮版本)上。**

### Phase 4b: DOFManager版本 ✅ 已完成, 並已升格為主要求解器
在`frame2d/dofmanager.py`實作了完全獨立的第二套求解路徑(不呼叫solve.py
任何邏輯): 讓每個release端擁有自己專屬、不共用的轉角自由度, 而不是靠
靜力凝縮修剪局部勁度矩陣。這樣就可以直接用「標準」(未修改過的)局部
勁度矩陣跟固定端反力公式, 鉸接的物理效果完全靠DOF不共用來體現(解聯立
方程組時, 那個自由度沒有其他東西連著, 自然浮現M=0)。

這在數學上跟靜力凝縮是等價的(高斯消去法, 只是消去時機不同), 兩者答案
在Gerber梁跟門型鋼架案例上都精確一致(誤差~1e-13~1e-15, 浮點精度等級) ——
交叉驗證抓的是「兩邊各自實作有沒有bug」, 不是驗證哪個物理模型才對。

**額外好處**: 因為不需要「release專屬固定端反力公式」, DOFManager版本
可以直接支援靜力凝縮版本不支援的組合(桿件內部集中力加在release桿件上),
已用節點分割法獨立驗證過正確性。見`tests/test_dofmanager_vs_condensation.py`。

**升格記錄**(2026-08-29): `frame2d/__init__.py`的`solve`這個公開名字
現在指向`solve_dofmanager`(補上cable鬆弛迭代+跟`solve_condensation()`
統一輸出格式`SolveResult`/`MemberResult`後才升格, 這兩件事之前是缺的,
不是單純改名)。原本的靜力凝縮實作改名`solve_condensation()`, 保留當
參考/回歸測試實作。**升格後全部27個既有測試(含所有SW FEA逐點比對案例)
重跑一次確認零回歸**, 這等於是拿DOFManager重新過一次專案裡最嚴格的
驗證關卡。

同一批改動也順手把`dofs_of()`底層從`3*node_id`直接硬編碼, 改成
`node_id -> 緊湊索引`的對照表(`Frame2D.node_index()`), 這樣node_id
不用連續、不用從0開始(例如10,25,99這種id, 3個節點只佔用9個DOF, 不是
300個)——因為`dofs_of()`本來就是設計成呼叫端不用管底層怎麼編號的抽象
邊界(見model.py開頭的原始設計註解), 這次換掉底層實作完全不用改任何
呼叫端。

順便發現並修正一個真的問題: 17份測試檔案裡有16份是寫成「直接執行的
腳本」風格(print+assert), 不是pytest的`def test_*()`函式, 導致單純
執行`pytest`(不加參數)只會收集到11個測試, 漏掉16個檔案的驗證。已加
`tests/test_zz_all_script_style_tests.py`當wrapper(用subprocess跑每支
腳本、檢查exit code), 讓`pytest`能收集並執行全部27個測試。

## 暫緩/刻意不做的項目

- 斜張橋以外的其他特殊橋型/結構型式: frame+truss+cable已經證明夠通用,
  不需要再做專屬功能證明
- 勁度/強度劣化、捏縮(pinching)、等向硬化: 循環塑鉸(D7)目前只有雙線性運動硬化(等同 Steel01 預設),
  這幾項還沒做
- 動力分析初期(D1~D5)不做: 多支承激振、顯式時間積分、非古典阻尼、cable
  參與動力分析、EqualDOF 參與動力分析(見下方「既有限制」表格)
- 3D、土-結構互制: 不在範疇內

## 動力分析路線 (Dynamics Roadmap)

> 狀態: **規劃中, 尚未開始實作。** 開發分支: `feature/dynamics`(從 main 開出)。
> 制定日期 2026-09-21, 依據是對 v1.1 程式碼逐檔盤點的結果。階段編號用 D0~D8,
> 避免跟上面 Load System v2 的 Phase 編號混淆。

### 定位

frame2d 不做「小型 OpenSeesPy」。OpenSees 是成熟的通用非線性有限元素框架;
frame2d 的價值是**透明**: `K → M → C → eigen → RSA → transient → nonlinear
transient`, 每一步的矩陣與演算法都看得到、都有解析解或第三方可以對照。
動力分析建立在既有的靜力與非線性核心之上, 不另開新 repo。

「線性/非線性」是模型性質, 不是分析類型:

| | 線性 | 非線性 |
|---|---|---|
| 靜力 | `solve()` | `solve_pdelta()` / pushover / corotational Newton (已完成) |
| 模態 | eigen (D2) | 不做(需先定義線性化的基準狀態) |
| 反應譜 | RSA (D3) | 不適用 |
| 時程 | 線性時程 (D4, D5) | 非線性時程 (D6~D8) |

反應譜分析算結構動力分析, 是完全建立在模態上的線性方法, 給的是各模態最大反應
的估計, 不是時間歷程。

### 設計原則

1. **只加不改。** 沿用 `analyze_pushover()` 的作法: 新增函式、新增檔案, 不修改
   既有求解器。(D0 原本規劃要重構 `_solve_once_dofmanager`, 實際檢查後發現不需要,
   改成公開薄包裝, 所以目前沒有任何例外。)
2. **明確拒絕, 不靜默退化。** 動力分析遇到尚未支援的組合(cable、EqualDOF、
   非零指定位移等)一律 `raise ValueError`, 跟現有風格一致。
3. **每個階段三件套:** (a) 解析解或手算基準 (b) 跟 OpenSeesPy 交叉驗證
   (c) 既有全部測試零回歸。
4. **檔案維持平坦。** 不預先拆 `static/`、`dynamics/` 子資料夾, 新功能放新檔案
   (`assembly.py`、`mass.py`、`modal.py`、`spectrum.py`、`damping.py`、
   `newmark.py`、`excitation.py`), 等 D4~D5 穩定後再評估是否整理目錄。
5. **數值判斷不能依賴「恰好為零」。** 機構/奇異的偵測一律用「以彈性(或凝縮前)對角縮放後的
   特徵值」判斷, 不靠 LAPACK 遇到恰好零主元才丟例外。x86 上剛好會丟, Termux(arm64、
   Python 3.14)上同樣的運算留下 1e-13 的殘餘就不丟了。新的機構偵測都要有「加雜訊」的迴歸測試
   (見 `test_force_control_roundoff.py`、`test_modal.py` 層6)。
6. **只依賴 numpy。** Termux / Pydroid3 環境裝 scipy 不一定順利。集中質量用
   `M^(-1/2) K M^(-1/2)` 加 `np.linalg.eigh`, 一致質量用 Cholesky 化成標準
   特徵值問題。教學尺度的稠密矩陣完全夠用。

### 動力分析前必須處理的既有限制

| # | 問題 | 現況位置 | 處理方式 | 何時 |
|---|---|---|---|---|
| 1 | K 沒有**公開**的組裝入口 | `pushover._assemble_stiffness_with_hinges` 已經是純組裝函式(不含載重/邊界條件), 但是私有、名字綁 pushover; `dofmanager._solve_once_dofmanager` 內嵌另一份; `newton._assemble_global` 是 corotational 版 | 新增 `assembly.assemble_K` 薄包裝, **不改任何既有求解器**; 測試守住兩份組裝不分歧 | D0 ✅ |
| 2 | Section 沒有質量 | `model.Section` 只有 E, I, A | 加 `Section.rho`(質量密度)與 `Frame2D.add_mass(node, mx, my, Iz)`, 預設 None 不影響既有行為 | D1 ✅ |
| 3 | 單位制 | 核心與網頁後端 | **核心與單位無關, 只要求一致: 質量單位 = 力單位·s²/長度單位。** 網頁後端固定 SI(Pa, m, N → 質量 kg); kN、m 手算配 ton。網頁「單位設定」新增質量/轉動慣量/密度/加速度/速度選單。質量單位寫錯不會報錯, 只會讓頻率差常數倍, 所以每個動力測試都要有解析解量級檢核 | D1 ✅ |
| 4 | 無質量自由度 | `build_dof_map` 的 release 專屬轉角 DOF; 集中質量沒有轉動慣量; truss 節點轉角 | 靜力凝縮 `K_eff = K_dd − K_dm K_mm⁻¹ K_md`(對無質量 DOF 是精確的) | D1~D2 |
| 5 | EqualDOF 是懲罰法 | `dofmanager._apply_equal_dof`(1e6 倍最大對角項) | D1~D5: 有 `equal_dofs` 就明確報錯; 之後改成精確的 master-slave 消去。懲罰彈簧會造成虛假高頻模態與條件數惡化 | D1 / 之後 |
| 6 | Cable 鬆弛迭代是狀態相依的 | `dofmanager` 的鬆弛迭代 | D1~D5: 有 cable 就明確報錯。模態分析要先定義基準狀態(例如重力下拉緊); 時程中鬆弛需要事件處理 | D6 之後 |
| 7 | 非零指定位移 | `Support` 指定值 | 動力分析只接受 `None` / `0.0`; 地震輸入一律用等效力 `−M r a_g`, 不用支承位移 | D1 |
| 8 | 塑鉸單向 | `hinge.HingeState.check_yield` | 新增循環塑鉸類別 `CyclicHingeState`(`cyclic.py`), 不改 `HingeState` | D7 ✅ |
| 9 | Newton 路徑是 pushover 驅動且不支援 release | `newton.py`(`_check_no_releases`) | 新增時間積分外層迴圈, 重用 `corotational.py` 與 `_assemble_global`; 是否補 release 支援在 D6 決定 | D6 |
| 10 | 文件落後 | README「尚未支援」段落、`BENCHMARK_SUITE.md` 的測試數量 | 同步更新 | D0 |

### 階段規劃

依賴關係: `D0 → D1 → D2 → {D3, D4 → D5} → D6 → D8`; **D7(循環塑鉸)只需要靜力求解器, 不依賴 D3~D6, 已提前完成**, D6/D8 直接使用它。D3(反應譜)與 D4(線性
時程)彼此獨立, 可以互換順序。

| Stage | 內容 | 預計新檔 | 狀態 |
|---|---|---|---|
| D0 | 前置: 公開 `assemble_K`(薄包裝), 同步文件 | `assembly.py` | ✅ |
| D1 | 質量矩陣(集中 → 一致) + 動力單位 | `mass.py` | ✅ |
| D2 | 特徵值/模態分析 + 模態性質 | `modal.py` | ✅ |
| D2b | 網頁: 節點質量輸入、模態分析端點、結果表與振型顯示 | webapi | ✅ |
| D3 | 反應譜分析 (SRSS / CQC) | `spectrum.py` | ✅ |
| D3b | 網頁: 「反應譜」分析、規範/自訂反應譜、反應譜曲線圖、模態/位移/桿件內力結果表 | webapi | ✅ |
| D4 | 線性時程(Newmark)+諧和/脈衝/任意力輸入 | `newmark.py`, `excitation.py` | ✅ |
| D5 | Rayleigh阻尼 + 地面加速度輸入 | `damping.py`, `excitation.py`(擴充) | ✅ |
| D6 | 非線性時程(Newmark + 循環塑鉸, event-to-event) | `nonlinear_newmark.py` | ✅ |
| D7 | 循環塑鉸(遲滯) + 準靜態反覆載重(**提前**, 見下) | `cyclic.py` | ✅ |
| D7b | 網頁: 「循環」分析、遲滯迴圈圖、塑鉸 M-θp 圖、每級耗能與等效阻尼比 | webapi | ✅ |
| D8 | 非線性地震反應一站式入口 + 能量平衡診斷 | `seismic.py` | ✅ |
| D6b | 網頁: 非線性地震反應(選地震歷程、地震動畫播放、時程/遲滯迴圈/塑鉸M-θp圖) | webapi | ✅ |
| D9 | 真實地震紀錄支援(PEER NGA .AT2格式) | `ground_motion_io.py` | ✅ |
| D10 | 網頁匯出Markdown: 模態/反應譜/循環/非線性地震 | webapi/static/index.html | ✅(Markdown; PDF規劃中) |
| D11 | 網頁匯出PDF: 模態(反應譜/循環/非線性地震規劃中) | webapi/pdf_export.py | ✅(模態; 其餘規劃中) |

#### D0 前置: 公開 assemble_K ✅ 已完成

- 新增 `frame2d/assembly.py`: `assemble_K(frame, hinge_states=None, axial_forces=None)`
  回傳 `Assembly(K, member_dofs, member_T, member_L, n_node_dof, n_extra_dof)`
  (NamedTuple, 可用欄位名稱或直接拆包)。內部只是呼叫
  `pushover._assemble_stiffness_with_hinges()` 加 `build_dof_map()`, 沒有新的組裝邏輯
- **規劃更正**: 原本寫「K 沒有獨立的組裝函式, 要重構 `_solve_once_dofmanager`」,
  這個說法太重。`_assemble_stiffness_with_hinges` 早就是純組裝函式, 缺的只是公開入口。
  所以改成薄包裝, 沒有動任何既有求解器, 也就不需要逐位元回歸 baseline
- 代價: 專案裡有兩份獨立寫的 K 組裝(`assemble_K` 呼叫的那份, 跟
  `_solve_once_dofmanager` 內嵌的那份)。分歧風險由 `tests/test_assembly.py` 守住
- 驗證分兩層, **兩層都需要**:
  - 層1(不依賴任何求解器): 無支承、無 release 的剛架, K 剛好 3 個零特徵值, 且
    `K @ 剛體位移 = 0`
  - 層2(對照既有求解器): 用 `assemble_K` 自己組 K、自己劃分邊界條件、自己解,
    對照 `solve()`、`_solve_once_dofmanager(member_axial=)`、`solve_with_hinges()`;
    涵蓋一般剛架、release(額外 DOF)、桁架、equalDOF、P-Delta 幾何勁度、含塑鉸勁度,
    位移差為 0.0(逐位元一致)
- 測試有效性(突變檢查): 故意讓 pushover 那份組裝的 A 偏 1%, 層2 抓到(位移相對差
  1.5e-4), 層1 抓不到(剛體模態不受 A 影響)——所以兩層缺一不可
- 已知範圍: 不處理鬆弛 cable(一律當 taut); equalDOF 仍是懲罰法。兩者在 D1~D5 由
  動力分析入口明確拒絕
- 文件同步: README 的範疇與結構、`BENCHMARK_SUITE.md` 的測試數量與未編號測試

#### D1 質量矩陣 ✅ 已完成

- `Section.rho`(質量密度, 預設 None=無分佈質量)、`Frame2D.add_mass(node, mx, my, Iz)`
  (同節點可多次呼叫, 會加總)
- `frame2d/mass.py`: `assemble_M(frame, kind='lumped' | 'consistent')`, DOF 編號與
  `assemble_K` 完全相同; 另有 `influence_vector(frame, 'x'|'y')`、`total_mass()`、
  `consistent_mass_local()`
  - 集中質量: 桿件質量各半加到兩端 ux, uy; 不計轉動慣量; 矩陣對角
  - 一致質量: Euler-Bernoulli 局部 6×6(軸向線性、橫向 Hermite 三次), 桁架用 4×4
    平動質量; release 端的專屬 DOF 放標準 6×6 項(一致有限元素對鉸接端的正確處理)
- 明確拒絕(`ValueError`): cable、equal_dof、非零指定支承位移、非法 kind、負密度/負質量、
  add_mass 指到不存在的節點
- **單位(規劃更正)**: 原本寫「決定 kN-m-ton-s」, 但網頁後端其實固定用 SI(Pa, m, N),
  質量單位應該是 kg, 不是 ton。核心對單位不做任何假設, 只要求一致(質量 = 力·s²/長度)。
  網頁「單位設定」新增: 質量(kg, t)、質量轉動慣量(kg·m², t·m²)、密度(kg/m³, t/m³)、
  加速度(m/s², g, gal)、速度(m/s, cm/s, mm/s), 預設 t、t·m²、t/m³、g、cm/s(跟預設的
  kN 自洽: 1 kN = 1 t·m/s²)。斷面管理面板新增 ρ 輸入; 節點質量與加速度/速度單位目前
  還沒有輸入或顯示的地方, 留給 D2(模態結果)與 D3~D5(反應譜、時程)接上
- 順手修正一個既有的網頁問題: 斷面管理面板開著時換單位(E/I/A/密度), 畫面上的數字與標籤
  不會重畫, 但按「套用」時會用新單位換算, 會靜默存成錯的值。現在換單位會重畫該面板
- 驗證(`tests/test_mass.py`):
  - 層1 用 sympy 從形函數積分獨立推導局部一致質量矩陣, 相對差 6e-17
  - 層2 不依賴係數的解析不變量: 總質量 rᵀMr = ΣρAL + 節點質量(含斜桿、桁架、release);
    繞原點剛體轉動的極慣量, 一致質量精確等於 ΣρAL(|c|²+L²/12), 集中質量等於
    ΣρAL(|c|²+L²/4)(剛體轉動位移場在形函數空間內, 所以一致質量是精確的, 不是近似)
  - 層3 對稱、半正定、集中質量對角
  - 層4 頂端質量懸臂柱靜力凝縮後 ω²=3EI/(mL³); 同一個物理結構分別用 SI(Pa, kg) 與
    kN·m·ton 兩套單位算, ω 完全相同; 反面案例(SI 力單位配 2 而不是 2000 kg)差 √1000 倍
  - 突變檢查: 一致質量係數 13L 抄成 12L 被層1 抓到; 座標轉換方向反了被層2 極慣量抓到
- 單位表同步(`tests/test_dynamic_units.py`): 單位換算表同時存在於兩份 index.html 與兩份
  pdf_export.py, 這支測試守住四份一致、選單/unitOf/持久化/匯出 payload 都掛上、
  單位系統自洽(kN↔t、N↔kg、7.85 t/m³=7850 kg/m³), 並驗證 rho 從網頁 JSON 抵達
  `Section.rho`

#### D2 特徵值與模態性質 ✅ 已完成(核心; 網頁部分見 D2b)

- `frame2d/modal.py`: `eigen(frame, n_modes=None, mass='lumped'|'consistent', zero_tol=1e-10)`
  回傳 `Modal`: `omega`(rad/s)、`period`、`frequency`、`phi`(全DOF、對M正規化)、
  `gamma`、`eff_mass`、`cum_ratio`、`cum_ratio_total`、`mass_free`、`mass_total`,
  以及 `node_shape()`、`modes_needed(ratio, direction)`、`table()`
- 做法(只用 numpy): 對無質量DOF靜力凝縮(精確, 不是Guyan近似) → 對角縮放 →
  Cholesky 化成標準特徵值問題 → 回填無質量DOF。全部細節寫在 `modal.py` 開頭
- **對角縮放的價值**: 結構勁度常常「梯度式」(A 預設 1e8 讓軸向勁度比彎曲大 1e13 倍)。
  不縮放時 eigh 的絕對誤差 ~ eps·‖K‖ 會吃掉低階模態; 縮放後 ω² 跨度 4e12 的案例
  彎曲與軸向頻率都精確到 2e-16
- **機構偵測不能用 ω² 的相對大小**: 我最初用「最小 ω² / 最大 ω²」判零頻率, 結果把上述
  合法但跨度大的結構誤判成機構。改成看「縮放後(單位對角)勁度矩陣的最小特徵值」,
  健康結構是 O(1), 機構/剛體模態趨近 0
- **質量參與的定義**: Γ = φᵀ M r 用完整地面運動向量 r(含支承DOF), 一致質量與支承DOF
  的耦合項也算進去。`mass_free` 定義為取全部模態時 ΣM* 的精確值(動力完備性), 集中質量下
  等於總質量扣掉支承DOF上的節點質量。`cum_ratio` 相對可動質量(取全部模態必為1),
  `cum_ratio_total` 相對總質量(規範檢核「≥90% 總質量」用這個, 較保守)
- 驗證(`tests/test_modal.py`, 只依賴 numpy/sympy 之外的解析解):
  - 頂端質量懸臂 ω²=3EI/(mL³); 兩質量懸臂(無質量梁)的精確柔度解, 頻率到 1e-11、|Γ| 到 1e-10
  - 均質懸臂 (β_nL)²√(EI/ρAL⁴): 一致質量從上方收斂、第1模態誤差比 15.7→15.9→16.8
    (理論 16 = O(h⁴)); 集中質量從下方收斂、誤差比 3.91→3.98→3.99(理論 4 = O(h²))
  - 簡支梁 (nπ)²√(EI/ρAL⁴); 懸臂第1模態有效質量比 0.61308(精確積分)對有限元素 6.6e-10
  - 正交性 φᵀMφ=I、φᵀKφ=diag(ω²) 到 1e-15; 完備性 ΣM*=可動質量, 且與獨立的 Schur 補算法一致
  - 靜力凝縮 vs 給每個轉角加極小轉動慣量 ε(質量正定、不需凝縮的獨立路徑): 誤差隨 ε 線性
    趨近 0(4e-4 → 4e-6 → 8e-8)
  - release 端專屬DOF: 「固定+兩端release」與「鉸支承+滾支承」特徵值相同(集中 7e-15、一致 2e-11)
  - 單位一致性: 同一結構 SI(Pa, kg) 與 kN·m·ton 的 ω 相對差 2e-13
  - 突變檢查: 拿掉縮放回復、凝縮符號反、Γ 用只含自由DOF的 r、不回填無質量DOF, 四個都被抓到
- **對 OpenSeesPy 交叉驗證**(`tests/test_modal_vs_openseespy.py`, 選用, 沒裝就 SKIPPED,
  OpenSeesPy 3.8.0): 兩個模型(8 元素懸臂、含斜桿/桁架/轉動慣量的不對稱混合剛架)×
  集中/一致質量: 特徵值相對誤差 3e-13、前 6 個模態 MAC = 1.0、總質量一致; 集中質量的可動質量、
  |Γ|、有效質量逐項一致
- **機構偵測必須耐得住捨入雜訊**(在 Termux / arm64 / Python 3.14 上實際踩到同類問題後補上):
  理論上「恰好為零」的勁度, 換平台(FMA、BLAS 版本)會變成 1e-14 量級的殘餘。第一版縮放參考用
  凝縮後自己的對角, 結果鉸支承懸臂在 1e-16 雜訊下就回傳 ω=6.7e-7 的假模態而不是報錯——勁度掉到
  雜訊量級的 DOF 被自己的對角放大成單位對角, 看起來反而「健康」。現在縮放參考改用**凝縮前**的
  對角 K_dd,ii。`test_modal.py` 層6 把 assemble_K 換成加雜訊版本(1e-16~1e-13, ∝√(對角_i·對角_j)),
  三種機構(無支承、鉸支承懸臂、雙鉸接質量節點)×5 個雜訊量級×4 組全部正確報錯(60/60), 合法結構
  (含 A=1e8 的軸向極剛懸臂)頻率擾動 ≤ 5e-11; 舊縮放版本會被這個測試抓到
- 限制: 一致質量的參與係數與 OpenSees `modalProperties` 不逐項一致(第1模態差 0.3~0.5%):
  它的可動質量在一致質量下等於集中質量的值, 定義與這裡不同, 我沒有完整重現。這裡的定義有
  連續體解(0.61308)背書, 差異只列出、不宣稱相等

#### D2b 網頁模態分析 ✅ 已完成

- **使用方式**: 斷面管理輸入密度 ρ、點節點在屬性面板輸入質量 mx/my/Iz(用「單位設定」選的質量/
  轉動慣量單位, 存成 SI 的 kg/kg·m²; 有質量的節點旁邊會標紫色 m) → 上方「分析」選「模態」→ 設定
  模態數與質量矩陣(集中/一致) → Solve。結果: 結果表(T、f、ω、有效質量、質量比、累積質量比)、
  點選一列切換畫面上的振型。質量比 = 相對可動質量; 「累積(總)」= 相對總質量, 規範檢核「≥90%」看這欄
- 核心: `Modal.member_curves(mode, scale, n)` 用桿端局部位移(含 release 端專屬 DOF)+ Hermite 三次
  形函數畫每根桿件的變形曲線(truss 維持直線); `modal_to_dict()` 轉成 JSON, **兩個後端共用**, 不各自
  重複實作。後端: FastAPI 與 stdlib 各一個 `POST /modal`(`NodeIn` 新增 mx/my/Iz、`FrameIn` 新增
  modal_n_modes/modal_mass_kind), 錯誤(無支承、沒有質量、cable、equalDOF、非零支座沉陷)都是清楚的
  400 訊息; `/solve` 完全不受影響
- 驗證:
  - `tests/test_modal_shape.py`: 剛體平移/轉動位移場(含斜桿與桁架)下曲線精確等於解析位置
    (8.9e-16 / 2.2e-15, 同時驗證局部↔全域座標轉換); 端點=節點位置+位移且共用節點連續; 16 元素懸臂
    第 1 模態曲線對精確 cosh/cos 形狀最大差 4.7e-7; release 端斜率=專屬 DOF 轉角(不是節點轉角);
    JSON 無 nan、放大倍率正確。突變檢查(轉換符號、Hermite 係數、用節點轉角取代專屬 DOF)三個都抓到
  - `tests/test_web_modal_api.py`: stdlib 後端 `/modal` 與核心直接算逐項相同; FastAPI 與 stdlib 回傳
    完全相同; 各種錯誤情況; 舊呼叫端(沒有質量欄位)不受影響
  - `tests/test_web_modal_e2e.py`(選用, 需要 node + jsdom): 用 jsdom 載入**真正由後端提供**的
    index.html, `window.fetch` 接到真的 Python 後端, 像使用者一樣操作: 節點質量輸入(5 t → 5000 kg、
    負數被拒絕)→ 選「模態」→ Solve → 結果表 3 列、振型曲線 3 條、週期與核心逐項一致(1e-9)→ 點選
    切換模態 → 切換質量單位(有效質量 10.831 t → 10831 kg)→ 拿掉支承得到清楚的機構訊息且不破壞舊結果
    → 新建/載入重置。突變檢查(切單位不重畫、用錯單位種類、節點面板忘了存質量)三個都抓到
- **限制(誠實)**: jsdom 不是真的瀏覽器——沒有版面/CSS 渲染、沒有真實觸控事件, 手機瀏覽器上的外觀與
  操作手感要你實際開網頁確認。振型目前是靜態顯示(沒有動畫); 模態結果尚未納入 PDF 匯出

#### D3 反應譜分析 ✅ 已完成(核心; 網頁 D3b 之後再做)

- `frame2d/spectrum.py`: `response_spectrum(modal, spectrum, direction='x'|'y', damping=0.05,
  combine='SRSS'|'CQC', n_modes=None)`, `spectrum` 只吃 `S_a(T)` 的 callable(單位跟 modal 一致,
  核心不做單位換算, 台灣規範譜的查表/公式做成 adapter, 放使用端不進核心)
- 做法: 每個模態的位移 `D_i = Γ_i · S_d(T_i) · φ_i` 精確滿足 `K D_i = f_i`(等效靜力
  `f_i = Γ_i · S_a(T_i) · M φ_i`, 這是把 `K φ_i = ω_i² M φ_i` 代入 `K D_i` 直接得到的恆等式,
  **不需要另外解一次線性方程式**); 桿件內力用局部勁度矩陣直接算(純線性彈性); 反力
  `R_i = K D_i − f_i`(跟 `SolveResult.reactions` 同一套 `R=Ku-F` 慣例)。**每個反應量(位移、
  桿件內力、反力、基底剪力)各自算出各模態的訊號值(有正負號)後再組合**, 不是先組合位移
  再反推內力/反力——`RSAResult` 的 `modal_*` 系列欄位就是組合前的訊號值
- 總基底剪力有不用組裝 K 的簡潔公式(標準結果, Chopra): `V_i = Γ_i² · S_a(T_i) = M*_i · S_a(T_i)`
  (D2 的有效模態質量直接可用); `response_spectrum()` 內部每次呼叫都會做一次自我一致性檢查
  (自由DOF上 `K@D_i-f_i` 必須精確為0), 不一致直接報錯而不是吃案
- SRSS: `combined=√(Σvᵢ²)`。CQC(Der Kiureghian & Rosenblueth): `ρᵢⱼ` 用單一阻尼比(全模態同一個,
  等 D5 Rayleigh 阻尼才有逐模態的阻尼比); `combine_cqc`/`combine_srss`/`cqc_rho` 也公開給外部
  客製化組合用
- 驗證(`tests/test_spectrum.py`, 只依賴 numpy):
  - SDOF 解析解(頂端質量懸臂柱): `S_d=A0/ω²`、基底剪力 `=m·A0`、反力慣例, 全部到 1e-12
  - 兩質量懸臂(無質量梁, 沿用 D2 測試的獨立柔度矩陣解析解): 各模態訊號位移、SRSS/CQC 組合位移
    對手算值(含 CQC 用同一組 `cqc_rho` 手算的相關係數矩陣), 到 1e-10
  - 反力/基底剪力的簽名恆等式 `Σ(反力·r) = -V_i`, 在多支承斜桿桁架混合剛架(集中/一致質量、
    x/y 兩方向, 共 4 組 × 5 個模態)驗證, 最大相對誤差 2.8e-13
  - 獨立 K 線性方程式交叉驗證: 不用 `K D_i=f_i` 的恆等式捷徑, 直接對自由DOF解一次線性方程式,
    確認等於 `D_i`(驗證「不用解方程式」這個捷徑本身是對的, 不是自我一致的巧合), 到 6.9e-14
  - CQC 性質: 對角=1、對稱, 對 Der Kiureghian 公式的**手算數值**(不是鬆的區間檢查)到 1e-12,
    頻率分離良好時 CQC 收斂到 SRSS(相對差 3.8e-4)
  - 突變檢查 6 個(等效靜力用 M·r 而非 M·φ、CQC 分母漏項、SRSS 漏平方、基底剪力用 Γ 而非 Γ²、
    S_d 公式漏平方、寬鬆區間檢查蓋不住的 CQC 分母錯誤), 全部被抓到——第一次寫的區間檢查
    (`ρ>0.3`)真的漏掉了一個分母錯誤, 換成手算數值後才抓到, 過程記錄在測試裡當提醒
- **對 OpenSeesPy 交叉驗證**(`tests/test_spectrum_vs_openseespy.py`, 選用, 沒裝就 SKIPPED):
  單一模態(1、2、3 個模態各自單獨算)的位移與反力, 位移/反力相對差 1e-13。**範圍限制**:
  `ops.responseSpectrumAnalysis` 的多模態 SRSS/CQC 組合語法在這個環境裡沒能摸清楚(試了幾種
  `-mode` 傳法, 結果都像只用了其中一個模態, 不像做了組合), 所以只比對單一模態, 多模態組合
  改由上面的兩質量懸臂獨立解析解驗證。集中質量下逐項相同; 一致質量沿用 D2 記錄過的同一個
  限制(參與係數定義跟 OpenSees `modalProperties` 不完全一致), 這裡的比對用集中質量
- 限制(誠實): `spectrum` 是純量函式(單一方向、單一譜), 沒有多方向組合(例如 30% 法則)、
  沒有考慮扭轉耦合以外的三維效應(本來就是 2D 框架); CQC 的阻尼比是全模態同一個純量

#### D3b 網頁反應譜分析 ✅ 已完成

- **使用方式**: 斷面填密度 ρ、節點填質量(跟模態分析一樣)→ 上方「分析」選「**反應譜**」→ 設定
  方向、組合方法(SRSS/CQC)、阻尼比、模態數(留空=全部)、質量矩陣 → 選反應譜來源:
  「規範(簡化)」填 S<sub>DS</sub>/S<sub>D1</sub>/T<sub>L</sub>(四段式標準形狀, **不是**精確
  工址查表版本), 或「自訂」用可增減列的週期-S<sub>a</sub> 表格(分段線性內插, 順序不拘會自動
  排序)→ Solve。結果: 反應譜曲線圖(紅點標模態落點)+ 結果面板(每模態的 T/f/S<sub>a</sub>/Γ/
  質量比/累積(總)/模態基底剪力、節點位移、桿件端點內力), 單位跟「單位設定」一致
- 核心新增(`frame2d/spectrum.py`): `taiwan_code_spectrum(SDS, SD1, TL)`(簡化四段式形狀,
  明確標注不是精確查表版)、`custom_spectrum(points)`(分段線性內插, 自動排序, 範圍外夾在端點
  值不外插)、`spectrum_analysis()`(網頁一站式入口: 建反應譜 callable → eigen() → 
  response_spectrum() → 取樣反應譜曲線給前端畫圖, 不用在 JS 重新實作反應譜公式)、
  `rsa_to_dict()`(JSON)
- 後端: FastAPI 與 stdlib 各一個 `POST /rsa`(`FrameIn` 新增 rsa_direction / rsa_damping /
  rsa_combine / rsa_n_modes / rsa_mass_kind / rsa_spectrum_type / rsa_code_sds / rsa_code_sd1 /
  rsa_code_tl / rsa_custom_points); 錯誤都是清楚的 400(缺 SDS/SD1、缺自訂點、沒有支承等, 後者
  直接沿用 `eigen()` 的訊息)
- 驗證:
  - `tests/test_spectrum_web.py`: `taiwan_code_spectrum` 四段式形狀(三個轉角週期連續、平台段
    恆等 S<sub>DS</sub>·g、長週期段 ∝1/T²)、`custom_spectrum` 對 `np.interp` 逐點比對與邊界外
    夾住(不外插)、`spectrum_analysis`+`rsa_to_dict` 跟直接呼叫 `eigen()+response_spectrum()`
    逐項相同、curve 取樣範圍涵蓋 T<sub>L</sub> 與模態週期、明確拒絕。突變檢查 4 個全部被抓到
  - `tests/test_web_rsa_api.py`: stdlib 後端對核心逐項相同(規範與自訂反應譜都測); FastAPI 與
    stdlib 回傳完全相同; 各種錯誤訊息
  - `tests/test_web_rsa_e2e.py`(選用, 需要 node + jsdom): 走「匯入JSON」載入
    `examples/portal_rsa_demo.json` → 選「反應譜」→ 輸入驗證(負 SDS、阻尼比超出範圍)→ Solve →
    網頁結果與核心逐項一致(含每個模態的參與係數 Γ)→ 曲線圖(300 點取樣、無 NaN)→ 結果表(**檢查
    渲染出來的表格文字, 不是只比對底層 JS 資料物件**——見下面的踩坑記錄)→ 單位切換 → 自訂反應譜
    的新增/編輯/刪除(含最少 2 點的前端即時檢查, 不用等後端回應)→ 拿掉支承得到清楚的機構訊息且
    不破壞舊結果 → 重置。突變檢查 4/4 被抓到
  - **踩坑記錄**: 第一版的 e2e 測試只比對 `rsaResult`(前端存的 JSON 資料)有沒有跟核心一致,
    沒有檢查表格實際渲染出來的文字。結果一個「Γ 欄位畫錯成模態編號」的突變沒被抓到——因為
    底層資料物件本身沒錯, 只有渲染那一步把欄位對錯了。加上「比對渲染出來的 `<td>` 文字」之後
    才抓到。這件事本身也提醒: 純粹比對資料物件不能取代檢查畫面上實際顯示的內容
  - 另外一次修測試檔時因為文字替換操作沒對齊, 把整個 `rsa_e2e.js` 意外清空成 0 bytes, 靠重新
    整份寫入救回來——這是純粹的工具操作失誤, 不是程式邏輯問題, 但記錄下來提醒自己以後改大檔案
    優先用鎖定唯一片段的小範圍替換, 而不是整段回填
- **限制(誠實)**: jsdom 不是真的瀏覽器, 手機上的外觀請實際確認; 反應譜曲線的橫軸目前是均勻
  線性取樣(不是對數), 密集短週期段可能不夠平滑; 沒有 PDF 匯出; `taiwan_code_spectrum` 只是
  標準形狀, 精確查表要接 taiwan-seismic-code-calc 那類工具

#### D4 線性時程(Newmark)✅ 已完成

- 新增 `frame2d/newmark.py`: `newmark_integrate(frame, dt, n_steps, force=None, mass_kind='lumped',
  damping_matrix=None, initial_disp=None, initial_vel=None, beta=0.25, gamma=0.5)` -> `NewmarkResult`
  (`t/u/v/a` 時間序列、`dof_history()`、`member_force_history()`)。預設平均加速度法
  (無條件穩定、無數值阻尼)
- 新增 `frame2d/excitation.py`: `point_force_vector()`(建空間力型態)、`harmonic/step/pulse/ramp`
  (純量時間函式)、`force_series_from_pattern()`(組成 `F(t)=p·f(t)` 給 `newmark_integrate` 用)
- **無質量自由度一定要先靜力凝縮**(跟 `modal.eigen()` 同一個公式 `K_eff=K_dd-K_dm K_mm⁻¹K_md`),
  Newmark 只在有質量的自由度上跑, 跑完用同一組線性關係回填無質量DOF的 u/v/a。**這是修正過的
  設計**: 最早的版本沒有凝縮, 直接把完整系統丟進 Newmark 遞迴——K_hat 的 solve 本身沒問題, 但
  使用者給的非零初始位移如果在無質量DOF(例如集中質量下的轉角)上不是靜力平衡一致的, 從第一步
  就會錯; 在頂端質量懸臂柱的自由振動測試量到「誤差比訊號本身還大 6 倍」才抓到。詳見
  `newmark.py` 模組開頭的完整說明
- 目前限制(刻意, 見模組說明): 不支援力直接施加在無質量DOF上; 阻尼矩陣在無質量DOF上必須整列
  整行是0(D5 的 Rayleigh 阻尼要建立在凝縮後的 K_eff 上, 不是原始 K, 才不會遇到 DAE 複雜度)
- 驗證(`tests/test_newmark.py`, 只依賴 numpy):
  - SDOF自由振動(無阻尼)解析解, Δt 減半時誤差比值 3.94~4.00(平均加速度法的 O(Δt²) 收斂階數)
  - SDOF自由振動(有阻尼)解析解(誤差/振幅 3e-3), 對數遞減率反推 ζ=0.03999(給定0.04)
  - SDOF諧和強迫振動穩態振幅對動力放大係數解析公式(相對差 3.3e-5)、相位對解析公式
  - SDOF無阻尼階躍載重: 峰值=2倍靜位移(相對差2e-3)、峰值時刻=半週期(相對差1e-2)
  - 兩質量懸臂對脈衝載重: 對獨立 Duhamel 積分(數值求積, 不呼叫 frame2d 任何模組)——第一版
    參考解算錯了(誤用地震參與係數公式 ΓᵢᵀM 而不是點力該用的 φᵢᵀp, 差了整整2倍, 用「單位點力
    靜位移」反推才抓到, 過程記在測試檔裡), 修正後相對峰值誤差5e-3
  - 無質量DOF回填的獨立靜力平衡殘差檢查(K@u(t)在無質量DOF上應精確為0), 相對殘差6.3e-14
  - 突變檢查 6 個(Δa公式加回錯誤係數、Δv係數搞錯、凝縮符號錯、有效力漏阻尼項、無質量DOF回填
    符號錯、Khat漏質量項), 全部被抓到
- **對 OpenSeesPy 交叉驗證**(`tests/test_newmark_vs_openseespy.py`, 選用, 沒裝就 SKIPPED):
  SDOF 與 6元素懸臂梁(MDOF)在「從靜止開始、無阻尼、諧和力」情境下跟 OpenSeesPy 的 Newmark
  遞迴一致到機器精度(SDOF 3e-13、MDOF 1e-11)。**範圍限制**(過程中踩到、記錄在測試檔開頭):
  (1) 非零初始位移的自由振動不比——`ops.setNodeDisp(...,'-commit')` 設非零初始位移時, OpenSees
  不會反解一個平衡一致的初始加速度 a0, 就是留在0(用 `ops.nodeAccel()` 直接驗證過), 這跟
  frame2d 的作法(一定解 Ma0=F0-Cv0-Ku0)不同, 不是 frame2d 的 bug; (2) 有阻尼的情況不比——
  兜了一個 zeroLength+Viscous 材料的等效阻尼元素, 結果跟 frame2d 差了 30% 量級, 沒能追出原因,
  改用 test_newmark.py 的解析解驗證阻尼部分; (3) 步階載重不比——`Constant` timeSeries 加步階力
  一開始差 0.5%(Δt減半只降到0.1%, 收斂速率只有O(Δt)), 換成諧和力(`Trig` timeSeries)後立刻對到
  機器精度, 判斷是 OpenSees 對「力在t=0是否已完全作用」的慣例跟 frame2d 不同, 沒有追出 OpenSees
  端確切的機制, 改用 test_newmark.py 層4 的解析解(峰值=2倍靜位移)驗證
- 範例 `examples/newmark_portal_demo.py`: 門型鋼架分別受接近共振的諧和力(看到無阻尼的典型
  「拍音(beating)」振幅包絡)與矩形脈衝力(脈衝結束後在自然週期做無阻尼自由振動), 存位移
  時間歷程圖

#### D5 Rayleigh阻尼 + 地震輸入 ✅ 已完成

- 新增 `frame2d/damping.py`: `rayleigh_coefficients(wi, wj, zeta_i, zeta_j=None)`(反解
  α、β, `zeta_j=None`時退化成經典公式)、`rayleigh_damping_ratio(alpha, beta, omega)`(檢核
  其他模態的實際阻尼比)、`rayleigh_damping_matrix(frame, wi, wj, zeta_i, zeta_j=None,
  mass_kind='lumped')`(組出可以直接傳給 `newmark_integrate(damping_matrix=...)` 的全域矩陣)
- **關鍵設計決定, D4 已經預告過**: `C_dd = α·M_dd + β·K_eff` 建立在**凝縮後**的系統上, 不是
  原始 K。如果用原始 K, 無質量DOF那幾列因為 β·K 這一項不會是0(K在那裡有勁度貢獻, 只有 M 是
  0), newmark_integrate() 的「無質量DOF是純代數約束、不含阻尼」假設就會被打破, 變成
  DAE。`newmark.newmark_integrate` 因此重構出共用的 `condense_for_dynamics()`(`Condensation`
  dataclass), `damping.py` 呼叫同一個函式, 兩邊的凝縮公式永遠是同一個實作, 不會各自算一次
  導致以後改公式時漏改一邊
- `frame2d/excitation.py` 新增: `ground_motion_force(frame, direction, ag, mass_kind)`(等效
  地震力 `P_eff(t)=-M r ag(t)`; 這個力天生滿足「無質量DOF上必須是0」的要求, 因為半正定質量
  矩陣對角線為0則整列/整行為0, 不需要額外檢查)、`absolute_acceleration(result, direction,
  ag)`(重建絕對加速度 `a_abs=a_rel+r·ag(t)`, 給樓層加速度/設備需求用)
- 驗證(`tests/test_damping.py`, 只依賴 numpy):
  - `rayleigh_coefficients()` 解出的 α、β 代回 `rayleigh_damping_ratio()` 精確重現給定的
    ζᵢ、ζⱼ(含 ζᵢ≠ζⱼ 的一般情況), `zeta_j=None` 精確等於經典公式
  - SDOF: `rayleigh_damping_matrix()` 精確等於古典 c=2ζωm(代數恆等式)
  - MDOF(6元素懸臂梁)**模態投影檢核**: 用 `modal.eigen()` 獨立算出的模態形狀, 把 Rayleigh
    矩陣投影到全部 6 個模態(不只是拿來當控制頻率的那兩個), ζₙ=φₙᵀCφₙ/(2ωₙ) 對 `rayleigh_
    damping_ratio()` 的解析公式逐一比對, 相對差 <2e-14, 同時驗證兩個控制頻率之間的阻尼比確實
    低於目標值(Rayleigh 阻尼的已知特性)
  - SDOF對諧和地震輸入的穩態反應對動力放大係數解析公式(相對差 4.8e-5), 跟 D4 的諧和強迫振動
    用同一套公式, 只是這次的等效力是 `-m·ag(t)`
  - 絕對加速度殘差 `m·a_abs+c·v_rel+k·u_rel≈0`(運動方程式本身的恆等式, 獨立於怎麼重建
    a_abs), 相對殘差 1e-9
  - 突變檢查 5 個(Rayleigh矩陣公式係數搞反、C用K而非Keff繞過凝縮、阻尼比公式漏1/2、地震力
    少負號、絕對加速度重建誤用純量廣播), 全部被抓到——最後一個一開始沒抓到(SDOF模型只有
    一個動態DOF, 廣播剛好跟正確答案一樣), 加了一個「轉角DOF不該被地面加速度影響」的獨立檢查
    才抓到
- **對 OpenSeesPy 交叉驗證**(`tests/test_damping_vs_openseespy.py`, 選用, 沒裝就 SKIPPED):
  SDOF、有阻尼的諧和地震輸入, 對 `ops.rayleigh()+ops.pattern('UniformExcitation',...)` 到機器
  精度(5e-15)。**範圍限制**: 只比 SDOF, 不比 MDOF——`ops.rayleigh(alphaM,betaK,...)` 的 βK
  項是對每個元素自己的(未凝縮)勁度矩陣做比例阻尼, 不是對凝縮後的 K_eff, 這兩者在 SDOF(沒有
  無質量DOF)時完全相同, 但在 6 元素懸臂梁(集中質量, 轉角無質量)實測差了 50%(不是數值誤差,
  是方法論本身不同, 正好印證了上面「Rayleigh 阻尼要建立在 K_eff 上」的設計理由)。MDOF 的驗證
  改用 test_damping.py 的獨立模態投影, 不依賴 OpenSees
- 範例 `examples/damping_ground_motion_portal_demo.py`: 門型鋼架 + Rayleigh 阻尼(5%, 兩個
  控制頻率)+ 虛構的衰減正弦波地震脈衝(**不是真實地震紀錄**, 只是示範), 印各模態實際阻尼比、
  存地面輸入/屋頂相對位移/屋頂絕對加速度三張時間歷程圖

#### D6 非線性時程(Newmark + 循環塑鉸)✅ 已完成

- **規劃更正**: 原本寫「先用彈性corotational把時間迴圈本身驗證乾淨」。實際做的時候發現不需要
  這個過渡階段: D7 的雙線性運動硬化塑鉸(`CyclicHingeState`)已經用解析解跟 OpenSeesPy 驗證到
  機器精度, 而且 event-to-event 對分段線性系統是**精確解**, 不需要 Newton-Raphson 收斂——
  直接把它接上 Newmark 比另外做一個「先驗證迴圈本身」的過渡階段更省工、風險更低, 也是這次
  真正做的事: 把 D5 的 Newmark 遞迴改成「力控制」, 每個時間步都是一次 event-to-event 增量分析
  (D7 的 `run_cyclic()` 是位移控制; Newmark 每步的有效力增量是已知的, 反過來是力控制, 結構上
  更簡單), 最大化重用 D7 的元件(`_assemble_stiffness_with_hinges`、`_member_force_increments`、
  `_yield_events`、`_spring_rotation_increment`), 新增的部分只有「怎麼跟 Newmark 時間積分接
  起來」
- 新增 `frame2d/nonlinear_newmark.py`: `nonlinear_newmark_integrate(frame, hinge_states, dt,
  n_steps, force=None, mass_kind='lumped', damping_matrix=None, initial_cum_forces=None,
  beta=0.25, gamma=0.5)` -> `NonlinearNewmarkResult`(`t/u/v/a`、`hinge_M`/`hinge_theta_p`/
  `hinge_work`、`member_force_history`(直接從逐步累積的內力記錄, 對降伏後的桿件也正確,
  跟線性引擎不同不能用 k_local@u 反推)、`events`)
- **一個關鍵簡化(跟D4/D5不同)**: 不對無質量DOF做靜力凝縮——為了避免塑鉸狀態改變時必須
  重新凝縮的複雜度, 直接在完整系統上解, 但這要求動態擾動一律從 u=v=a=0 開始(相對於用
  `initial_cum_forces` 代入的重力預載狀態, 跟 D7 同一個用法), 不接受非零初始位移/速度
  (D4 已經踩過這個坑, 這裡乾脆不開放這個選項)
- 驗證(`tests/test_nonlinear_newmark.py`, 只依賴 numpy):
  - **彈性極限**(SDOF, Mp=∞永不降伏): 對照一個獨立手刻的線性 Newmark 迴圈, 用**同一份**塑鉸
    基礎的近剛接勁度矩陣(不是標準樑元素——兩者本來就不會完全相同, hinge.py 的「未降伏」彈簧用
    RIGID_FACTOR·EI/L=1e8·EI/L 這個很大但有限的數字近似剛接, 不是真的無限剛。第一次比對用了
    錯誤的參照物件, 誤差 1.4e-7 一度以為是bug, 換成同一份近似基準後對到 1.9e-13)
  - **能量平衡**(SDOF, 有阻尼、會降伏): 外力作功 = 動能 + 阻尼耗能 + 桿件內力作功(用桿件局部
    節點力對局部節點位移做功算, 這個定義同時涵蓋彈性儲能與塑性耗能, 不需要另外拆開算, 對任何
    非線性材料模型都成立), 殘差 4.4e-13
  - **單調(準靜態等效)載重下的第一個降伏事件**跟 D7 `run_pushover`(同一個 6 塑鉸門型剛架)
    一致: 位移相對差 5e-5, 力相對差 5e-4。原本想比對整條非線性路徑的最終狀態, 但對「多慢才算
    夠慢」很敏感(降伏後結構軟化, 有效週期變長, 斜坡相對「原始彈性週期」很慢, 仍可能在降伏後
    的路徑上有動態超越/暫時卸載, 最終停留狀態對不上), 改成只比對「第一個降伏事件」(這一段
    路徑上結構還是線性彈性, 單一明確的自然週期, 斜坡夠慢就不會有這個問題, 穩健得多)
  - 突變檢查 6 個(Δa公式、remaining更新比例、θp累積少取絕對值、卸載偵測符號、Khat漏阻尼項、
    member_force_history沒有累積), 全部被抓到; 其中「θp累積少取絕對值」一開始沒被其他層
    抓到(累積淨值仍然對, 只有「絕對值總量」錯), 額外加了一個反覆降伏案例的 θp(絕對值)對
    θp_signed(帶號淨值)比較(有反覆降伏時前者應該明顯更大)才抓到
- **對 OpenSeesPy 交叉驗證**(`tests/test_nonlinear_newmark_vs_openseespy.py`, 選用, 沒裝就
  SKIPPED): SDOF(單一塑鉸懸臂柱, 含降伏與卸載)對到機器精度(7.6e-13)。**範圍限制**: MDOF
  (6塑鉸門型剛架)沒有完全對上——純彈性階段抓到並修正了一個真實問題(OpenSees模型裡內部
  duplicate節點的平移約束設錯: 對支承節點的dup可以直接固定, 但對一般節點的dup只能靠
  equalDOF跟隨、不能額外固定, 修正後彈性階段對到1e-8), 但降伏後仍有成長中的偏差, 沒有找到
  根本原因就沒有繼續追。**這不代表 frame2d 的 MDOF 結果是錯的**——MDOF 規模已經有另外兩種
  跟 OpenSeesPy 無關的獨立驗證撐著(上面提到的能量平衡與 run_pushover 比對), 只是沒有再加上
  OpenSeesPy 這第三種
- 範例 `examples/nonlinear_seismic_portal_demo.py`: 門型鋼架(塑鉸容量刻意設得比彈性需求低)
  受一段虛構的衰減正弦波地震脈衝, 同一個模型分別跑 D5(彈性)跟 D6(非線性)對照, 印尖峰位移、
  降伏次數、殘餘塑性轉角, 存屋頂位移時間歷程(彈性vs非線性疊圖)、整體遲滯迴圈、塑鉸M-θp迴圈
  三張圖
- **限制(誠實)**: 準靜態、小位移(不含 P-Delta/corotational 幾何非線性), 跟 D7 一致; 沒有
  勁度/強度劣化、沒有捏縮; 阻尼矩陣是常數, 不隨塑鉸降伏而改變(業界另一種做法是降伏後改用
  切線剛度比例阻尼, 這裡沒有實作那個選項); 不支援非零初始位移/速度(只能從相對於某個平衡
  狀態的靜止開始)

#### D7 循環塑鉸與準靜態反覆載重 ✅ 已完成(提前)

- **為什麼提前**: 遲滯迴圈來自構件的**組成律**(塑鉸的彎矩-轉角關係)在反覆變形下的行為, 不是
  「動力」本身產生的; 時程分析只是提供一條隨時間變化的位移歷程, 餵進同一個非線性構件, 構件
  畫出來的就是遲滯迴圈。所以先用準靜態反覆載重把構件模型做對、對得上解析解與 OpenSeesPy,
  D6/D8 的非線性時程才有可靠的構件可用。這個步驟只需要靜力求解器, 不必等 D3~D6
- 新增 `frame2d/cyclic.py`(**不修改** `HingeState` 與任何既有求解器):
  - `CyclicHingeState`: 雙線性**運動硬化**塑鉸(彈性→降伏→卸載→反向降伏; 包辛格型, 反向降伏在
    `M = M_max − 2Mp` 而不是 `−Mp`), 等同 OpenSees Steel01(雙線性、`b = R_post/K_e`、不開等向硬化);
    `R_post_yield = 0` 就是彈塑性。狀態: `yielded / direction / alpha(背應力) / theta_p_signed /
    theta_p(∑|Δθp|) / work(累積塑性功)`; `CyclicHingeState.from_hinge_states()` 由既有的
    `initial_hinge_states()` 轉換
  - `run_cyclic(frame, hinge_states, prescribed_dofs, direction, protocol, d_nominal, base_reaction_dofs)`:
    位移控制反覆載重, `protocol` 是控制點的位移目標序列, `make_protocol([幅值...], n_cycles)` 產生標準的
    「幅值遞增、每級 n 圈」歷程。回傳 `CyclicResult`(力-位移、外力功、各鉸的 M / θp / 塑性功歷程、
    降伏與卸載事件)
  - 演算法: 沿用 pushover 的元件與 event-to-event。一段增量內勁度固定 ⇒ 反應對增量線性 ⇒ 降伏事件
    可用線性外插定位, **卸載只會發生在一段的開頭**(符號在段內不變): 用目前塑性勁度試解, 某塑性端的
    彈簧轉角增量與降伏方向相反就切回彈性、重新試解; 一個鉸卸載會改變其他鉸的增量, 所以這個
    符號一致性要疊代到穩定(同接觸問題的 active-set), 不收斂就明確報錯
- 驗證:
  - `tests/test_cyclic.py`: (1) 對獨立的 1D 雙線性運動硬化**回歸映射**參考實作, 每個記錄點比對力-位移
    (等幅多圈、幅值遞增(背應力累積)、不對稱與部分卸載、完全塑性 R=0、不同步長), 最大誤差 5.6e-8(來自
    K_e 有限); (2) 解析頂點: 反向降伏在 `F_a − 2F_y`; (3) 穩態迴圈面積: 解析鞋帶公式 = 數值 ∮F du =
    塑鉸累積塑性功(14.285714); (4) 單調載重與 `run_pushover(HingeState)` 逐點一致(1.4e-14, 6 個塑鉸);
    (5) 塑性狀態下 `M − R·θp` 恆在降伏面上(彎矩記帳與轉角記帳獨立算出來必須一致); (6) 6 塑鉸門型剛架
    反覆載重: 正反向峰值反對稱、穩態迴圈外力功 = 全部塑鉸塑性功; (7) 明確拒絕。**突變檢查 5/5 被抓到**
    (背應力不凍結、不偵測卸載、塑性功用端點值、θp 符號、降伏事件比例)
  - `tests/test_cyclic_vs_openseespy.py`(選用): 對 OpenSeesPy `zeroLength` 旋轉彈簧 + `Steel01`, 4 個案例
    (懸臂柱 2 種歷程、6 塑鉸門型剛架 2 種歷程, 共 54 次降伏/49 次卸載, 各鉸容量與硬化剛度都不同):
    力-位移逐點最大差 **1.5e-13(相對 2e-15)**、轉折點各鉸彎矩差 2e-13
  - 範例 `examples/cyclic_portal_demo.py`: 門型鋼架反覆載重, 印降伏順序、每圈耗能與等效黏性阻尼比,
    存遲滯迴圈圖。塑性崩塌載重 4Mp/h = 100 kN, 模擬平台值 102.5 kN(含硬化), 吻合
- **OpenSeesPy 交叉驗證的實務細節**: Newton 在「近剛接彈簧 + 極小的降伏後剛度」(b = R/K_e ~ 1e-5)下,
  多個塑鉸同時改變狀態會來回震盪、100 次疊代不收斂(K_e ≥ 1000·EI/L 就會; ModifiedNewton、
  NewtonLineSearch、KrylovNewton 都一樣)。比對時兩邊用**同一個**有限的 K_e(把 `hinge.RIGID_FACTOR`
  暫時設成 50), 模型仍然完全相同。另外 OpenSees 的底剪力不能讀 `nodeReaction`(equalDOF 的約束力不在
  該節點上), 要用載重係數
- 限制(誠實): 準靜態、小位移(沒有 P-Delta / 幾何更新); 沒有劣化、捏縮、等向硬化; 只有位移控制;
  網頁尚未接上(D7b)

#### D7b 網頁遲滯迴圈 ✅ 已完成

- **使用方式**: 桿件屬性面板填塑鉸容量 Mp 與降伏後硬化剛度 R(跟 Pushover 一樣)→ 上方「分析」選「**循環**」→
  設定列填「控制節點、方向、位移幅值(逗號分隔)、每級圈數、步長」→ Solve。結果有兩個視圖 +
  一個結果面板: 「遲滯迴圈」(底剪力-控制點位移, 紅點 = 塑鉸降伏)、「塑鉸 M-θp」(各鉸的彎矩-累積塑性轉角小圖
  網格, 依降伏先後排序)、面板裡的「每級幅值一個完整迴圈的耗能與等效黏性阻尼比」與「塑鉸降伏順序」。
  單位跟「單位設定」一致(位移、力、彎矩; 耗能用彎矩單位); 模型裡的載重當重力預載
- 核心新增(`frame2d/cyclic.py`): `loop_summary()`(每個幅值的穩態迴圈耗能與 ξ_eq; 只走 1 圈時在載入到下一級的途中
  內插閉合, 沒閉合就回 None 不亂給數字)、`cyclic_analysis()`(網頁一站式入口, 兩個後端共用: 建塑鉸狀態、控制點與底反力
  自由度、重力預載、protocol、步數上限 40000 快速失敗)、`cyclic_to_dict()`(JSON, 只附曾降伏的鉸)。
  `run_cyclic` 新增: 起始狀態(重力預載)已超過 Mp 明確報錯; 矩陣奇異時若有 R=0 的塑鉸, 訊息指出原因與建議
- **實測發現**: 門型剛架 6 個鉸全部 R=0 時, 兩個 R=0 的鉸夾住的節點轉角沒有勁度, 矩陣奇異(單一鉸的構件沒問題)。
  這是物理上的不定, 不是 bug; 建議填小的正值(0.01·EI/L 時平台力 100.15 kN, 對照崩塌載重 4Mp/h = 100 kN)。
  網頁說明文字與錯誤訊息都已寫明
- 後端: FastAPI 與 stdlib 各一個 `POST /cyclic`(`FrameIn` 新增 cyclic_control_nodes / cyclic_weights /
  cyclic_direction / cyclic_amplitudes / cyclic_n_cycles / cyclic_step); 錯誤都是清楚的 400
- 驗證:
  - `tests/test_cyclic.py` 層8: `loop_summary` 對解析迴圈面積(6 個幅值等級, 含 1 圈內插閉合與未閉合)、ξ_eq 公式、
    重力預載(降伏底剪力 = (Mp ∓ M0)/L, 鉸在降伏事件時彎矩恰為 ±Mp)、預載超過 Mp、11 種錯誤輸入
  - `tests/test_web_cyclic_api.py`: stdlib 後端對核心逐項相同; FastAPI 與 stdlib 回傳完全相同; 錯誤訊息; 重力預載
  - `tests/test_web_cyclic_e2e.py`(選用, 需要 node + jsdom): 用 jsdom 載入真正由後端提供的 index.html, 走「匯入
    JSON」載入 `examples/portal_cyclic_demo.json` → 選「循環」→ 輸入驗證(全空、負幅值、步數過多、節點格式)→ Solve
    → 網頁結果與核心逐項一致(642 點、耗能 0 / 1902.7 / 6134.6 / 12035.5 J、最大底剪力 102460.7 N、6 個鉸降伏順序)→
    迴圈圖(1 條折線、58 個降伏紅點)、M-θp 圖(6 張小圖)、結果表 → 單位切換(15 mm → 0.015 m, 55.818 kN →
    55818 N)→ 拿掉塑鉸容量得到清楚訊息且不破壞舊結果 → 重置。突變檢查 5/5 被抓到(幅值單位轉換錯、紅點沒畫、
    換單位不重畫、浮點雜訊顯示成科學記號、塑鉸視圖畫成整體迴圈)
  - 另外把頁面實際的 SVG 轉成圖片肉眼檢查過版面(軸、刻度、標題、紅點、小圖網格)
- **限制(誠實)**: jsdom 不是真的瀏覽器——沒有版面/CSS 渲染與真實觸控, 手機上的外觀請實際確認; 目前沒有動畫、
  沒有 PDF 匯出; 準靜態、小位移、無劣化與捏縮(見 D7)

#### D8 非線性地震反應一站式入口 ✅ 已完成

- 新增 `frame2d/seismic.py`: `seismic_analysis(frame, ag, direction='x', dt=None, n_steps=None,
  zeta=0.05, damping_modes=(0,2), mass_kind='lumped', apply_gravity_loads=True)` ->
  `SeismicResult`。**這個檔案不新增任何求解邏輯**, 純粹是把 D2(模態)/D5(Rayleigh阻尼+地震
  輸入)/D6(非線性時程)/D7(重力預載沿用 `pushover.apply_gravity()`)串成一次呼叫, 跟
  `cyclic.cyclic_analysis()`(D7 給網頁用的一站式入口)是同一種角色
- `SeismicResult` 提供的便利方法: `peak_displacement(node, direction)`、
  `peak_drift_ratio(node_top, node_bottom, height, direction)`(層間位移角)、
  `absolute_acceleration()`(轉發 D5 的 `excitation.absolute_acceleration()`)、
  `energy_balance()`(算出 KE/阻尼耗能/桿件內力作功/外力作功四條時間序列)
- **能量平衡**: 外力作功(等效地震力對相對位移做的功)= 動能 + 阻尼耗能 + 桿件內力作功(用
  桿件局部節點力對局部節點位移做功算, 同時涵蓋彈性儲能與塑性耗能), 這條式子在地震輸入下
  跟一般外力沒有本質差異, 一樣精確成立(見 `seismic.py` 模組說明); 有重力預載時驗證過預載
  狀態不影響這條平衡式(重力在純相對動態運動下不作淨功)
- 驗證(`tests/test_seismic.py`, 只依賴 numpy): 能量平衡(無/有重力預載, 殘差 1e-10~1e-12
  量級)、`seismic_analysis()` 對照手動組裝 D2/D5/D6 逐項相同(位移完全相同, 加速度相對差
  <1e-9)、Rayleigh 阻尼比對解析公式(含全部模態的模態投影檢核)、`peak_displacement`/
  `peak_drift_ratio` 對照直接陣列索引、明確拒絕。突變檢查 5 個, 4 個被抓到, 1 個
  (Rayleigh控制模態索引wi/wj寫反)在ζᵢ=ζⱼ時數學上就是等價的(不是bug, 兩個模態目標阻尼比
  相同時交換順序對解無影響, 這個「突變」本來就該存活)
  - **這一層抓到一個真實 bug**: `energy_balance()` 第一版把每根桿件的功增量累加後忘記對時間
    做 `np.cumsum()`(Wext/Wdamp 都有 cumsum, Wint 少了), 導致殘差高達外力功尺度的 82%。
    追查時先用單一構件(D6測試的SDOF)重跑一次確認底層公式沒問題, 才發現是這個檔案自己的
    accumulate邏輯漏了一步——這也是為什麼即使底層元件都個別驗證過, 整合層本身還是需要
    自己的測試
- **沒有另外做 OpenSeesPy 交叉驗證**: `seismic_analysis()` 內部呼叫的
  `nonlinear_newmark_integrate()` 已經在 D6 對 OpenSeesPy 驗證過, `ground_motion_force()`/
  `rayleigh_damping_matrix()` 已經在 D5 驗證過; 這裡新增的只是「怎麼把它們接起來」, 用「手動
  組裝結果逐項相同」這個測試(上面提到的那個)驗證組裝本身沒接錯, 比另外重跑一次 OpenSeesPy
  更對症(能真正測到「這個檔案自己寫的代碼」, 而不是重複驗證已經驗證過的底層數學)
- 範例 `examples/seismic_analysis_demo.py`: 跟 D6 demo 同一個模型跟同一段地震歷程, 改用一行
  `seismic_analysis()` 呼叫, 畫出屋頂相對位移+絕對加速度雙軸圖、能量平衡堆疊圖(動能/阻尼
  耗能/桿件內力作功疊加, 對照外力作功的虛線, 視覺確認平衡)
- **限制(誠實)**: 繼承 D6 的所有限制(準靜態小位移、無勁度/強度劣化、常數阻尼矩陣、動態
  擾動只能從相對靜止開始); `energy_balance()` 的 `Wext`/`Wint` 個別不保證單調(可能因為彈性
  能量的可逆振盪而暫時下降), 只有它們的**組合**保證平衡, 用時間上的最大值(而不是終值)當
  比例尺度時要留意這一點

#### D6b 網頁非線性地震反應 ✅ 已完成

- **使用方式**: 桿件填塑鉸容量 Mp(跟Pushover/循環一樣)、斷面填密度 ρ 或節點填質量 → 上方
  「分析」選「**非線性地震**」→ 設定控制節點、方向、阻尼比、Rayleigh控制模態(1起算, 例如
  1,3)、質量矩陣、是否先做重力預載、Δt/總步數(留空=自動) → 選地震歷程來源: 「衰減脈衝」
  (填峰值加速度/頻率/衰減係數, **不是真實地震紀錄**, 示範用)或「自訂」(週期-加速度資料表,
  單位g, 可增減列)→ Solve。結果有 4 個視圖: **地震動畫**(結構隨時間變形, 有播放/暫停/滑桿,
  直線畫法不是精確Hermite曲線——動畫幀多時逐幀重算曲線太貴, 跟D2b模態視圖同一種取捨)、
  **地震時程**(地面加速度+控制節點位移雙圖, 紅點標目前播放位置)、**遲滯迴圈**(控制節點位移
  vs「有效慣性力合力」, 用完整時間序列畫、不是動畫抽稀後的幀, 紅點標目前位置)、**塑鉸
  M-θp**(降伏過的塑鉸各一張小圖)。單位跟「單位設定」一致
- 核心新增(`frame2d/seismic.py`): `sine_pulse_ground_motion(amplitude_g, freq_hz, decay)`
  (簡化衰減正弦波)、`custom_ground_motion(points_g)`(分段線性內插, 自動排序)、
  `nonlinear_seismic_web_analysis()`(網頁一站式入口: 建地震歷程 callable → `seismic_analysis()`
  → 動畫影格抽稀, 均勻取樣到 `max_frames`, 首尾保留, 不超過就不抽稀)、`seismic_to_dict()`
  (JSON; 塑鉸M-θp/時程/遲滯迴圈用完整未抽稀的時間序列, 只有「每個時間步存全部節點位移」這件
  事——動畫要用的——才抽稀, 因為那個資料量遠大於其他曲線)
- 後端: FastAPI 與 stdlib 各一個 `POST /nonlinear_seismic`(`FrameIn` 新增 seismic_control_node/
  seismic_direction/seismic_zeta/seismic_damping_modes/seismic_mass_kind/
  seismic_apply_gravity_loads/seismic_dt/seismic_n_steps/seismic_ground_motion_type/
  seismic_pulse_amplitude_g/seismic_pulse_freq_hz/seismic_pulse_decay/
  seismic_custom_points_g); 步數上限 4000(沿用 `nonlinear_seismic_web_analysis` 的
  `MAX_WEB_STEPS`), 超過時明確錯誤訊息建議加大Δt
- 驗證:
  - `tests/test_seismic_web.py`: `sine_pulse_ground_motion`/`custom_ground_motion` 公式逐點
    比對(含自動排序、邊界外夾住)、`nonlinear_seismic_web_analysis`+`seismic_to_dict` 跟直接
    呼叫 `seismic_analysis()` 逐項相同(含JSON化後能量平衡仍然成立)、動畫抽稀正確(首尾保留、
    小案例不抽稀)、明確拒絕。突變檢查 5 個, 一開始只抓到 3 個——「有效慣性力合力算錯(用M而非
    M@r)」跟「hinges降伏次數篩選條件顛倒」兩個測試沒抓到, 各自補上獨立算出的期望值比對
    (前者對照獨立算的「總質量×地面加速度」、後者直接檢查hinges清單內容跟n_yield>0)才抓到
  - `tests/test_web_seismic_api.py`: stdlib 後端對核心逐項相同(脈衝與自訂地震歷程都測);
    FastAPI 與 stdlib 回傳完全相同; 各種錯誤訊息
  - `tests/test_web_seismic_e2e.py`(選用, 需要 node + jsdom): 匯入範例模型 → 選「非線性地震」
    → 輸入驗證 → Solve → 網頁結果與核心逐項一致(含能量平衡最終值)→ 地震動畫(播放/暫停/
    滑桿都測過, 含**變形放大倍率對照獨立算出的期望值**——這個檢查是補測出來的, 見下面的
    踩坑記錄)→ 地震時程圖(雙圖+紅點)→ 遲滯迴圈(確認用完整時間序列、不是抽稀後的動畫幀)→
    塑鉸M-θp小圖網格 → 結果表 → 單位切換 → 自訂地震歷程編輯 → 拿掉塑鉸容量得到清楚訊息且
    不破壞舊結果 → 重置。突變檢查 5/5 被抓到(播放不前進、滑桿拖動不更新影格、遲滯迴圈誤用
    抽稀後的資料、重置沒清乾淨、變形放大倍率算錯)
  - **踩坑記錄**: 第一版的 e2e 測試沒有檢查「變形放大倍率」的實際數值, 只檢查畫面上有沒有畫出
    線、線的端點座標有沒有 NaN——結果一個「忘記除以 maxDisp」的突變(倍率變成只有正確值的
    1/45, 動畫看起來幾乎不會動, 但仍然是個有效數字, 線也確實畫出來了)沒被抓到。加上「獨立
    重算倍率該有的值, 逐一比對」之後才抓到。這是這個檔案在 D3b/D7b 之後第三次遇到「畫面上有
    東西、座標也不是NaN, 但數值本身錯」這種只有比對實際數字才抓得到的錯誤類型
  - 另外把頁面實際的 SVG 轉成圖片肉眼檢查過全部 4 個視圖(動畫變形形狀、時程雙圖+紅點、
    衰減螺旋狀遲滯迴圈+紅點、塑鉸M-θp), 版面跟數值都正常
- **限制(誠實)**: jsdom 不是真的瀏覽器, 手機上的外觀請實際確認; 動畫用直線畫變形形狀
  (不是精確的Hermite三次曲線); 沒有動畫匯出(GIF/影片)、沒有PDF匯出; 繼承 D6/D8 的所有限制
  (準靜態小位移、無勁度/強度劣化、常數阻尼矩陣); 「有效慣性力合力」不是嚴格的基底反力
  (nonlinear_newmark沒有算反力), 前端文案已經明確標註這件事, 不叫「基底剪力」

#### D9 真實地震紀錄支援(PEER NGA) ✅ 已完成

- 新增 `frame2d/ground_motion_io.py`: `parse_peer_nga(text)` 解析 PEER NGA(太平洋地震工程
  研究中心強震資料庫)的 .AT2 文字格式, 回傳 (dt, 加速度陣列[g])。**這個檔案只負責讀檔解析,
  不含任何新的求解邏輯**——解析結果接上 D5 已經驗證過的 `custom_ground_motion()`, 就能餵進
  D6/D8/D6b 整條非線性地震反應管線
- **格式容錯設計**: 不同資料庫/工具匯出的 .AT2 檔案在「每行幾筆數值」「有沒有逗號分隔」
  「NPTS/DT那行結尾有沒有"SEC"字樣」這些地方常有差異。解析器不假設固定的每行筆數或欄寬,
  只用正規表示式找出 NPTS/DT 那一行, 之後把剩下所有文字直接依空白切開轉成浮點數, 取前 NPTS
  個——對這類格式差異有容錯能力, 同一份測試套件驗證過 4 種排版變體都能正確解析
  (`tests/test_ground_motion_io.py`)
- `peer_nga_to_points(text, max_points=None)`: 轉成 `custom_ground_motion()` 要的
  `[(t, ag_g), ...]` 點列表, `max_points` 選用抽稀(真實強震紀錄常有上萬個取樣點, 網頁
  傳輸/內插不必要地肥大時用, 均勻覆蓋、頭尾保留)
- `seismic.nonlinear_seismic_web_analysis()` 新增 `ground_motion_type='peer_nga'`
  (`peer_nga_text`=檔案文字內容, `peer_nga_max_points`=抽稀上限預設2000), 網頁後端
  `/nonlinear_seismic` 新增對應欄位(`seismic_peer_nga_text`/`seismic_peer_nga_max_points`)
- 前端: D6b 的「非線性地震」設定列, 地震歷程來源多一個「PEER NGA檔」選項, 用瀏覽器檔案選取
  (跟既有的「匯入JSON」同一種 FileReader 讀檔模式)讀取 .AT2 檔案文字內容, 送給後端解析
  (不在前端重新實作解析邏輯, 維持解析只有 Python 這一份實作)
- 驗證:
  - `tests/test_ground_motion_io.py`(只依賴 numpy): 標準格式逐點比對; 4種格式變體(每行
    8筆/沒有逗號/沒有SEC字樣/每行1筆)都能正確解析; `peer_nga_to_points()` 的 t=i·dt 逐點
    驗證、抽稀頭尾保留、原始點數不夠時不抽稀; 明確拒絕(找不到NPTS/DT、NPTS或DT不是正數、
    數值筆數不夠)。突變檢查 4 個, 一開始只抓到 3 個——「取數值時多取1個」這個 off-by-one
    錯誤在乾淨的測試資料(檔案剛好只有NPTS個數值, 沒有多的可以誤取)下完全沒有影響, 加一個
    「檔案末尾多一個數值」的案例才抓到
  - `tests/test_seismic_web.py`/`tests/test_web_seismic_api.py`: `nonlinear_seismic_web_
    analysis(ground_motion_type='peer_nga')` 對照「獨立呼叫 `custom_ground_motion(
    peer_nga_to_points(...))` 再接 `seismic_analysis()`」的路徑逐項相同; stdlib/FastAPI
    兩後端一致; 明確拒絕(缺文字、格式錯誤)
  - `tests/web/seismic_e2e.js`(選用, 需要 node + jsdom): 上傳一個合成的 .AT2 檔案(用
    `File`+`FileReader`, 跟「匯入JSON」同一套 jsdom 模擬方式), 確認讀到的內容跟原始文字
    完全一致、送出分析成功。突變檢查 3 個(payload漏送文字、面板切換沒隱藏、沒選檔案時
    驗證被拿掉)全部被抓到
  - 另外把網頁實際跑出來的地震時程圖轉成圖片肉眼確認: 上傳的合成紀錄(8秒長)確實驅動了
    結構反應, 分析時長(預設約10秒, 比紀錄本身長)超出紀錄範圍後正確夾在最後一個值(不是
    外插到0或報錯), 這是 `custom_ground_motion()` 文件裡已經講清楚的行為, 不是bug
- **限制(誠實)**: 只支援 PEER NGA .AT2 格式(業界最常見的公開格式之一, 但不是唯一格式,
  例如日本 K-NET/KiK-net、歐洲 ESM 資料庫用不同格式, 目前不支援); 假設檔案裡的加速度單位
  是 g(PEER NGA 慣例上就是, 但不會主動驗證檔案內容是否真的符合這個慣例); 沒有自動判斷
  南北/東西分量或多分量檔案的邏輯(.AT2 檔案本身通常已經是單一分量, 這符合這裡的雙向水平
  分析假設)

#### D10 網頁匯出 Markdown: 模態/反應譜/循環/非線性地震 ✅ 已完成(Markdown; PDF 規劃中)

- 之前「匯出Markdown」「匯出PDF」兩個按鈕只認 `lastResult`(只有線性/P-Delta/Pushover 會設定這
  個變數), 模態/反應譜/循環/非線性地震四種分析結果完全沒有匯出功能, 按鈕點了沒反應(顯示
  「請先按Solve才有結果可匯出」)。這次補上這四種的 Markdown 匯出
- **設計原則(對應使用者的需求: 拿到報告的人要能自己重建模型、驗證所有結果)**: 每份報告除了
  結果, 一定包含完整的模型輸入(節點/斷面/桿件含塑鉸容量/支承/載重, 用共用函式
  `buildInputDataMd()`)跟質量設定(`buildMassDataMd()`, D2以後的分析都建立在質量矩陣上,
  沒有這段還原不出動力分析要用的質量), 再加上這次分析實際用的設定參數(阻尼比、地震歷程
  定義等)。**這些內容一律從 `lastXxxPayload`(每次求解實際送給後端的那份 payload, 求解成功
  後存下來)讀, 不是從匯出當下的即時 `model`/DOM 狀態讀**——這樣報告永遠精確對應「畫面上這組
  結果」是拿哪個模型、哪組設定算出來的, 不會因為匯出前又動了模型或設定列而對不上
- 地震歷程來源是 PEER NGA 上傳檔案時, 把**完整的原始 .AT2 檔案內容**寫進報告的程式碼區塊裡,
  可以直接存成檔案重新上傳, 完全重現同一段地震歷程
- 循環/非線性地震還多附了**完整的逐點數值表**(力-位移曲線、時間歷程), 不是只有摘要統計——
  拿到報告的人可以逐點核對畫面上看到的曲線, 不用相信任何摘要數字
- 驗證: 四個 e2e 測試(`tests/web/modal_e2e.js`/`rsa_e2e.js`/`cyclic_e2e.js`/`seismic_e2e.js`)
  各自加了匯出檢查, 確認標題/輸入資料/分析設定/結果表逐項正確(包含跟 `modalResult`/
  `rsaResult`/`cyclicResult`/`seismicResult` 的實際數值比對, 不是只檢查「有沒有東西」)。
  突變檢查 5 個, 一開始只抓到 4 個——「模態表 Γx 欄位被誤填成 Γy 的值」這個錯誤因為只檢查
  了列數沒檢查數值而沒被抓到, 加了一個對照 `modalResult.modes[0].gamma_x` 實際數值的檢查
  才抓到
- **規劃中**: PDF 匯出(需要重新畫圖, 工作量比 Markdown 大很多, 這次先不做)

#### D11 網頁匯出 PDF: 模態 ✅ 已完成(反應譜/循環/非線性地震規劃中)

- 新增 `webapi/pdf_export.py` 的 `build_modal_pdf_report()`: 週期/頻率/參與係數/有效質量比表
  + 振型圖(每頁最多6個模態, **直接重用 `modal_to_dict()` 已經算好的變形曲線座標, 不重算**
  ——`_mode_shape_fig()` 對這份資料是精確的 pass-through)+ 質量設定頁(新增
  `build_mass_data_page()`, D2以後的分析都建立在質量矩陣上, 沒這頁還原不出質量設定)+
  完整輸入資料頁(重用既有的 `build_input_data_pages()`)
- `/export/pdf` 依 `analysis_type` 分派(FastAPI 與 stdlib 兩後端都接): `modal` 重新跑一次
  `eigen()` 產生報告; `rsa`/`cyclic`/`seismic` 目前回傳清楚的「還在做, 目前只能匯出
  Markdown」400 訊息, 不是壞掉的PDF或不明錯誤; `FrameIn.analysis_type` 的允許值從
  `linear`/`pdelta`/`pushover` 擴充到含 `modal`/`rsa`/`cyclic`/`seismic`(只有
  `/export/pdf` 認得後四種, `/solve` 還是只認前三種)
- 驗證(`tests/test_pdf_export_modal.py`): **不依賴pymupdf等PDF解析套件**——直接檢查
  `_mode_shape_fig()`/`build_mass_data_page()` 回傳的 matplotlib Figure 物件本身(座標軸
  標題文字、線段資料點、表格儲存格文字), 這些內容之後會被 `PdfPages` 原封不動存成PDF頁面,
  檢查Figure本身等於檢查了最終PDF長什麼樣子。分頁邏輯用monkeypatch監聽
  `build_modal_pdf_report()`內部實際呼叫`_mode_shape_fig()`時傳入的模態編號區間, 不是隔著
  PDF猜頁數。突變檢查4個, 一開始只抓到1個——「振型圖標題只檢查週期/頻率沒檢查參與係數」
  「分頁大小改變但總頁數不變」「質量表某一欄漏轉換單位但檢查用的是全表模糊比對」這三種都
  是「檢查得不夠精確」型態的測試盲點, 逐一補上精確比對後4個全部抓到。第4層(FastAPI跟
  stdlib文字內容比對)選用, 需要pymupdf(不是這個專案原有的依賴, 沒裝就SKIPPED)
  - **這裡也踩到一個測試撰寫本身的坑**: 第一版把「選用pymupdf」跟後面「還需要fastapi」的
    程式碼混在同一個módule層級的try/except裡, 但只catch了ImportError一次——結果在沒有
    fastapi、但有pymupdf的環境下(這個專案的「Termux風格」測試環境正是這樣), `import
    pymupdf`成功進到`else`分支, 接著`from fastapi.testclient import TestClient`才真正
    炸開, 而且是在pytest**收集(collection)** 階段就炸開(這種純腳本風格的測試檔案沒有
    `if __name__=="__main__"`包住主要邏輯時, pytest匯入檔案就等於直接執行全部內容),
    導致整個測試套件在Termux風格環境下完全跑不起來、不是只有這個檔案的問題。修法是把
    pymupdf 跟 fastapi 的 import 合併在同一個 try 區塊裡, 讓兩者只要缺任何一個都會被同一個
    except ImportError 擋下來
- **對 PDF 匯出補的測試涵蓋率缺口(誠實記錄)**: 這個專案原本的線性/P-Delta/Pushover PDF
  匯出(`build_pdf_report()`/`build_pushover_pdf_report()`)完全沒有自動化測試, 這次沒有
  回頭補——只針對新增的模態PDF匯出建立測試, 沒有擴大範圍去補舊功能的測試債
- 範例: 網頁選「模態」分析、Solve、按「匯出PDF」

