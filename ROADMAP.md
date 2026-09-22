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
| D4 | 線性時程 (Newmark) + 諧和/脈衝/任意力輸入 | `newmark.py`, `excitation.py` | ⬜ |
| D5 | Rayleigh 阻尼 + 地面加速度輸入 | `damping.py` | ⬜ |
| D6 | 非線性時程框架 (Newmark + Newton, 先用彈性 corotational) | `transient.py` | ⬜ |
| D7 | 循環塑鉸(遲滯) + 準靜態反覆載重(**提前**, 見下) | `cyclic.py` | ✅ |
| D7b | 網頁: 「循環」分析、遲滯迴圈圖、塑鉸 M-θp 圖、每級耗能與等效阻尼比 | webapi | ✅ |
| D8 | 非線性地震時程 + 能量平衡 | — | ⬜ |

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

#### D4 線性時程

- Newmark-β, 預設平均加速度法 (β=1/4, γ=1/2, 無條件穩定、無數值阻尼),
  有效勁度 `K̂ = K + γ/(βΔt) C + 1/(βΔt²) M`
- 輸入: 諧和力、脈衝/階躍、任意時間序列力
- 驗證: 單自由度自由振動(阻尼頻率、對數衰減率)、諧和穩態振幅對照動力放大係數
  `1/√((1−r²)² + (2ζr)²)`、無阻尼階躍載重峰值為靜位移的 2 倍
- 內部交叉驗證: 同一模型「Newmark 直接積分」對照「D2 模態疊加」
- Δt 收斂測試, 記錄週期伸長誤差隨 `Δt/T` 的變化

#### D5 阻尼與地震輸入

- Rayleigh 阻尼 `C = αM + βK`; 兩個模態取相同阻尼比 ζ 時
  `α = 2ζω_iω_j / (ω_i+ω_j)`, `β = 2ζ / (ω_i+ω_j)`; 輸出其他模態的實際阻尼比
  當檢核(Rayleigh 在兩個控制頻率之外會偏離)
- 地面加速度輸入: `P_eff(t) = −M r a_g(t)`; 輸出同時提供相對位移與絕對加速度,
  正負號慣例寫清楚
- 驗證: 單自由度對地面加速度的反應對照 Duhamel 積分; 對很多個週期跑單自由度時程
  取峰值, 得到的反應譜就是 D3 的基準

#### D6 非線性時程框架

- Newmark + Newton-Raphson, 殘餘力
  `R = P_eff − M a − C v − f_int(u)`, 切線勁度 `K̂_T`
- 狀態機制: trial / commit, 迭代不收斂時能還原; 重力預載重用 `newton.py` 既有作法
- **第一版材料為彈性**, 只開 corotational 幾何非線性, 用來把時間迴圈本身驗證乾淨
- 驗證: 小振幅時必須重現 D4 線性結果(振幅趨近 0 時相對誤差趨近 0); 大振幅懸臂
  自由振動的週期伸長對照 OpenSeesPy corotational
- 不收斂時如實回報(沿用 pushover 一貫作法), 可選擇性支援時間步長細分

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

#### D8 非線性地震時程

- 整合 D5~D7, 加上能量平衡檢核: 輸入能 = 動能 + 阻尼耗能 + 應變能 + 遲滯耗能
- pushover 與時程對照: 容量曲線 vs 時程峰值位移與塑鉸形成順序; P-Δ 效應開/關比較
- 對 OpenSeesPy 完整模型(多層框架)比對, 可延伸現有 pushover 的比對案例

### 對 OpenSeesPy 交叉驗證的注意事項

- **質量定義要一致。** `elasticBeamColumn` 有 `-mass`(集中)與 `-cMass`(一致)選項
  (以你安裝版本的文件為準); 轉動慣量是否計入要兩邊明確指定, 否則頻率對不上
- **特徵向量比 MAC**, 不比數值; 實測特徵值相對誤差 3e-13、MAC = 1.0, 門檻設 1e-9
- **OpenSeesPy 的 `ops.eigen(n)`(預設 ARPACK)在有無質量DOF時會失敗**(集中質量下的
  轉角自由度): `ArpackSolver::Error ... Could not build an Arnoldi factorization`。
  要用 `ops.eigen('-fullGenLapack', n)`(很慢, 小模型夠用)
- **OpenSees `modalProperties` 的參與係數**: 集中質量與 frame2d 逐項一致; 一致質量定義不同,
  只比特徵值與 MAC
- **阻尼要用同一組 (α, β)**, 不要各自從阻尼比反算
- **時程要同 Δt、同 Newmark 參數、同輸入序列**; 注意 OpenSees `timeSeries` 的
  因子與正負號慣例
- **單位先統一**再比: 兩邊用同一套一致單位(例如 N-m-kg-s, 跟網頁後端一致), 質量單位 = 力單位·s²/長度單位

### 分支與版本

- 舊標記不動: `v1.0-linear-elastic`(線彈性)、`v1.1-nonlinear-static`(非線性靜力)
- 動力分析在 `feature/dynamics` 進行, 每個 Stage 至少一個 commit; D0 的回歸測試
  通過後才開始 D1
- 每個 Stage 完成後可選擇打 tag(例如 D2 完成打 `v1.2-modal`), 並在此表把 ⬜ 改成 ✅,
  補上實作記錄與遇到的 bug, 沿用上面 Phase 1~4b 的記錄風格
