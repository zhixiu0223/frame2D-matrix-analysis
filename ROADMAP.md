# frame2d 開發規劃

這份文件記錄「接下來要做什麼、為什麼」。README.md 專心回答「現在能幹嘛、
怎麼用」,兩份文件分工。

## 目前狀態 (已完成)

```
frame2d
├── Model:    Node / Section / Member(frame/truss/cable) / Support(bool)
├── Elements: 樑柱(彎曲+軸向) / 桁架(僅軸向,可拉可壓) / 纜線(僅軸向,只受拉+自動鬆弛迭代)
├── Loads:    節點集中力/力矩(point_load) / 桿件全長均佈載重(distributed_load)
├── Solver:   線性靜力, 直接勁度法
├── Postprocess: N/V/M圖、變形圖、反力
└── Validation: 解析解 / sd_framework+anastruct(11案例46數值) / SW FEA(third-party) /
                桁架節點法手算 / 懸臂梁斜張橋(位移法vs力量法交叉驗證)
```

這已經是一個能處理「frame + truss + cable 混合結構」的通用矩陣位移法核心,
不再是單純的驗證腳本集合。**接下來不再繼續擴張「結構種類」**(不做第二種
斜張橋案例、不做更多bridge專屬功能),把力氣收回來讓「載重/邊界條件/元素/
後處理」這四個核心抽象更完整、更穩定。

## 下一步: Load System v2

SW FEA 的載重介面對照 frame2d 現況:

| SW FEA 類型 | frame2d 現況 |
|---|---|
| Nodal Point Load | ✅ `point_load(node, fx, fy, m)` |
| Distributed(全長, 垂直桿件) | ✅ `distributed_load(member, w_start, w_end)` |
| Point Load(桿件中間任意位置) | ❌ 未做 |
| Moment(桿件中間任意位置) | ❌ 未做 |
| Distributed(局部段, 不用整根桿件) | ❌ 未做 |
| Support Displacement(強制位移/沉陷) | ❌ 未做 |
| (角度/局部-全域座標系選項) | ❌ 未做, 暫緩(見上方說明) |

**優先順序:**

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

### Phase 2: 局部段均佈載重
`distributed_load` 的 `DistributedLoad` dataclass 這次直接加上 `x_start`,
`x_end` 兩個欄位(預設None=整根桿件, 向下相容現有呼叫方式), 是Phase 1
點載重公式沿桿長積分的推廣, 做完Phase 1後這項不難。
驗證: 簡支梁局部段均佈載重, 解析解比對。

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
反力公式沒處理(等同truss但沒有內部載重公式)。

## 暫緩/刻意不做的項目

- 均佈載重的任意角度/座標系選項: 沒有具體驗證場景前不加欄位
- 斜張橋以外的其他特殊橋型/結構型式: frame+truss+cable已經證明夠通用,
  不需要再做專屬功能證明
- 非線性(大變形/材料非線性): 目前定位是線性靜力分析核心, 不在範疇內
