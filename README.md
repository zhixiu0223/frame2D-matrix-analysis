# frame2d — 2D 矩陣位移法(直接勁度法)通用框架分析

從 [slope_deflection_framework](https://github.com/zhixiu0223/slope_deflection_framework)
分出的新專案:傾角變位法 repo 保留「每一步推導透明」的教學價值,這裡專注「任意拓樸都能解」。

下一步規劃見 [ROADMAP.md](ROADMAP.md)。

## 目前範疇

- 2D Euler-Bernoulli 樑柱元素 (frame),每節點3自由度 (ux, uy, rot)
- 桁架元素 (truss):兩端鉸接、僅軸向勁度,可拉可壓
- 纜線元素 (cable):跟truss共用軸向公式,但只受拉,受壓時自動判定鬆弛、
  移除勁度貢獻、重新求解,反覆迭代至收斂(見下方「纜線元素」章節)
- **尚未支援**: 局部段均佈載重(distributed_load目前只能整根桿件)、支承強制位移(沉陷)、
  內部鉸接(internal hinge)、幾何非線性(P-Delta)、材料非線性——這些等基本
  框架穩定後再視需求加入,不要為了還沒出現的需求先付架構成本(詳細優先順序
  見ROADMAP.md)

## 結構

```
frame2d/
  model.py    — Node / Section / Member / Support / PointLoad / DistributedLoad / Frame2D (建模API)
  elements.py — 局部6x6勁度矩陣、座標轉換矩陣、均佈載重固定端反力公式
  solve.py    — 組裝、邊界條件(partition method)、求解、桿端內力回代
tests/
  test_cantilever.py            — 懸臂梁點載重 vs 解析解 (機器精度)
  test_simply_supported_udl.py  — 簡支梁均佈載重 vs 解析解 (機器精度, 驗證分佈載重公式)
  test_case08_vs_swfea.py       — 兩層兩跨鋼架 vs SW FEA第三方工具報告 (rel_err<0.2%)
```

## 使用範例

```python
from frame2d import Frame2D, solve
from frame2d.plotting import plot_all
import matplotlib.pyplot as plt

f = Frame2D()
f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
f.add_section('sec', E=200e9, I=8e-5, A=1e-2)
f.add_member(0, node_i=0, node_j=1, section='sec')
f.add_member(1, node_i=1, node_j=2, section='sec')
f.add_member(2, node_i=3, node_j=2, section='sec')
f.fix(0).fix(3)
f.point_load(1, fx=12.0)

result = solve(f)
print(result.reactions)
print(result.member_results[0].end_forces_local)  # [Fx1,Fy1,M1,Fx2,Fy2,M2]

fig = plot_all(f, result)   # 六合一: 結構圖/受力圖/變形圖/軸力圖/剪力圖/彎矩圖
fig.savefig('output.png')   # Pydroid3上可改用 plt.show() 直接跳出檢視畫面
```

更多範例見 `examples/demo_plots.py`。

## 後處理/視覺化 (frame2d.plotting)

```python
from frame2d.plotting import plot_structure, plot_loads, plot_diagram, plot_deformed, plot_all
```

- `plot_structure(frame)` — ① 結構尺寸圖 (節點/桿件編號、支承符號)
- `plot_loads(frame)` — ② 受力圖 (點載重箭頭+分佈載重箭頭)
- `plot_diagram(frame, result, kind='N'|'V'|'M')` — ③④⑤ 軸力/剪力/彎矩圖
- `plot_deformed(frame, result)` — ⑥ 變形圖 (自動抓合理放大倍率)
- `plot_all(frame, result)` — 六張一次畫在 2x3 網格

繪圖只是把 `SolveResult` 的數字換一種視角呈現,不會重新計算任何力學——
所有圖表的數值來源都是 `solve()` 算出來的同一份結果。

**注意事項:**
- 標題用英文,避免不同環境(手機/Termux/伺服器)字型缺 CJK 字形時變成方框
- 變形圖用 Hermite cubic 內插,節點值精確;若某根桿件跨中有分佈載重
  且只用單一元素代表整根桿件,內插出來的跨中撓度會略微低估真實下垂量
  (真實解在均佈載重下是四次多項式),想要更精確可以把該桿件切成多個元素

## DOF 設計note

`Frame2D.dofs_of(node_id)` 目前用 `3*node_id` 直接算,假設節點id從0連續編號。
這是刻意的MVP簡化——真正的 DOFManager(支援不連續id、truss少一個自由度、
release/hinge)留到「第二階段: 加入 Truss member」時才需要,現在硬做只是
「為了還沒出現的需求先付架構成本」。呼叫端一律透過 `dofs_of()`,不要自己算
`3*node_id`,將來換底層實作時呼叫端不用改。

## 支承 — fixed/pin/roller/沉陷 統一API

`Support(node, ux=None/0.0/數值, uy=..., rot=...)`——`None`=該方向自由,
`0.0`=固定在原位,非零數值=強制位移(沉陷/施工誤差分析)。`fix()`/`pin()`/
`roller_y()`是這個的簡寫(內部都是0.0), 通用的`f.support(node, ux=, uy=,
rot=)`可以直接設定沉陷量。核心公式從`K_ff u_f = F_f`(劃掉拘束自由度)
改成`K_ff u_f = F_f - K_fc u_c`(u_c放已知的指定位移值), fix/pin/roller
因為u_c全部是0.0, 結果完全不變(已跑過全部既有測試確認零回歸)。
驗證: 靜定結構沉陷應該零內力(結構學基本性質), 一次靜不定梁沉陷對照
傾角變位法經典公式 M=-3EIΔ/L², 見`tests/test_support_displacement.py`。

**注意**: 這次改動過程中抓到一個真的bug——`plotting.py`畫支承符號時原本用
`if support.rot and support.ux and support.uy`判斷"是否為固定端", 但
`Support`從布林值改成數值後, `fix()`設的是`0.0`, 在Python布林判斷裡`0.0`
是假值, 導致固定端符號會被誤判成別種支承圖示。已修正成`is not None`判斷,
並確認三種支承符號(fixed/pin/roller)視覺上都正確。

## 桿件內部集中力/力矩 (不在節點上)

`f.member_point_load(member, a, fx=0.0, fy=0.0, m=0.0)`——`a`是距離
`node_i`沿桿軸的局部座標(0<=a<=L),`fx`(局部軸向)、`fy`(局部橫向,跟
`distributed_load`的w同一套正負號慣例)、`m`(逆時針為正)三個分量都是選填。
**不會新增節點**,用等效節點載重(work-equivalent, 跟`distributed_load`同一套
機制)處理,固定端反力公式用sympy從梁的微分方程重新推導(不是憑記憶抄書),
並用簡支梁的決定性反力(純靜力學, 不依賴公式本身對不對)交叉驗證過,見
`tests/test_member_point_load.py`。

## 局部段均佈載重 (不用整根桿件都有)

`distributed_load(member, w, w_end=None, x_start=None, x_end=None)`——不給
`x_start`/`x_end`時維持原本「整根桿件都有」的行為(向下相容); 給定的話
只在桿件局部座標`[x_start, x_end]`範圍內有載重。公式推導方式: 對已經驗證過的
`fixed_end_forces_point_load()`在`[x_start,x_end]`區間做高斯積分(6點高斯,
對這種低次多項式被積函數是機器精度的數值精確解, 不是近似), 不手動謄寫
龐大的封閉式展開式(降低抄寫出錯風險)。退化情況(x_start=0, x_end=L)已驗證
跟既有全長UDL公式精確一致, 局部段案例用簡支梁決定性反力交叉驗證, 見
`tests/test_partial_udl.py`。

## 分佈載重方向慣例 — 重要note

`distributed_load(member, w)` 的 w 是「沿桿件局部 +y 方向」為正,局部座標系
由 `elements.member_geometry()` 用 `atan2(dy,dx)` 算出(標準右手系,angle=0時
局部y=全域y)。w<0 就是物理上的「向下」,直接照直覺用即可。

這個慣例已經三方驗證過,而且**全部直接吻合、不需要任何翻轉或特例**:
- sd_framework/anastruct(11個案例、46個數值,見下方驗證狀態表)
- SW FEA(Android app)的 Case-08 報告(純點載重版 + 含均佈載重版,共18個反力
  分量,見`test_case08_vs_swfea.py`)——P1/P2 水平點載重沿全域 +x(向右),
  均佈載重沿全域 -y(向下),反力 Rx/Ry/M 三個分量、節點位移 dX/dY 全部直接對上。

開發過程中一度誤判過 SW FEA 的方向慣例(見git歷史), 原因是: (1) 一開始用
app 對話框裡轉盤圖示的視覺猜測角度方向, 沒有實際根據; (2) PDF反力數字抄錄時
Ry/M 的正負號打反。後來直接照 app 畫面上實際畫出來的箭頭方向重建、並且用
「節點位移方向」(dX/dY, 比反力更不容易受慣例混淆的物理量)交叉確認後,才發現
三方本來就是同一套標準慣例, 之前的"SW FEA內部不一致"是誤判, 已更正。

## 驗證狀態

| 案例 | 方法 | 結果 |
|---|---|---|
| 懸臂梁點載重 | 解析解 (PL³/3EI) | rel_err ~1e-16 |
| 簡支梁均佈載重 | 解析解 (5wL⁴/384EI, wL²/8) | rel_err ~1e-16 |
| Case-01~08 (共11案例) | 使用者自己的 sd_framework.py 本體直接執行 | rel_err < 1e-3 (46個獨立數值) |
| Case-08 兩層兩跨鋼架 (純點載重+含UDL) | SW FEA 第三方工具 | 直接吻合, 零翻轉 (18個反力分量) |
| 桁架(truss)單桿軸力 | 解析解 (PL/EA) | rel_err ~1e-10 |
| 桁架(truss)對稱雙桿撐架 | 節點法(method of joints)獨立手算 + SW FEA第三方app | rel_err ~1e-6, 反力/軸力皆吻合(含受壓大小方向) |
| 纜線(cable)受拉/受壓自動鬆弛 | 跟truss結果交叉比對(受拉時) + 整體力平衡檢查 | 完全一致, 零殘差 |
| 懸臂梁斜張橋(frame+truss/cable混合) | 使用者自己的力量法筆記本(含塔柱彎曲修正版) | rel_err < 0.01 (5條纜索張力, 含受壓) |
| 桿件中間集中力(簡支梁反力) | 純靜力學決定性反力 | rel_err ~1e-16 |
| 桿件中間集中力(懸臂梁撓度) | 解析解 (Pa²(3L-a)/6EI) | rel_err ~1e-16 |
| 桿件中間集中力矩(簡支梁反力) | 純靜力學決定性反力 | rel_err ~1e-16 |
| 桿件中間軸向點載重 | 彈簧串聯解析解 | rel_err ~1e-16 |
| 局部段均佈/梯形載重(簡支梁反力) | 純靜力學決定性反力 | rel_err ~1e-6 |
| 局部段公式退化(c=0,d=L) | 既有全長UDL公式 | rel_err ~1e-14 |
| Phase2局部段均佈載重(11個位置逐點) | SW FEA 第三方工具(SF/BM/dY) | 最大誤差 SF~0, BM~0, dY~0.004mm |
| 靜定梁支承沉陷 | 結構學基本性質(靜定結構沉陷零內力) | abs_err ~1e-15 |
| 一次靜不定梁支承沉陷 | 傾角變位法經典公式 (M=-3EIΔ/L²) | rel_err ~1e-16 |

## 桁架(truss)元素

`f.add_truss(id, node_i, node_j, section)` 建立兩端鉸接、只傳軸力的桁架元素
(跟預設的 `add_member` 樑柱元素共用同一個 `Frame2D`/`Section` API, 差別只在
`member_type='truss'`, 局部勁度矩陣的彎曲block是零矩陣)。純桁架節點(只連接
truss桿件、沒有frame桿件)的轉角自由度沒有任何勁度貢獻, solver會自動偵測並
跳過(不會變成奇異矩陣), 但如果對這種節點外加彎矩會直接報錯(沒有任何桿件
能抵抗該彎矩)。

**分佈載重(distributed_load)不能加在truss/cable桿件上**——兩者都沒有彎曲勁度,
`fixed_end_forces_udl` 假設的是frame元素的彎曲能力, 硬加會直接報錯。要模擬
自重, 改成在兩端節點各加一半重量的 `point_load`。

繪圖時 `plot_structure` 會自動把truss/cable桿件畫成虛線(frame桿件是實線),
支承符號請用 `pin()`(鉸接三角形), 不要用 `fix()`(固定端牆面排線)——桁架/
纜線節點本來就沒有轉角勁度, 用pin在視覺上才符合物理意涵, 也才會跟SW FEA
這類工具畫出來的鉸接符號一致。`plot_deformed` 會在每個有位移的節點旁標出
實際位移量(Δ=...), 方便快速比對材料/斷面設定是否正確, 不用只看變形形狀。

## 纜線(cable)元素 — 只能受拉, 自動處理鬆弛

`f.add_cable(id, node_i, node_j, section)` 建立纜線元素——跟truss一樣兩端
鉸接、只傳軸力, 但**只能受拉**。如果某條cable在某個載重組合下該受壓
(物理上代表它會鬆弛退出作用), `solve()` 會自動偵測、移除該桿件的勁度貢獻、
重新求解, 反覆直到沒有cable受壓為止(經典的tension-only member迭代解法)。
`SolveResult.slack_cables` 是最終判定鬆弛的cable id集合, 每個
`MemberResult.slack` 也會標記該桿件這次是否被判定鬆弛(鬆弛時
`end_forces_local`全部是0)。

如果移除所有鬆弛cable之後結構變成機構(無法承受載重), 或超過預設20次迭代
仍未收斂, 會拋出`RuntimeError`並說明原因——這通常代表模型設計本身有問題
(例如某個節點的所有cable在這個載重下都會鬆弛), 不是求解器的bug。

truss跟cable的差別: **truss可以同時承受拉力跟壓力**(適合撐架、桁架的斜撐桿),
**cable只能受拉**(適合真正的纜線、吊索——這是斜張橋的纜線該用的元素類型,
不能直接用truss, 否則載重方向不對時會算出不合理的壓力)。

`test_all_slope_deflection_cases.py` 把 slope_deflection_framework 的 Case-01~08(含 04.5/06.5/07.5,共11個案例)全部轉成 Frame2D 模型,幾何/支承/載重都是照 `samples/model_*.py` 的 `draw_geometry()` 讀出來的真實定義,不是憑印象猜的;正確答案是直接執行使用者的 `sd_framework.py` 本體(`SlopeDeflectionSolver._solve_core()`)算出來的,不是讀 notebook/PDF 截圖轉錄。

下一步可以把 slope_deflection_framework 的 Case-01~08 全部轉成 Frame2D 模型,
系統化地跑一輪回歸測試——這個框架已經證明了同時處理點載重+分佈載重+
多節點多桿件的正確性,可以正式當作那個對照工作的求解器核心。
