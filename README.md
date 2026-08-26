# frame2d — 2D 矩陣位移法(直接勁度法)通用框架分析

從 [slope_deflection_framework](https://github.com/zhixiu0223/slope_deflection_framework)
分出的新專案:傾角變位法 repo 保留「每一步推導透明」的教學價值,這裡專注「任意拓樸都能解」。

## 目前範疇 (刻意設限,避免無限擴張)

- 2D Euler-Bernoulli 樑柱元素 (frame member),每節點3自由度 (ux, uy, rot)
- **尚未支援**: truss member(純軸力元素)、內部鉸接(internal hinge)、
  幾何非線性(P-Delta)、材料非線性——這些等基本框架穩定後再視需求加入,
  不要為了還沒出現的需求先付架構成本

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
```

## DOF 設計note

`Frame2D.dofs_of(node_id)` 目前用 `3*node_id` 直接算,假設節點id從0連續編號。
這是刻意的MVP簡化——真正的 DOFManager(支援不連續id、truss少一個自由度、
release/hinge)留到「第二階段: 加入 Truss member」時才需要,現在硬做只是
「為了還沒出現的需求先付架構成本」。呼叫端一律透過 `dofs_of()`,不要自己算
`3*node_id`,將來換底層實作時呼叫端不用改。

## 分佈載重方向慣例 — 重要note

`distributed_load(member, w)` 的 w 是「沿桿件局部 +y 方向」為正,局部座標系
由 `elements.member_geometry()` 用 `atan2(dy,dx)` 算出(標準右手系,angle=0時
局部y=全域y)。**這跟 SW FEA (Android app) 的方向慣例相反**——比對 Case-08 時
發現要用正值(w=+12000)才能跟SW FEA報告吻合,推測是SW FEA內部用螢幕座標
(Y向下為正)。這不是bug,只是不同工具的輸入慣例不同,使用本套件時請自行
依「local +y 實際指向哪個全域方向」決定w的正負號(可以用
`member_geometry(node_i,node_j)` 印出 angle 來確認)。

同樣地,`SolveResult.reactions` 的整體正負號慣例(定義為 `K@u - F`,即「結構
受到的外力」)也跟SW FEA報告的Rx/Ry/M相反號——比對時記得統一。

## 驗證狀態

| 案例 | 方法 | 結果 |
|---|---|---|
| 懸臂梁點載重 | 解析解 (PL³/3EI) | rel_err ~1e-16 |
| 簡支梁均佈載重 | 解析解 (5wL⁴/384EI, wL²/8) | rel_err ~1e-16 |
| Case-08 兩層兩跨鋼架 | SW FEA 第三方工具 | rel_err < 0.2% (9個反力分量) |

下一步可以把 slope_deflection_framework 的 Case-01~08 全部轉成 Frame2D 模型,
系統化地跑一輪回歸測試——這個框架已經證明了同時處理點載重+分佈載重+
多節點多桿件的正確性,可以正式當作那個對照工作的求解器核心。
