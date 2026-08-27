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

### Phase 1: 桿件中間的集中力 + 集中力矩
最常用、也是驗證「element load → equivalent nodal load」這個矩陣位移法
核心模式的最佳題目。做法照現有 `distributed_load` 的模式:算出固定端反力
→ 等效節點載重疊加進F → 回代時扣回來。**不新增節點**(不切割member,
不污染拓樸),用等效節點載重表示,這樣M(x)/V(x)公式也可以直接沿用現有的
`member_internal_forces`架構(局部段落分開處理)。
驗證: 懸臂梁受任意位置集中力/力矩, 跟解析解(標準懸臂梁公式)比對, 機器精度。

### Phase 2: 局部段均佈載重
`distributed_load` 的 `DistributedLoad` dataclass 這次直接加上 `x_start`,
`x_end` 兩個欄位(預設None=整根桿件, 向下相容現有呼叫方式), 是Phase 1
點載重公式沿桿長積分的推廣, 做完Phase 1後這項不難。
驗證: 簡支梁局部段均佈載重, 解析解比對。

### Phase 3: Support改成「指定值」而非布林值
`Support(node, ux=True/False)` → `Support(node, ux=None/0.0/指定值)`,
`None`=自由, `0.0`=固定在原位, 非零值=強制位移(沉陷/施工誤差分析)。
這個改動統一了fixed/pin/roller/settlement, 不用另外做一個
`SupportDisplacement` API。核心公式: 邊界條件從「劃掉」變成
`K_ff u_f = F_f - K_fc u_c`(u_c為已知的支承位移向量), 改動集中在
`solve.py`, 不影響其他模組。
驗證: 簡支梁其中一端沉陷已知量, 跟解析解比對。

### Phase 4 (較晚, 需要先做DOFManager才能開始): Internal hinge / element release
桿件端點內鉸(例如 A──────○B, B端不傳彎矩)。這個會需要先把
`Frame2D.dofs_of()`的DOF系統從現在的MVP版(`3*node_id`固定)升級成真正
的DOFManager, 因為release會讓同一個節點在不同桿件眼中看到不同的自由度
耦合關係。**先不要在還沒做release前提早做DOFManager升級**——沒有具體
需求驅動的架構升級容易做錯方向。

## 暫緩/刻意不做的項目

- 均佈載重的任意角度/座標系選項: 沒有具體驗證場景前不加欄位
- 斜張橋以外的其他特殊橋型/結構型式: frame+truss+cable已經證明夠通用,
  不需要再做專屬功能證明
- 非線性(大變形/材料非線性): 目前定位是線性靜力分析核心, 不在範疇內
