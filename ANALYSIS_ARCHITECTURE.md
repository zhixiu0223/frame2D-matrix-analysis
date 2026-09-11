# frame2D-matrix-analysis 架構現況

這份文件只記錄**現況**(哪個檔案的哪個函式實際在做什麼),不是理想架構
的規劃書——理想架構的討論放在文末「規劃中」那一節,並且明確標示哪些
還沒做。每一個「已驗證」的陳述, 都對應到 `tests/` 底下實際存在、會被
CI 跑到的測試, 不是憑印象寫的。

## 一句話總結現況

現在**有**一個統一的 `analyze_pushover(geometry=, solver=)` 入口
(`frame2d/analyze.py`),但它是**純dispatch層**,底下仍然是四個各自
獨立驗證過的既有函式(`run_pushover()`、`run_pushover_converged()`、
`run_pushover_newton()`、`run_pushover_corotational_oneshot()`)——
`analyze_pushover()`本身沒有新的物理或求解邏輯,只是把「選幾何處理
方式」跟「選求解策略」這兩個維度,從原本「綁死在挑三個求解器名稱
之一」,變成可以獨立組合的兩個參數。另外還有一個更早、不屬於
pushover 的線性 P-Delta 求解函式(`solve_pdelta()`),目前**還沒**
被納入這個 dispatch 層。

## 四層模型(用來分類現有程式碼,不是新設計)

```
1. MODEL(frame2d/model.py)
   Frame2D / Node / Member / Section / 各種載重(PointLoad,
   DistributedLoad, MemberPointLoad)/ 支承條件

2. PHYSICS / FORMULATION(材料 + 幾何非線性怎麼處理)
   材料: 彈性(預設)、雙折線塑鉸(frame2d/hinge.py)
   幾何: 線性(不管幾何非線性)
         Global P-Delta(elements.local_geometric_stiffness(), 套在
           "絕對座標"v/θ上——frame2d/dofmanager.py, frame2d/pushover.py)
         幾何更新(每步用變形後幾何重算L/角度, 仍是小角度公式——
           frame2d/pushover.py的geometry_update選項)
         Co-rotational(精確拆開剛體轉動+局部變形——
           frame2d/corotational.py)
         Local P-Delta(套在co-rotational的自然座標θ_def上, 跟
           Global P-Delta同一條公式elements.local_geometric_
           stiffness(), 只是套用的座標系不同——frame2d/corotational.py
           的use_pdelta參數)

3. SOLVER / STRATEGY(怎麼把方程式解出來)
   直接單次解(np.linalg.solve一次, 不疊代——frame2d/solve.py,
     event-to-event每一步的作法)
   Picard(不動點)疊代(疊代到幾何/軸力自洽, 不是殘餘力收斂——
     solve_pdelta()、run_pushover_converged())
   Newton-Raphson(疊代到殘餘力f_ext-f_int收斂——
     run_pushover_newton()裡的_newton_iterate())

4. POSTPROCESS(frame2d/postprocess.py, webapi/diagrams.py)
   N/V/M圖、變形圖、反力、塑鉸狀態、容量曲線、收斂旗標
```

## 現有函式對照表(這是本文件最重要的部分, 全部對照過原始碼)

| 函式 | 檔案 | 材料 | 幾何 | 求解策略 | 备注 |
|---|---|---|---|---|---|
| `solve()` | frame2d/solve.py | 彈性 | 線性 | 直接單次解 | 最基本的線性分析 |
| `solve_pdelta()` | frame2d/dofmanager.py | 彈性 | Global P-Delta | Picard疊代 | 不是pushover, 是「分析」選單裡的「P-Delta」線性分析選項 |
| `run_pushover()` | frame2d/pushover.py | 塑鉸 | 線性/Global P-Delta/幾何更新(`use_pdelta`、`geometry_update`兩個**互相獨立**的bool) | 直接單次解(event-to-event, 每個塑鉸狀態凍結區間內只解一次) | 三種pushover求解器裡最快的 |
| `run_pushover_converged()` | frame2d/pushover.py | 塑鉸 | 同上(`use_pdelta`、`geometry_update`一樣互相獨立) | Picard疊代(疊代到幾何/軸力自洽) | **驗證過**: 關掉`use_pdelta`跟`geometry_update`時, 結果跟`run_pushover()`逐位元一致(相對誤差0.0)——因為疊代的對象(幾何/軸力)根本沒有東西在變, 疊代等於沒疊代 |
| `run_pushover_newton()` | frame2d/newton.py | 塑鉸 | Co-rotational(內建, 不能關)+ Local P-Delta(`use_pdelta`, 獨立開關) | Newton-Raphson疊代 | 目前**不支援**release端; 支援重力預載(內部先疊代解一次純重力平衡, 不是天真假設u=0) |
| `run_pushover_corotational_oneshot()` | frame2d/newton.py | 塑鉸 | Co-rotational(內建)+ Local P-Delta(`use_pdelta`) | 直接單次解(event-to-event式) | **驗證性函式**(見「規劃中」第1點), 不修改corotational.py; 重用newton.py既有函式; 已驗證跟run_pushover_newton()同步長收斂到一致結果(誤差隨步長縮小) |

## 已經做到的「分層」證據(不是理想化的說法, 有程式碼可以查)

`frame2d/corotational.py`的import只有`numpy`、`hinge_bending_stiffness`
(材料模型)、`local_geometric_stiffness`(P-Delta公式)——**完全不知道
Newton疊代迴圈的存在**。也就是說, 對co-rotational這條路徑,「物理層」
(corotational.py)跟「求解層」(newton.py的`_newton_iterate()`)已經是
乾淨分開的兩個檔案, `newton.py`只是corotational.py眾多可能的呼叫端之一
——這是2025-09-11對話裡確認過的事實, 不是規劃。

`run_pushover()`跟`run_pushover_converged()`這條舊路徑, 物理(P-Delta/
幾何更新怎麼組裝勁度矩陣)目前**還跟**求解迴圈(event-to-event的降伏
事件精確定位邏輯)寫在同一個檔案、同一組函式裡(`_assemble_stiffness_
with_hinges()`直接被pushover.py的主迴圈呼叫, 不是一個可以被其他求解
策略重用的獨立模組)——這部分**還沒有**做到跟corotational.py同等程度
的分離。

## 現有的9種(理論上)組合, 哪些已驗證存在

| # | 求解函式 | Global/Local P-Delta | 幾何更新/Co-rotational | 疊代到? |
|---|---|---|---|---|
| 1 | run_pushover | 關 | 關 | (不疊代) |
| 2 | run_pushover | 開(global) | 關 | (不疊代) |
| 3 | run_pushover | 開(global) | 開(幾何更新) | (不疊代) |
| 4 | run_pushover_converged | 開(global) | 開(幾何更新) | 幾何/軸力自洽 |
| 5 | run_pushover_converged | 開(global) | 關 | 幾何/軸力自洽 |
| 6 | run_pushover_converged | 關 | 開(幾何更新) | 幾何/軸力自洽(意義有限, 見下) |
| 7 | run_pushover_newton | 關(local) | co-rotational(內建) | 殘餘力收斂 |
| 8 | run_pushover_newton | 開(local) | co-rotational(內建) | 殘餘力收斂 |
| 9 | run_pushover_corotational_oneshot | 關/開(local) | co-rotational(內建) | (不疊代, event-to-event式) |

組合6意義有限: 沒有P-Delta回饋進勁度矩陣, 疊代檢查的只是「幾何有沒有
自洽」, 沒有軸力效應可以疊代, 實務上少見, 但機制上存在、不會crash。

**組合9已完成、已驗證**(見下方「規劃中」第1點的完整記錄)——介於
組合3(幾何更新, 小角度, 不疊代)跟組合7(co-rotational, 疊代)之間,
用co-rotational的精確幾何, 但求解策略換成跟event-to-event同等級的
「一次到位」。這是驗證「物理層/求解層能不能自由組合」最直接的試金石:
不修改`corotational.py`一行程式碼, 只靠重用`newton.py`既有的
`_assemble_global()`等函式就做出這個新求解策略, 證明了這個分離
架構是真的撐得住的, 不是巧合。

## 規劃中(還沒做, 明確標示)

1. **`run_pushover_corotational_oneshot()`**(對應上面第9種組合)——
   **已完成, 已驗證**(2025-09-11): 用跟`run_pushover_newton()`完全同一套
   co-rotational元素公式(`_assemble_global()`/`_member_moments()`等既有
   函式, **沒有修改`corotational.py`一行程式碼**), 求解策略換成
   event-to-event式的"一次到位"。

   驗證過程中抓到一個真實bug(狀態更新順序): 原本先把塑鉸標記成"已
   降伏", 才用(已經變軟的)新勁度回頭去算"這一步累積了多少彎矩",
   把這一步裡其實還是彈性的那一段也錯誤地套用降伏後的軟化勁度去算,
   導致彎矩系統性低估、且隨步數累積——已修正成"先用還沒標記降伏的
   狀態把這一步該有的彎矩算完, 再標記降伏"。

   修正後, 依系統性隔離協議驗證(純彈性→有HingeState但永不降伏→
   跨越降伏, 逐層排除變因): 純彈性案例縮小步長呈現乾淨的一階收斂;
   有HingeState但Mp設超大永不降伏時, 結果跟純彈性相對誤差2.77e-6
   (幾乎完全一致); 真正跨越降伏時, 用**同樣細的步長**分別跑
   `run_pushover_corotational_oneshot()`跟`run_pushover_newton()`
   比較, 步長縮小5倍, 兩者差距從0.44%縮小到0.11%——**沒有殘留的
   系統性偏差**。

   除錯過程中的重要教訓(值得記錄, 避免重蹈覆轍): 一開始誤以為存在
   "11%的殘留系統性偏差"(縮小步長, 兩者差距反而穩定在11%附近不再
   縮小), 後來發現是比較方式本身有問題——拿`run_pushover_
   corotational_oneshot()`(細步長, 已經收斂)去對照`run_pushover_
   newton()`(**粗步長**, 沒有收斂完全)這個當作"基準"的參照值, 而
   Newton本身的降伏事件偵測是在疊代"內部"做(見newton.py docstring
   的已知限制), 精確度同樣隨步長縮小——粗步長的Newton自己都還沒收斂
   到真正解, 拿來當基準自然會有誤導性的"殘留偏差"。教訓: **比較兩個
   求解器的精確度時, 要用同樣細的步長, 不能把粗步長的某一方當成
   "已知正確"的基準**。

2. **`analyze_pushover()`統一入口**——**已完成, 已驗證**(2025-09-11):
   `frame2d/analyze.py`, 純粹的dispatch層, **沒有任何新的物理或求解
   邏輯**, 依`geometry=`(`'linear'`/`'updated'`/`'corotational'`)、
   `solver=`(`'direct'`/`'picard'`/`'newton'`)兩個獨立參數, 呼叫
   `frame2d.pushover`/`frame2d.newton`裡對應的既有函式。**沒有修改、
   也沒有取代**任何一個現有函式——這是刻意的安全網設計: `tests/
   test_analyze_pushover.py`驗證過, 上面表格裡9種組合當中7種可以
   直接對照(組合6意義有限, 沒有寫對應測試; 組合7/8/9額外測過), 每
   一種都跟直接呼叫底層函式的結果逐位元/浮點數誤差為0一致, 加上
   kwargs正確轉傳、不支援的組合(例如`geometry='corotational'`配
   `solver='picard'`)正確拒絕、不是靜默退化。
3. **完全獨立於pushover之外的統一入口**——`solve()`/`solve_pdelta()`
   目前也還沒被納入這個dispatch層, 未來如果要做真正一致的
   `analyze()`, 這兩個也要一起考慮進去。
4. **動力分析(Newmark/HHT等時間積分)**——完全還沒開始, 但如果第1、2
   點做成功, `_newton_iterate()`這種"疊代到殘餘力收斂"的核心邏輯,
   理論上可以被動力分析的每個時間步重用(Newton法不是pushover專屬的,
   只是"用來解非線性平衡方程式的數值方法")——這是為什麼現在不建議把
   `newton.py`寫死成"pushover solver"的原因, 但這一點目前只是設計
   上預留空間, 沒有任何動力分析的程式碼。
