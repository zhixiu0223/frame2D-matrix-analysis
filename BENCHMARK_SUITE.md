# Benchmark Suite 索引

frame2d 的正確性不是靠單一權威來源背書,是靠多種互相獨立的驗證方式
交叉確認。這份文件把 `tests/` 底下 30 個驗證案例依「驗證對象」分類整理,
方便快速找到「某個功能是用什麼方法驗證過的」,不用逐一打開每個檔案的
docstring。

跑法: `PYTHONPATH=. pytest tests/ -q`(全部 40 個測試, 含
`test_zz_all_script_style_tests.py` 用 subprocess 執行的 script 風格
測試), 或直接 `PYTHONPATH=. python tests/test_xxx.py` 單獨執行某一個
(會印出詳細比對數值方便人眼檢查)。

## 驗證方法分類

frame2d 使用四種互相獨立的驗證方式, 沒有任何一個功能只靠單一方式背書:

| 方法 | 特性 | 用在哪裡 |
|---|---|---|
| **解析解** | 純數學推導, 不依賴任何程式碼 | 懸臂梁、簡支梁等靜定/簡單靜不定結構 |
| **純靜力學(力平衡)** | 不依賴 fixed_end_forces 公式本身對不對 | 桿件內部載重、局部段均佈載重、支承沉陷 |
| **第三方獨立工具**(slope_deflection_framework、SW FEA) | 完全獨立的另一套實作 | 側移剛架、內部鉸接、支承沉陷等多數案例 |
| **兩套獨立內部實作互相驗證**(靜力凝縮 vs DOFManager) | 同一個團隊寫的兩套完全不共用組裝邏輯的程式碼 | 內部鉸接系列 |

## 案例索引(按主題分組)

### 基礎元素(案例 1-2, 6-8)
| 案例 | 檔案 | 驗證對象 | 方法 |
|---|---|---|---|
| 1 | `test_cantilever.py` | 懸臂梁自由端點載重 | 解析解 |
| 2 | `test_simply_supported_udl.py` | 簡支梁均佈載重 | 解析解 |
| 6 | `test_truss.py` | 桁架(truss)元素 | 解析解 + 節點法手算 |
| 7 | `test_cable.py` | 纜線(cable)元素, 鬆弛自動判定 | 純力平衡(含機構偵測) |
| 8 | `test_cable_stayed_bridge.py` | 懸臂梁斜張橋(5條纜索) | 解析解(位移法 vs 力法交叉驗證) |

### 側移剛架(案例 3-5)
| 案例 | 檔案 | 驗證對象 | 方法 |
|---|---|---|---|
| 3 | `test_case08_vs_swfea.py` | 兩層兩跨側移鋼架(九節點/十桿件) | SW FEA(純點載重版+含UDL版, 18個反力分量) |
| 4 | `test_case04_vs_sdframework.py` | 側移單跨剛架 Case-04/04.5 | slope_deflection_framework(notebook結果) |
| 5 | `test_all_slope_deflection_cases.py` | Case-01~08 全部11個案例 | slope_deflection_framework(直接執行對方原始碼, 不是讀PDF) |

### Load System v2: 桿件內部載重(案例 9-14)
| 案例 | 檔案 | 驗證對象 | 方法 |
|---|---|---|---|
| 9 | `test_member_point_load.py` | 桿件內部集中力/力矩 | 簡支梁決定性反力(純靜力學) |
| 10 | `test_partial_udl.py` | 局部段均佈載重 | 簡支梁決定性反力(純靜力學) |
| 11 | `test_support_displacement.py` | 支承強制位移(沉陷) | 靜定結構性質(沉陷不產生內力) |
| 12 | `test_deformed_shape_accuracy.py` | 變形曲線精確度(M(x)雙重積分法) | 節點分割法 + 200段密網格 |
| 13 | `test_partial_udl_vs_swfea.py` | 局部段均佈載重 | SW FEA(phase-02案例) |
| 14 | `test_support_displacement_vs_swfea.py` | 支承強制位移 | SW FEA(phase-03-rot/sinking案例) |

### 內部鉸接(Internal Hinge / Release)系列(案例 15-26)
這是整個 benchmark suite 裡份量最重的一組, 因為過程中不斷發現 SW FEA
在這個功能上有數值邊界問題, 逼出了一整套「用短桿段模擬鉸接的數值安全
邊界」分析。

| 案例 | 檔案 | 驗證對象 | 方法 |
|---|---|---|---|
| 15 | `test_element_release.py` | 內部鉸接基礎功能(靜力凝縮) | Gerber梁純靜力學 |
| 16 | `test_element_release_vs_swfea.py` | 內部鉸接 | SW FEA(phase-04案例) |
| 17 | `test_dofmanager_vs_condensation.py` | 靜力凝縮 vs DOFManager兩套獨立實作 | 互相交叉驗證(誤差~1e-13浮點精度) |
| 18 | `test_two_story_release_vs_swfea.py` | 兩層樓框架, 三根桿件交界的4種釋放組合 | SW FEA |
| 19 | `test_combined_release_loads_vs_swfea.py` | 內部鉸接+桿件內部載重+局部UDL+桿件力矩全部疊加 | SW FEA(最複雜的單一組合案例) |
| 20 | `test_swfea_boundary_bug_note.py` | **SW FEA「鉸接距離恰好=0」邊界bug** 的發現與記錄 | SW FEA自我內部矛盾(F4軸力/剪力不連續) |
| 21 | `test_triple_release_swfea_bug.py` | **SW FEA「三桿件同時全釋放」bug**: 用pin2_4≡pin2_5的數學證明抓出 | 純數學證明(必然相等) + SW FEA |
| 22 | `test_offset_sensitivity_analysis.py` | 短桿段offset敏感度分析(第一版, 過程中抓到自己的模型bug) | 自我比對(節點分割 vs release_i/j) |
| 23 | `test_pin2_1_to_5_clean_comparison.py` | 乾淨版本(不做節點分割)對SW FEA的殘差規律 | SW FEA |
| 24 | `test_double_member_vs_clean_release.py` | 雙桿件模擬(精確offset) vs 乾淨版本, 哪個更準 | SW FEA(結論: 乾淨版本更準) |
| 25 | `test_meaningful_length_split_vs_swfea.py` | 雙桿件模擬, 用有意義長度(0.8~1.2m)驗證模型設定無誤 | SW FEA |
| 26 | `test_tiny_offset_reproduction.py` | 用已驗證雙桿件方法重新套用SW FEA原始極小offset | SW FEA + 自我比對 |
| 27 | `test_a_over_L_sweep_boundary.py` | **精確定位a/L數值安全邊界**(懸崖在a/L≈4~5e-5) | 自我比對(雙桿件 vs 乾淨release, 20+點掃描) |
| 28 | `test_a_over_L_1e1_to_1e4_vs_swfea.py` | a/L=1e-1~1e-4(安全區內)逐點對照 | SW FEA |

### 其他功能(案例 29-30)
| 案例 | 檔案 | 驗證對象 | 方法 |
|---|---|---|---|
| 29 | `test_result_api.py` | Result API便利查詢介面(純介面整理) | 對照底層函式手動查詢+解析解 |
| 30 | `test_sloped_roof_global_udl.py` | `distributed_load(direction='global_y')`(斜屋頂重力/雪載重) | SW FEA(slop-roof案例, 4桿件x11點N/V/M) |

## 這套 benchmark suite 意外抓到的 SW FEA 問題

開發過程中, 有兩次獨立發現 SW FEA(第三方參照工具)本身在特定邊界情況
下有問題, 不是 frame2d 的模型設錯:

1. **鉸接距離恰好等於 0**(案例20): F4桿件在distance=0處設鉸接, SW FEA
   算出的軸力/剪力沿全長不連續(物理上不可能, 該桿件沒有任何內部載重)。
   偏移0.0001m後恢復正常, 而且精確吻合frame2d(不用offset)的答案。

2. **同一節點三根桿件同時全釋放**(案例21): 用數學證明"F2+F6+F5三個都
   釋放"必然等於"只釋放F2+F5(梁自然變成唯一DOF使用者)", frame2d驗證
   誤差~1e-12(浮點精度), 但SW FEA報告的兩個案例反力差異高達~20, 而且
   節點轉角出現-101 rad的荒謬數字。後續(案例26)用同款offset手法在
   frame2d自己身上也重現了同量級的"應該相等卻不相等"現象, 確認這是
   "短桿段模擬鉸接"這個建模手法本身的數值脆弱性, 不是SW FEA獨有的
   神秘bug。

這兩個發現最終收斂成案例27(a/L敏感度掃描), 精確定位出數值安全邊界在
a/L≈4~5×10⁻⁵, 並寫入README.md正式的使用準則: 鉸接在桿件端點用
release_i/release_j(不受這個問題影響); 鉸接在桿件中段需要節點分割時,
offset要用有意義的長度(建議>=桿長的1%), 不要用極小值逼近端點模擬。
