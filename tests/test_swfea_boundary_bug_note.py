"""
驗證案例20: SW FEA「鉸接位置恰好=0」的邊界bug, 以及frame2d為什麼不受影響

背景: 使用者發現SW FEA的pin2_1案例(F4梁在distance=0處設鉸接, 剛好在
節點5)算出的結果內部矛盾(F4沒有任何桿件內部載重, 但軸力/剪力沿全長
應該是常數卻不是, 詳見test_two_story_release_vs_swfea.py的原始發現)。
使用者測試把鉸接位置從distance=0改成distance=0.0001(偏移0.1mm), 這個
不連續消失了, 而且新的反力/逐點資料精確吻合frame2d原本(在distance=0
沒有offset)算出的答案。

**這證實了兩件事**:
1. SW FEA在「鉸接distance恰好=0(剛好在節點上)」這個邊界情況有實作bug,
   偏移一點點(0.0001m)就正常了。這不是我們冤枉SW FEA——F4軸力/剪力
   不連續是物理上不可能發生的事, 不管哪一種鉸接語意解讀都說不通,
   而偏移後的結果精確吻合frame2d, 雙重確認。
2. frame2d的release_i/release_j機制不受這類bug影響, 而且架構上就不會有
   ——因為release_i/release_j是直接指定"哪個端點"的布林旗標(離散的
   member_id + i端/j端選擇), 根本沒有"距離"這種連續參數, 不存在
   "剛好等於0"跟"差一點點"的模糊地帶可以出問題。

這對後續驗證的意義: SW FEA仍然是有用的外部參照工具, 但如果要設定鉸接
剛好在節點上, 建議之後都用0.0001這種小offset, 避免踩到這個邊界bug。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}


def build():
    f = Frame2D()
    for nid, (x, y) in NODES.items():
        f.add_node(nid, x, y)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 3, 's')
    f.add_member(2, 3, 1, 's')
    f.add_member(3, 4, 2, 's')
    f.add_member(4, 5, 4, 's', release_i=True)   # 梁在節點5端釋放(=distance 0, 不用offset)
    f.add_member(5, 5, 3, 's')
    f.fix(0)
    f.fix(1)
    f.point_load(2, fx=10.0)
    f.point_load(4, fx=15.0)
    return f


print("=== frame2d(release_i, 剛好在節點5端) vs SW FEA pin2_1_1(offset 0.0001m) ===")
f = build()
r = solve(f)
sw_reactions = {0: (-11.589, -15.354, 32.752), 1: (-13.411, 15.354, 35.125)}
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    Rx, Ry, M = r.reactions[ux], r.reactions[uy], r.reactions[rot]
    swv = sw_reactions[n]
    print(f"  node{n}: f2d=({Rx:.3f},{Ry:.3f},{M:.3f})  sw(offset)={swv}")
    assert abs(Rx - swv[0]) < 0.01 and abs(Ry - swv[1]) < 0.01 and abs(M - swv[2]) < 0.01
print("  PASS: frame2d(不用offset)精確吻合SW FEA(offset 0.0001m後的結果)\n")

# 逐點BM比對(F0,F2,F3,F4,F5,F6六根桿件x11點, 用內插對齊)
sw_BM = {
    0: [-32.752, -28.116, -23.481, -18.845, -14.210, -9.574, -4.939, -0.303, 4.333, 8.968, 13.604],
    2: [-18.519, -13.155, -7.790, -2.426, 2.939, 8.303, 13.668, 19.032, 24.396, 29.761, 35.125],
    3: [-20.941, -16.564, -12.186, -7.809, -3.431, 0.946, 5.324, 9.702, 14.079, 18.457, 22.834],
    4: [0.000, -2.094, -4.188, -6.282, -8.376, -10.470, -12.565, -14.659, -16.753, -18.847, -20.941],
    5: [-0.000, 1.622, 3.245, 4.867, 6.489, 8.112, 9.734, 11.357, 12.979, 14.602, 16.224],
    6: [36.438, 29.320, 22.202, 15.084, 7.965, 0.847, -6.271, -13.389, -20.507, -27.625, -34.743],
}
member_L = {0: 4.0, 2: 4.0, 3: 4.0, 4: 6.0, 5: 4.0, 6: 6.0}
# frame2d的member id跟SW FEA的F編號對照: F0->0, F2->2, F3->3, F4->4(=frame2d member4), F5->5, F6->1
id_map = {0: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 1}
from frame2d.postprocess import member_internal_forces
max_err = 0
for sw_id, bm_list in sw_BM.items():
    mid = id_map[sw_id]
    L = member_L[sw_id]
    x_full, N, V, M = member_internal_forces(f, r, mid, n=201)
    x_query = np.array([frac * L for frac in np.linspace(0, 1, 11)])
    M_interp = np.interp(x_query, x_full, M)
    err = np.max(np.abs(M_interp - np.array(bm_list)))
    max_err = max(max_err, err)
    print(f"  SW F{sw_id}(frame2d member{mid}): 11點BM最大誤差={err:.4f}")
assert max_err < 0.02, f"BM誤差過大: {max_err}"
print(f"  PASS: 六桿件x11點BM全數吻合(整體最大誤差{max_err:.4f})\n")


# ---- frame2d架構上不受"距離恰好=0"這類邊界bug影響的說明性驗證 ----
print("=== frame2d的release機制架構上不存在這類邊界bug ===")
print("  release_i/release_j是直接指定端點的布林旗標, 沒有連續的距離參數,")
print("  不存在'剛好等於0'跟'差一點點'的模糊地帶。")
print("  (對照: 桿件內部力矩的a位置雖然是連續參數, 但a=0跟a=0.0001的答案")
print("  平滑連續, 沒有離散跳躍, 詳見獨立驗算)")

print("\nPASS: SW FEA邊界bug記錄+frame2d驗證完成")
