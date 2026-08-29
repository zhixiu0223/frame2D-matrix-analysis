"""
驗證案例19: 單層門型鋼架, 桿件內部點載重+局部段均佈載重+桿件力矩+內部鉸接
全部組合在一起(對照 phase-04-internal-pin_3_1.frame)

這是目前最複雜的單一模型組合測試: F1(梁)同時有桿件內部點載重、局部段
均佈載重、J端內部鉸接、以及I端(a=0, 剛好在節點2位置)的力矩; F2(柱)
另外加了桿件內部力矩。

**單位換算的發現**: .frame檔案裡 point_loads/distributed_loads/moments/
internal_hinges 這幾張表的 start_dist/end_dist 欄位都是"桿件長度的比例"
(0~1之間), 不是絕對公尺數。PDF報告裡均佈載重跟力矩兩張表有正確換算成
公尺顯示, 但**點載重那張表疑似沒有乘上桿長, 直接把比例值當公尺數顯示**
(顯示0.333, 但真實物理位置是0.3333*6=2.0m)。

**重要更正(2026-08-30)**: 這個檔案原本記錄「桿件內部力矩, PDF/UI標示
「逆時針」對應到的方向, 跟我們認知的CCW正相反」——**這個結論是錯的,已經
撤回**。真正的原因是 frame2d 自己的 `fixed_end_forces_point_moment()`
公式當時有正負號bug(把point_load推導時"sympy原始結果要整體取負號"這個
規律誤套用到力矩案例, 但兩者的sympy推導約定不是同一回事)。這個bug讓
frame2d內部兩種表示同一個力矩的方式(member_point_load的m參數 vs 真正
的節點point_load的m參數)彼此不一致, 用「懸臂梁自由端加外力矩」(固定端
反力矩應為-M0, 梁彎矩應為常數+M0, 這是不依賴這條公式本身的物理事實)
跟「節點分割法」(在力矩位置插入真實節點, 改用已驗證過的point_load的m)
互相佐證抓到並修正, 見 tests/test_member_point_load.py 案例C的更正記錄。
修正後, SW FEA UI標示的「Counter-Clockwise」直接對應m=+10(不用反號),
跟我們的CCW正慣例完全一致, 之前的"UI標示相反"是誤判, 已更正。

**兩套求解器的差異(這是這次題目的重點, 使用者原本問"凝縮模型和
DOFmanager都可以嗎?")**: F1桿件同時有release_j跟(桿件內部點載重+局部段
均佈載重+桿件力矩), 這個組合solve_condensation()(靜力凝縮, 參考實作)
會直接報錯(這是已知、有文件記錄的限制), 只有solve()(=solve_dofmanager(),
主要求解器)可以直接處理。
"""
import numpy as np
from frame2d import Frame2D, solve, solve_condensation
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01


def build():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')                     # F0: 左柱
    f.add_member(1, 2, 3, 's', release_j=True)     # F1: 梁, J端(node3)鉸接
    f.add_member(2, 3, 1, 's')                     # F2: 右柱
    f.fix(0)
    f.fix(1)
    f.member_point_load(1, a=2.0, fy=-10.0)                  # F1: a=2.0m, 10kN向下
    f.distributed_load(1, w=-10.0, x_start=2.0, x_end=4.0)   # F1: [2,4]這段, 10kN/m向下
    f.member_point_load(2, a=2.0, m=10.0)                    # F2: a=2.0m, 10kNm逆時針
    f.member_point_load(1, a=0.0, m=10.0)                    # F1: a=0(=node2位置), 10kNm逆時針
    return f


print("=== 主要求解器 solve() (DOFManager) ===")
f = build()
r = solve(f)
sw_reactions = {0: (3.765, 20.434, -2.457), 1: (-3.765, 9.566, 5.060)}
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    Rx, Ry, M = r.reactions[ux], r.reactions[uy], r.reactions[rot]
    swv = sw_reactions[n]
    print(f"  node{n}: f2d=({Rx:.4f},{Ry:.4f},{M:.4f})  sw={swv}")
    assert abs(Rx - swv[0]) < 0.01 and abs(Ry - swv[1]) < 0.01 and abs(M - swv[2]) < 0.01
print("  PASS: 反力精確吻合SW FEA\n")

# 逐點BM比對(3桿件x11點, 用內插對齊, 因為桿件內部有點載重/局部段均佈載重
# 時member_internal_forces()會自動插入額外取樣點, 不能直接按陣列位置比對)
sw_BM = {
    0: [2.457, 0.951, -0.555, -2.061, -3.567, -5.073, -6.579, -8.085, -9.591, -11.097, -12.603],
    1: [-22.603, -10.343, 1.918, 14.178, 21.638, 23.698, 22.159, 17.219, 11.479, 5.740, 0.000],
    2: [0.000, 1.506, 3.012, 4.518, 6.024, -2.470, -0.964, 0.542, 2.048, 3.554, 5.060],
}
member_L = {0: 4.0, 1: 6.0, 2: 4.0}
max_err = 0
for mid, bm_list in sw_BM.items():
    L = member_L[mid]
    x_full, N, V, M = member_internal_forces(f, r, mid, n=201)
    x_query = np.array([frac * L for frac in np.linspace(0, 1, 11)])
    M_interp = np.interp(x_query, x_full, M)
    err = np.max(np.abs(M_interp - np.array(bm_list)))
    max_err = max(max_err, err)
    print(f"  member{mid}: 11點BM最大誤差={err:.4f}")
assert max_err < 0.01, f"BM誤差過大: {max_err}"
print(f"  PASS: 三桿件x11點BM全數吻合(整體最大誤差{max_err:.4f})\n")

m_hinge = r.member_results[1].end_forces_local[5]
assert abs(m_hinge) < 1e-9
print(f"  鉸接端(F1的J端)彎矩: {m_hinge:.2e} (精確為0)\n")


# ---- 額外驗證: 桿件a=0的力矩 跟 真正節點力矩, 現在應該完全等價 ----
print("=== 額外驗證: 桿件a=0力矩 vs 真正節點力矩, 應該完全一致 ===")
fA = Frame2D()
fA.add_node(0, 0, 0); fA.add_node(1, 6, 0); fA.add_node(2, 0, 4); fA.add_node(3, 6, 4)
fA.add_section('s', E=E, I=I, A=A)
fA.add_member(0, 0, 2, 's'); fA.add_member(1, 2, 3, 's', release_j=True); fA.add_member(2, 3, 1, 's')
fA.fix(0); fA.fix(1)
fA.member_point_load(1, a=0.0, m=10.0)
rA = solve(fA)

fB = Frame2D()
fB.add_node(0, 0, 0); fB.add_node(1, 6, 0); fB.add_node(2, 0, 4); fB.add_node(3, 6, 4)
fB.add_section('s', E=E, I=I, A=A)
fB.add_member(0, 0, 2, 's'); fB.add_member(1, 2, 3, 's', release_j=True); fB.add_member(2, 3, 1, 's')
fB.fix(0); fB.fix(1)
fB.point_load(2, m=10.0)
rB = solve(fB)

for n in [0, 1]:
    ux, uy, rot = fA.dofs_of(n)
    assert abs(rA.reactions[ux] - rB.reactions[ux]) < 1e-9
    assert abs(rA.reactions[uy] - rB.reactions[uy]) < 1e-9
    assert abs(rA.reactions[rot] - rB.reactions[rot]) < 1e-9
print("  PASS: 桿件a=0的力矩 跟 真正節點力矩(point_load的m), 現在完全等價(修正前是相反號)\n")


print("=== 靜力凝縮版本 solve_condensation() (參考實作) ===")
try:
    solve_condensation(f)
    raise AssertionError("預期應該報錯(release桿件上有內部點載重+局部段均佈載重), 但沒有")
except ValueError as e:
    print(f"  確認: 報錯(已知限制, 符合預期): {str(e)[:50]}...")
print("  回答使用者的問題「凝縮模型和DOFManager都可以嗎?」:")
print("  DOFManager(主要求解器solve())可以直接處理這個組合。")
print("  靜力凝縮(solve_condensation(), 參考實作)不行, 因為release桿件上同時有")
print("  桿件內部點載重跟局部段均佈載重, 這兩種載重都還沒推導release專屬公式。\n")

print("PASS: phase-04-internal-pin_3 完整組合驗證完成")
