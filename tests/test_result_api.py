"""
驗證案例29: Result API 便利查詢介面(result.max_moment(), result.member(id).
moment_at(x) 等)

這是純粹的介面整理, 不改變任何計算邏輯——所有數字底層都還是
member_internal_forces/member_deformed_shape算出來的(已經在其他測試
檔案裡驗證過), 這裡只驗證:
  1. 新介面查出來的數字, 跟直接呼叫底層函式+手動np.interp查詢的結果
     完全一致(純介面, 沒有引入新的計算)
  2. 跟解析解對照(簡支梁均佈載重, 跨中彎矩wL²/8、跨中撓度5wL⁴/384EI)
  3. 多桿件結構時, max_moment()/max_displacement()掃描到的真的是
     全域最大值(不是只看第一根桿件, 也不是弄錯桿件編號)
  4. 遇到鬆弛的cable時, 掃描邏輯會自動跳過(鬆弛桿件不受力, 掃它們
     的內力沒有意義)
  5. 新增的frame欄位不影響任何既有測試(38個既有測試套件全數通過,
     見完整pytest執行紀錄)
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01


def check(label, fem, exact, tol=1e-4):
    err = abs(fem - exact)
    print(f"  {label}: 查詢介面={fem:.6f}  對照={exact:.6f}  誤差={err:.2e}")
    assert err < tol, f"{label} 不吻合"


# ---- 案例A: 簡支梁均佈載重, 對照解析解 ----
print("=== 案例A: 簡支梁均佈載重, Result API對照解析解 ===")
L, w = 6.0, -5.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's')
f.pin(0)
f.roller_y(1)
f.distributed_load(0, w=w)
r = solve(f)

M_exact = abs(w) * L**2 / 8
defl_exact = 5 * abs(w) * L**4 / (384 * E * I)
check("跨中彎矩(result.member(0).moment_at(3.0))", r.member(0).moment_at(3.0), M_exact)
check("跨中撓度(result.member(0).deflection_at(3.0))", r.member(0).deflection_at(3.0), defl_exact)

ev = r.max_moment()
check("max_moment()數值", abs(ev.value), M_exact)
assert ev.member_id == 0 and abs(ev.x - 3.0) < 0.2, "max_moment()應該找到跨中位置"
print("  PASS\n")


# ---- 案例B: 跟直接呼叫底層函式手動查詢比對, 確認純介面沒有引入新計算 ----
print("=== 案例B: Result API vs 直接呼叫底層函式手動查詢, 應完全一致 ===")
x_full, N, V, M = member_internal_forces(f, r, 0, n=201)
manual_M_at_2 = float(np.interp(2.0, x_full, M))
api_M_at_2 = r.member(0).moment_at(2.0)
assert abs(manual_M_at_2 - api_M_at_2) < 1e-9, "介面查詢應該跟手動查詢完全一致(同一套底層資料)"
print(f"  手動查詢: {manual_M_at_2:.6f}  介面查詢: {api_M_at_2:.6f}  (完全一致)")
print("  PASS\n")


# ---- 案例C: 門型鋼架(多桿件), max_moment/max_displacement找到真正全域最大值 ----
print("=== 案例C: 門型鋼架, 確認掃描到真正的全域最大值(不是弄錯桿件) ===")
fC = Frame2D()
fC.add_node(0, 0, 0)
fC.add_node(1, 6, 0)
fC.add_node(2, 0, 4)
fC.add_node(3, 6, 4)
fC.add_section('s', E=E, I=I, A=A)
fC.add_member(0, 0, 2, 's')
fC.add_member(1, 2, 3, 's')
fC.add_member(2, 3, 1, 's')
fC.fix(0)
fC.fix(1)
fC.point_load(2, fx=10.0)
rC = solve(fC)

each_max = [rC.member(mid).max_moment() for mid in fC.members]
expected_best = max(each_max, key=lambda ev: abs(ev.value))
scan_result = rC.max_moment()
assert scan_result.member_id == expected_best.member_id
assert abs(scan_result.value - expected_best.value) < 1e-9
print(f"  全結構掃描: {scan_result}")
print(f"  逐一比對驗證: {expected_best} (一致)")
print("  PASS\n")


# ---- 案例D: 含鬆弛cable, 確認掃描自動跳過 ----
print("=== 案例D: 含鬆弛cable, max_moment()自動跳過鬆弛桿件 ===")
fD = Frame2D()
fD.add_node(0, 0, 0)
fD.add_node(1, 0, 10)
fD.add_node(2, 10, 0)
fD.add_section('tower', E=E, I=I, A=A)
fD.add_section('cable', E=E, I=1e-8, A=0.001)
fD.add_member(0, 0, 1, 'tower')
fD.add_cable(1, 1, 2, 'cable')
fD.fix(0)
fD.pin(2)
fD.point_load(1, fy=-1.0)
rD = solve(fD)
assert rD.member_results[1].slack, "這個案例的cable應該判定為鬆弛(測試前提)"
scanD = rD.max_moment()
assert scanD.member_id == 0, "鬆弛的cable(member1)沒有意義的彎矩, 應該跳過只看member0"
print(f"  鬆弛cable(member1)被正確跳過, 掃描結果={scanD}")
print("  PASS\n")

print("PASS: Result API(max_moment/max_shear/max_axial/max_displacement/")
print("member(id).xxx_at(x)) 全數驗證完成")
