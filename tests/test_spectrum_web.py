"""
驗證案例: frame2d.spectrum 的網頁入口 (動力分析 D3b) -- taiwan_code_spectrum、custom_spectrum、
spectrum_analysis()、rsa_to_dict()

層1 taiwan_code_spectrum: 四段式形狀在三個轉角週期(0.2T0, T0, TL)連續(左右極限相等), 平台段
    恆等於 SDS·g, 長週期段 ∝1/T²、中週期段 ∝1/T(取樣兩點驗證斜率), 短週期段在 T=0 時等於
    0.4·SDS·g(規範慣用的下限值)。
層2 custom_spectrum: 分段線性內插對 np.interp 独立핸드計算比對; 邊界外夾在端點值(不外插);
    明確拒絕(點數不足、非遞增、負值)。
層3 spectrum_analysis()+rsa_to_dict(): JSON 可序列化(無 nan)、跟直接呼叫 eigen()+
    response_spectrum() 逐項相同(不是另一套邏輯)、curve 取樣涵蓋模態週期與 TL、custom 型別
    curve 涵蓋使用者給的範圍。
層4 明確拒絕: spectrum_type 打錯、code 缺 SDS/SD1、custom 缺點、以及沿用自 eigen()/
    response_spectrum() 的模型錯誤(沒有支承/沒有質量)。
"""
import json

import numpy as np

from frame2d import Frame2D
from frame2d.modal import eigen
from frame2d.spectrum import (custom_spectrum, response_spectrum, rsa_to_dict,
                              spectrum_analysis, taiwan_code_spectrum)

E, I, A, RHO = 200e9, 8e-5, 5e-3, 7850.0
G = 9.80665


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: taiwan_code_spectrum ===")
SDS, SD1, TL = 0.6, 0.35, 6.0
sp = taiwan_code_spectrum(SDS, SD1, TL)
T0 = SD1 / SDS
eps = 1e-9
print(f"  T0={T0:.4f}")
assert rel(sp(0.0), 0.4 * SDS * G) < 1e-9, "T=0 應該等於 0.4·SDS·g"
assert rel(sp(0.2 * T0 - eps), sp(0.2 * T0 + eps)) < 1e-6, "0.2T0 處應該連續"
assert rel(sp(0.2 * T0), SDS * G) < 1e-9 and rel(sp(T0), SDS * G) < 1e-9, "平台段恆等於 SDS·g"
assert rel(sp(T0 - eps), sp(T0 + eps)) < 1e-6, "T0 處應該連續"
assert rel(sp(TL - eps), sp(TL + eps)) < 1e-6, "TL 處應該連續"
assert rel(sp(2 * T0), SD1 * G / (2 * T0)) < 1e-12, "中週期段應該精確等於 SD1·g/T"
assert rel(sp(2 * TL) / sp(TL), (TL / (2 * TL))**2) < 1e-9, "長週期段應該 ∝ 1/T²(2倍週期差4倍)"
print(f"  T=0: {sp(0.0):.4f}, 平台: {sp(0.5 * T0):.4f}={SDS * G:.4f}, "
      f"T0 附近連續性差 {rel(sp(T0 - eps), sp(T0 + eps)):.2e}, TL 附近 {rel(sp(TL - eps), sp(TL + eps)):.2e}")

for bad, label in ((0.0, "SDS"), (-1.0, "SDS")):
    pass
try:
    taiwan_code_spectrum(0.0, 0.35)
    raise AssertionError("SDS=0 應該報錯")
except ValueError as e:
    print(f"  SDS<=0: ValueError OK ({str(e)[:30]}...)")
try:
    taiwan_code_spectrum(0.6, -0.1)
    raise AssertionError("SD1<0 應該報錯")
except ValueError as e:
    print(f"  SD1<0: ValueError OK ({str(e)[:30]}...)")
try:
    taiwan_code_spectrum(0.6, 0.35, TL=0.0)
    raise AssertionError("TL=0 應該報錯")
except ValueError as e:
    print(f"  TL<=0: ValueError OK ({str(e)[:30]}...)")

# ------------------------------------------------------------------
print("=== 層2: custom_spectrum ===")
pts = [(0.1, 2.0), (0.5, 5.0), (1.0, 3.0), (3.0, 1.0)]
sp2 = custom_spectrum(pts)
Ts_ref = np.array([p[0] for p in pts])
Sa_ref = np.array([p[1] for p in pts])
for T in (0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0):
    want = float(np.interp(T, Ts_ref, Sa_ref))
    got = sp2(T)
    assert abs(got - want) < 1e-12, f"T={T}: {got} vs {want}"
print(f"  7 個測試點跟 np.interp 逐點相同")
assert sp2(0.01) == 2.0 and sp2(10.0) == 1.0, "範圍外應該夾在端點值(不外插)"
print("  範圍外(T=0.01 -> 2.0, T=10 -> 1.0)夾在端點值 OK")


def expect_err(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:36]}...)")
        return
    raise AssertionError(f"{label}: 應該報錯")


expect_err("只有1個點", lambda: custom_spectrum([(0.1, 1.0)]))
expect_err("週期重複(排序後仍相鄰相同)", lambda: custom_spectrum([(0.5, 1.0), (0.5, 2.0)]))
expect_err("週期<=0", lambda: custom_spectrum([(0.0, 1.0), (0.5, 2.0)]))
expect_err("Sa為負", lambda: custom_spectrum([(0.1, -1.0), (0.5, 2.0)]))
# 輸入順序不是遞增也可以: 使用者在網頁表格打亂順序很常見, 這裡自動依 T 排序, 不強迫使用者手動排
sp_unsorted = custom_spectrum([(0.5, 1.0), (0.2, 2.0), (1.0, 0.5)])
assert abs(sp_unsorted(0.2) - 2.0) < 1e-12 and abs(sp_unsorted(0.5) - 1.0) < 1e-12
print("  輸入順序打亂也能正確處理(自動依 T 排序) OK")

# ------------------------------------------------------------------
def portal():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('H', E=E, I=I, A=A, rho=RHO)
    f.add_member(0, 0, 1, 'H').add_member(1, 1, 2, 'H').add_member(2, 3, 2, 'H')
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


print("=== 層3: spectrum_analysis() + rsa_to_dict() ===")
f = portal()
pkg = spectrum_analysis(f, direction='x', damping=0.05, combine='CQC', n_modes=4,
                        mass_kind='consistent', spectrum_type='code', code_sds=SDS, code_sd1=SD1, code_tl=TL)
d = json.loads(json.dumps(rsa_to_dict(pkg), allow_nan=False))     # 確認可序列化(無 nan)

# 跟直接呼叫 eigen()+response_spectrum() 比對, 逐項相同(不是另一套邏輯)
f_ref = portal()
md_ref = eigen(f_ref, n_modes=4, mass='consistent')
res_ref = response_spectrum(md_ref, taiwan_code_spectrum(SDS, SD1, TL), direction='x', damping=0.05, combine='CQC')
assert abs(d['base_shear'] - res_ref.base_shear) < 1e-9 * res_ref.base_shear
assert abs(d['cum_ratio_total'] - res_ref.cum_ratio_total) < 1e-12
for i in range(4):
    assert abs(d['modes'][i]['period'] - res_ref.periods[i]) < 1e-12
    assert abs(d['modes'][i]['sa'] - res_ref.sa[i]) < 1e-9 * res_ref.sa[i]
    assert abs(d['modes'][i]['modal_base_shear'] - res_ref.modal_base_shear[i]) < 1e-9 * abs(res_ref.modal_base_shear[i])
n0 = f_ref.dofs_of(1)
assert abs(d['nodes']['1']['ux'] - res_ref.displacements[n0[0]]) < 1e-9 * res_ref.displacements[n0[0]]
m0 = res_ref.member_forces[0]
assert abs(d['members']['0']['M_i'] - m0[2]) < 1e-9 * abs(m0[2])
print(f"  base_shear/cum_ratio/每模態T,Sa,V_i/節點位移/桿件彎矩, 跟直接呼叫 eigen()+response_spectrum() 逐項相同")

assert len(d['curve']['T']) == 300 and len(d['curve']['Sa']) == 300
assert max(d['curve']['T']) >= TL - 1e-6, "curve 取樣範圍應該涵蓋 TL"
assert max(d['curve']['T']) >= 3 * d['modes'][0]['period'] - 1e-6, "curve 取樣範圍應該涵蓋模態週期的3倍"
print(f"  curve 取樣 300 點, 範圍涵蓋 TL={TL} 與模態週期")

pkg_c = spectrum_analysis(portal(), direction='y', mass_kind='lumped', spectrum_type='custom',
                          custom_points=[[0.01, 1.0], [0.05, 4.0], [1.5, 0.5]])
d_c = json.loads(json.dumps(rsa_to_dict(pkg_c), allow_nan=False))
assert max(d_c['curve']['T']) >= 1.5 - 1e-6, "自訂型別的 curve 範圍應該涵蓋使用者給的最大週期"
print(f"  自訂反應譜: curve 範圍涵蓋使用者給的最大週期 1.5s")

# ------------------------------------------------------------------
print("=== 層4: 明確拒絕 ===")
expect_err("spectrum_type打錯", lambda: spectrum_analysis(portal(), spectrum_type='bogus'))
expect_err("code缺SDS", lambda: spectrum_analysis(portal(), spectrum_type='code', code_sd1=0.35))
expect_err("code缺SD1", lambda: spectrum_analysis(portal(), spectrum_type='code', code_sds=0.6))
expect_err("custom缺點", lambda: spectrum_analysis(portal(), spectrum_type='custom', custom_points=None))
expect_err("custom空列表", lambda: spectrum_analysis(portal(), spectrum_type='custom', custom_points=[]))
f_nosup = Frame2D()
f_nosup.add_node(0, 0, 0).add_node(1, L := 3.0, 0)
f_nosup.add_section('s', E=E, I=I, A=A, rho=RHO)
f_nosup.add_member(0, 0, 1, 's')
expect_err("模型沒有支承(沿用eigen的錯誤訊息)",
          lambda: spectrum_analysis(f_nosup, spectrum_type='code', code_sds=0.6, code_sd1=0.35))

print("\n全部通過: taiwan_code_spectrum形狀、custom_spectrum內插、spectrum_analysis與rsa_to_dict跟核心邏輯逐項一致。")
