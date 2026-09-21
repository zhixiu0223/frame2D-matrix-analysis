"""
驗證案例: 模態形狀曲線 Modal.member_curves() 與網頁用的 modal_to_dict() (動力分析 D2b)

member_curves() 用桿端局部位移 (含release端專屬DOF) + Hermite三次形函數畫出每根桿件變形後的
全域座標, 網頁的模態形狀顯示直接用它。這支測試獨立驗證:

 1. 剛體位移場: 在模態向量裡塞進「x方向剛體平移」與「繞原點剛體轉動」(線性位移場, 剛好在
    形函數空間內), 每個點(含斜桿)必須精確等於解析位置 —— 同時驗證局部/全域座標轉換沒搞反。
 2. 端點: 曲線起訖點 = 節點位置 + 節點位移; 共用節點的兩根桿件在節點處連續。
 3. 精度: 16 元素懸臂梁第1模態, 曲線內部各點跟精確模態形狀 cosh/cos 公式比對。
 4. release: 鉸接端的斜率 = release專屬DOF的轉角(不是節點轉角), 鉸接處兩桿斜率不連續但位置連續。
 5. truss: 維持直線。
 6. modal_to_dict(): 可 json.dumps(allow_nan=False)、放大倍率使最大節點位移 = 對角線×0.15、
    無可動質量的方向比例為 None、曲線點數正確。
"""
import json
import math

import numpy as np

from frame2d import Frame2D
from frame2d.modal import eigen, modal_to_dict

E, I, A, RHO = 200e6, 8e-5, 1e-2, 7.85


def inclined_frame():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4.5).add_node(3, 6.5, 0).add_node(4, 3, 7.5)
    f.add_section('c', E=E, I=I, A=A, rho=RHO)
    f.add_section('t', E=E, I=I, A=1e-3, rho=RHO)
    f.add_member(0, 0, 1, 'c').add_member(1, 1, 2, 'c').add_member(2, 3, 2, 'c').add_member(3, 1, 4, 'c')
    f.add_truss(4, 4, 2, 't')
    f.fix(0).fix(3)
    f.add_mass(1, mx=3.0, my=1.0, Iz=0.5).add_mass(2, mx=2.0, my=2.5).add_mass(4, mx=1.5, my=1.0)
    return f


print("=== 1. 剛體位移場: 座標轉換與Hermite在線性場下精確 ===")
f = inclined_frame()
md = eigen(f, n_modes=3, mass='consistent')
n = md.n_dof
for label, setter in (
    ("x方向剛體平移", lambda nid, nd, v: (v.__setitem__(f.dofs_of(nid)[0], 1.0))),
    ("繞原點剛體轉動", lambda nid, nd, v: (v.__setitem__(f.dofs_of(nid)[0], -nd.y),
                                       v.__setitem__(f.dofs_of(nid)[1], nd.x),
                                       v.__setitem__(f.dofs_of(nid)[2], 1.0))),
):
    vec = np.zeros(n)
    for nid, nd in f.nodes.items():
        setter(nid, nd, vec)
    md.phi[:, 0] = vec
    scale = 0.7
    worst = 0.0
    for mid, (X, Y) in md.member_curves(0, scale=scale, n=11).items():
        m = f.members[mid]
        ni, nj = f.nodes[m.node_i], f.nodes[m.node_j]
        for k, xi in enumerate(np.linspace(0, 1, 11)):
            x0, y0 = ni.x + xi * (nj.x - ni.x), ni.y + xi * (nj.y - ni.y)
            if label.startswith("x方向"):
                ex, ey = x0 + scale * 1.0, y0
            else:
                ex, ey = x0 + scale * (-y0), y0 + scale * x0
            worst = max(worst, abs(X[k] - ex), abs(Y[k] - ey))
    print(f"  {label}: 所有桿件(含斜桿與桁架)各點與解析位置的最大差 = {worst:.2e}")
    assert worst < 1e-12, f"{label}: 曲線與解析位置不符(座標轉換錯了?)"

print("=== 2. 端點 = 節點位置+節點位移, 共用節點連續 ===")
md = eigen(inclined_frame(), n_modes=3, mass='consistent')
f = md.frame
for mode in range(3):
    scale = 2.5
    curves = md.member_curves(mode, scale=scale, n=9)
    node_pos = {}
    for mid, (X, Y) in curves.items():
        m = f.members[mid]
        for end, nid in ((0, m.node_i), (-1, m.node_j)):
            ux, uy, _ = md.node_shape(mode, nid)
            nd = f.nodes[nid]
            assert abs(X[end] - (nd.x + scale * ux)) < 1e-12 and abs(Y[end] - (nd.y + scale * uy)) < 1e-12, \
                f"模態{mode + 1} 桿件{mid} 端點沒有落在節點位移後的位置"
            node_pos.setdefault(nid, []).append((X[end], Y[end]))
    for nid, pts in node_pos.items():
        for p in pts[1:]:
            assert abs(p[0] - pts[0][0]) < 1e-12 and abs(p[1] - pts[0][1]) < 1e-12, f"節點{nid}處桿件曲線不連續"
print("  3個模態: 所有桿件端點位置正確、共用節點處連續")

print("=== 3. 精度: 16元素懸臂梁第1模態 vs 精確模態形狀 ===")
L = 3.0
N_EL = 16
fc = Frame2D()
for i in range(N_EL + 1):
    fc.add_node(i, L * i / N_EL, 0)
fc.add_section('s', E=E, I=I, A=A, rho=RHO)
for i in range(N_EL):
    fc.add_member(i, i, i + 1, 's')
fc.fix(0)
for i in range(1, N_EL + 1):
    fc.support(i, ux=0.0)
mdc = eigen(fc, n_modes=1, mass='consistent')
b = 1.8751040687
sig = (math.cosh(b) + math.cos(b)) / (math.sinh(b) + math.sin(b))
phi_ex = lambda xr: math.cosh(b * xr) - math.cos(b * xr) - sig * (math.sinh(b * xr) - math.sin(b * xr))
tip = mdc.node_shape(0, N_EL)[1]
worst = 0.0
for mid, (X, Y) in mdc.member_curves(0, scale=1.0, n=7).items():
    for x, y in zip(X, Y):
        worst = max(worst, abs(y / tip - phi_ex(x / L) / phi_ex(1.0)))
print(f"  曲線各點與精確模態形狀(尖端歸一)的最大差 = {worst:.2e}")
assert worst < 1e-4, "懸臂梁第1模態曲線與精確形狀差太多"

print("=== 4. release: 鉸接端斜率 = 專屬DOF轉角 ===")
fr = Frame2D()
fr.add_node(0, 0, 0).add_node(1, 3, 0).add_node(2, 6, 0)
fr.add_section('s', E=E, I=I, A=A, rho=RHO)
fr.add_member(0, 0, 1, 's', release_j=True).add_member(1, 1, 2, 's')
fr.fix(0)
fr.support(1, ux=0.0).support(2, ux=0.0).support(2, uy=0.0)      # 節點2: 只擋垂直(uy), 轉角自由
mdr = eigen(fr, n_modes=2, mass='consistent')
extra_dof = mdr.member_dofs[0][5]
assert extra_dof >= mdr.n_node_dof, "release端的第6個DOF應該是專屬(額外)DOF"
for mode in range(2):
    curves = mdr.member_curves(mode, scale=1.0, n=4001)
    X0, Y0 = curves[0]
    X1, Y1 = curves[1]
    slope_end0 = (Y0[-1] - Y0[-2]) / (X0[-1] - X0[-2])          # 桿件0(鉸接端)在節點1的斜率
    slope_start1 = (Y1[1] - Y1[0]) / (X1[1] - X1[0])            # 桿件1在節點1的斜率
    th_extra = mdr.phi[extra_dof, mode]
    th_node = mdr.node_shape(mode, 1)[2]
    print(f"  模態{mode + 1}: 桿0端斜率 {slope_end0:+.6f} (專屬DOF轉角 {th_extra:+.6f}), "
          f"桿1起點斜率 {slope_start1:+.6f} (節點轉角 {th_node:+.6f})")
    assert abs(slope_end0 - th_extra) < 5e-3 * max(abs(th_extra), 1e-3), "鉸接端斜率應該等於專屬DOF轉角"
    assert abs(slope_start1 - th_node) < 5e-3 * max(abs(th_node), 1e-3), "桿1起點斜率應該等於節點轉角"
    assert abs(X0[-1] - X1[0]) < 1e-12 and abs(Y0[-1] - Y1[0]) < 1e-12, "鉸接處位置必須連續"
    assert abs(th_extra - th_node) > 1e-6, "這個案例鉸接兩側的轉角應該不同(否則沒測到release)"

print("=== 5. truss 維持直線 ===")
md5 = eigen(inclined_frame(), n_modes=2, mass='consistent')
X, Y = md5.member_curves(0, scale=3.0, n=15)[4]
d = np.hypot(X[-1] - X[0], Y[-1] - Y[0])
dev = max(abs((X[k] - X[0]) * (Y[-1] - Y[0]) - (Y[k] - Y[0]) * (X[-1] - X[0])) / d for k in range(15))
print(f"  桁架桿件各點偏離端點連線的最大距離 = {dev:.2e}")
assert dev < 1e-12

print("=== 6. modal_to_dict() ===")
md6 = eigen(fc, n_modes=3, mass='consistent')
out = modal_to_dict(md6, n_curve=11)
json.dumps(out, allow_nan=False)                      # 不能有 nan/inf
assert out["analysis_type"] == "modal" and out["n_modes"] == 3 and len(out["modes"]) == 3
xs = [nd.x for nd in fc.nodes.values()]
ys = [nd.y for nd in fc.nodes.values()]
diag = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
for k, mode in enumerate(out["modes"]):
    assert mode["index"] == k + 1 and abs(mode["omega"] - md6.omega[k]) < 1e-12
    assert abs(mode["period"] * mode["frequency"] - 1.0) < 1e-12
    assert all(len(c["X"]) == 11 and len(c["Y"]) == 11 for c in mode["curves"].values())
    # 放大倍率: 最大節點位移放大後 = 對角線 * 0.15
    max_d = max(math.hypot(md6.phi[fc.dofs_of(nid)[0], k], md6.phi[fc.dofs_of(nid)[1], k]) for nid in fc.nodes)
    assert abs(mode["deform_scale"] * max_d - diag * 0.15) < 1e-9 * diag
    assert mode["ratio_x"] is None and mode["cum_x"] is None, "x方向被全部拘束(無可動質量), 比例應該是None"
    assert mode["ratio_y"] is not None and 0 <= mode["ratio_y"] <= 1
    assert abs(mode["gamma_y"]) == abs(md6.gamma['y'][k]) or abs(mode["gamma_y"] - md6.gamma['y'][k]) < 1e-12
assert abs(out["mass_total"]["y"] - md6.mass_total["y"]) < 1e-12
print("  JSON 可序列化(無nan)、放大倍率、無可動質量方向為None、點數 OK")

print("\n全部通過: 模態形狀曲線(座標轉換、端點連續、精度、release、桁架)與 JSON 序列化都吻合。")
