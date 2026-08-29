"""
驗證案例16: 內部鉸接(Phase 4) 對照 SW FEA 第三方工具

使用者上傳的 phase-04-internal-pin: 固接-固接門型鋼架, 節點0=(0,0)、
節點1=(6,0)為基礎, 節點2=(0,4)、節點3=(6,4)為柱頂。F0=左柱(剛接)、
F1=梁、F2=右柱(剛接), F1跟F2在節點3(右上角)是內部鉸接(不傳彎矩),
節點2(左上角)是一般剛接。節點2受水平力10kN向右(180°=向右, 已在
test_case08_vs_swfea.py確認過的角度慣例)。

斷面 E=200GPa, I=8.0E7mm^4(=8e-5m^4), A=0.01m^2, 跟本專案kN,m慣例
直接對應。
"""
import numpy as np
from frame2d import Frame2D, solve
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01


def _build():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')                    # F0: 左柱, 剛接
    f.add_member(1, 2, 3, 's', release_j=True)     # F1: 梁, J端(node3)鉸接
    f.add_member(2, 3, 1, 's')                     # F2: 右柱, 剛接
    f.fix(0)
    f.fix(1)
    f.point_load(2, fx=10.0)   # 10kN向右
    return f


f = _build()
r = solve(f)

print("=== 反力比對 ===")
sw_reactions = {0: (-6.671, -1.482, 17.792), 1: (-3.329, 1.482, 13.316)}
for nid, (Rx, Ry, M) in sw_reactions.items():
    ux, uy, rot = f.dofs_of(nid)
    print(f"  node{nid}: Rx={r.reactions[ux]:.3f}(sw={Rx})  Ry={r.reactions[uy]:.3f}(sw={Ry})  M={r.reactions[rot]:.3f}(sw={M})")
    assert abs(r.reactions[ux] - Rx) < 0.01
    assert abs(r.reactions[uy] - Ry) < 0.01
    assert abs(r.reactions[rot] - M) < 0.01
print("PASS\n")

print("=== 逐點BM比對 (每根桿件11個位置點) ===")
sw_members_BM = {
    0: [-17.792, -15.124, -12.455, -9.787, -7.118, -4.450, -1.782, 0.887, 3.555, 6.224, 8.892],
    1: [8.892, 8.003, 7.114, 6.224, 5.335, 4.446, 3.557, 2.668, 1.778, 0.889, 0.000],
    2: [-0.000, 1.332, 2.663, 3.995, 5.326, 6.658, 7.990, 9.321, 10.653, 11.984, 13.316],
}
max_err = 0
for mid, bm_list in sw_members_BM.items():
    x, N, V, M = member_internal_forces(f, r, mid, n=11)
    for i, bmsw in enumerate(bm_list):
        max_err = max(max_err, abs(M[i] - bmsw))
    print(f"  member{mid}: 最大誤差={max(abs(M[i]-bmsw) for i,bmsw in enumerate(bm_list)):.5f}")
assert max_err < 0.01, f"BM最大誤差過大: {max_err}"
print(f"整體最大誤差: {max_err:.5f}")
print("PASS\n")

# 鉸接端(F1的J端, node3)彎矩必須精確為0 -- 這是鉸接的定義
ef1 = r.member_results[1].end_forces_local
assert abs(ef1[5]) < 1e-9, f"鉸接端彎矩應為0, 得到{ef1[5]}"
print(f"鉸接端(F1的J端)彎矩: {ef1[5]:.2e} (精確為0, 符合鉸接定義)")

print("\nPASS: Phase 4(內部鉸接) 對照SW FEA門型鋼架案例, 驗證完成")
