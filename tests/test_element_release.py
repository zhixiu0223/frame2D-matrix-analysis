"""
驗證案例15: 內部鉸接(internal pin / element release), Load System Phase 4

用靜力凝縮(static condensation)處理桿端鉸接, 不改變DOF系統(dofs_of()
完全不用動)。這是這個功能不需要先做完整DOFManager升級就能實作的原因——
只需要修改單一桿件的局部勁度矩陣, 詳見elements.py的member_stiffness_local()
說明。

驗證策略: 用古典Gerber梁(兩跨連續梁, 中間支承處放內部鉸接)當最小驗證案例。
沒有鉸接時, 兩跨連續梁是一次靜不定; 加了鉸接後變成靜定結構, 可以直接用
純靜力學驗證(拆解成兩個獨立簡支梁分析), 不依賴frame2d自己的公式對不對。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01


def check(label, fem, exact, atol=1e-6):
    err = abs(fem - exact)
    print(f"  {label}: FEM={fem:.6f}  預期={exact:.6f}  誤差={err:.2e}")
    assert err < atol, f"{label} 不吻合"


# ---- 案例A: Gerber梁, UDL只加在第一跨(有鉸接那一跨) ----
print("=== 案例A: Gerber梁(release_j), UDL在第一跨 ===")
L1, L2, w = 6.0, 5.0, -10.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L1, 0)
f.add_node(2, L1 + L2, 0)
f.add_section('s', E=E, I=I, A=A)
f.add_member(0, 0, 1, 's', release_j=True)   # 鉸接在node1(第一跨的J端)
f.add_member(1, 1, 2, 's')
f.pin(0)
f.roller_y(1)
f.roller_y(2)
f.distributed_load(0, w=w)
r = solve(f)

# 純靜力學(拆成兩個獨立簡支梁): member0單獨簡支梁(pin-pin,兩端M=0),
# R0=R1=|w|*L1/2; member1完全無載重, 不貢獻任何內力
R_exact = abs(w) * L1 / 2
check("R0", r.reactions[f.dofs_of(0)[1]], R_exact)
check("R1", r.reactions[f.dofs_of(1)[1]], R_exact)
check("R2", r.reactions[f.dofs_of(2)[1]], 0.0)
ef0 = r.member_results[0].end_forces_local
ef1 = r.member_results[1].end_forces_local
check("鉸接端M2(member0)", ef0[5], 0.0)
for i, v in enumerate(ef1):
    assert abs(v) < 1e-6, f"member1(無載重段)應該完全無內力, 但index{i}={v}"
print("PASS: Gerber梁完全吻合純靜力學(拆解成獨立簡支梁分析)\n")


# ---- 案例B: release_i 的UDL固定端反力公式, 獨立用簡支梁純靜力學驗證 ----
print("=== 案例B: release_i的UDL公式, 獨立簡支梁驗證 ===")
L3, w3 = 6.0, -10.0
fB = Frame2D()
fB.add_node(0, 0, 0)
fB.add_node(1, L3, 0)
fB.add_section('s', E=E, I=I, A=A)
fB.add_member(0, 0, 1, 's', release_i=True)
fB.roller_y(0)
fB.pin(1)
fB.distributed_load(0, w=w3)
rB = solve(fB)
R_exact_B = abs(w3) * L3 / 2
check("R0", rB.reactions[fB.dofs_of(0)[1]], R_exact_B)
check("R1", rB.reactions[fB.dofs_of(1)[1]], R_exact_B)
ef_B = rB.member_results[0].end_forces_local
check("鉸接端M1(release_i)", ef_B[2], 0.0)
print("PASS: release_i的UDL公式獨立驗證通過(簡支梁等效, 純靜力學)\n")


# ---- 案例C: 同一個Gerber梁, 改成用"對面那根桿件"的release_i做同一個鉸接
#      (而不是member0的release_j), 物理上應該給出完全相同的結果 ----
print("=== 案例C: 同一個鉸接改用另一根桿件的release_i, 交叉驗證 ===")
fC = Frame2D()
fC.add_node(0, 0, 0)
fC.add_node(1, L1, 0)
fC.add_node(2, L1 + L2, 0)
fC.add_section('s', E=E, I=I, A=A)
fC.add_member(0, 0, 1, 's')                     # 這次不釋放
fC.add_member(1, 1, 2, 's', release_i=True)      # 改成這根釋放I端
fC.pin(0)
fC.roller_y(1)
fC.roller_y(2)
fC.distributed_load(0, w=w)   # 載重一樣加在第一跨(這次是剛接的那根)
rC = solve(fC)
check("R0(改用另一根桿件釋放)", rC.reactions[fC.dofs_of(0)[1]], R_exact)
check("R1(改用另一根桿件釋放)", rC.reactions[fC.dofs_of(1)[1]], R_exact)
check("R2(改用另一根桿件釋放)", rC.reactions[fC.dofs_of(2)[1]], 0.0)
print("PASS: 同一個鉸接不管哪根桿件負責釋放, 物理結果完全一致\n")

print("PASS: Phase 4(內部鉸接) 最小模型驗證完成")
