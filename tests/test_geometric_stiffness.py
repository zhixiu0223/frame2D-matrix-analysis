"""
驗證案例: 幾何剛度矩陣 local_geometric_stiffness() / P-Delta 疊加

案例A: P=0 完全退化 -- 疊加P=0時member_stiffness_global()結果跟原本
       (不傳P參數)完全一致, 鎖住"新增這個功能不影響任何既有結果"這件事。
案例B: 獨立推導驗證 -- 不信任移植過來的公式本身, 用一致剛度矩陣的定義
       (P * ∫ N_i'(x) N_j'(x) dx, N_i為Hermite三次形狀函數)直接用sympy
       從頭積分算一次, 跟local_geometric_stiffness()的輸出比對。
案例C: 物理合理性 -- 拉力(P>0)應該增加側向勁度, 壓力(P<0)應該降低側向
       勁度(這是P-Delta效應的物理本質: 壓力讓結構更容易側移)。
案例D: 帶內部鉸(release)時傳入非零P必須直接raise, 不能默默算出
       未驗證過的結果。
"""
import numpy as np
import sympy as sp
from frame2d import Frame2D, solve
from frame2d.elements import (
    local_geometric_stiffness,
    member_stiffness_local,
    member_stiffness_global,
)

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0

# ---- 案例A: P=0 完全退化 ----
print("=== 案例A: P=0 退化回原本結果 ===")
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, L, 0)
f.add_section('sec', E=E, I=I, A=A)
f.add_member(0, node_i=0, node_j=1, section='sec')

sec = f.sections['sec']
node_i, node_j = f.nodes[0], f.nodes[1]

k_no_P, *_ = member_stiffness_global(sec, node_i, node_j)
k_P_zero, *_ = member_stiffness_global(sec, node_i, node_j, P=0.0)
assert np.allclose(k_no_P, k_P_zero, atol=0), "P=0.0應與完全不傳P的結果逐位元一致"
print("PASS: P=0.0與不傳P參數的結果逐位元相同\n")

# ---- 案例B: 用sympy從Hermite形狀函數獨立積分驗證公式本身 ----
print("=== 案例B: 獨立sympy積分驗證公式 ===")
x, Ls, Ps = sp.symbols('x L P', positive=True)
N1 = 1 - 3 * (x / Ls) ** 2 + 2 * (x / Ls) ** 3
N2 = x - 2 * x ** 2 / Ls + x ** 3 / Ls ** 2
N3 = 3 * (x / Ls) ** 2 - 2 * (x / Ls) ** 3
N4 = -x ** 2 / Ls + x ** 3 / Ls ** 2
Ns = [N1, N2, N3, N4]
dNs = [sp.diff(N, x) for N in Ns]

P_test, L_test = 37.5, 4.0  # 任意取值, 跟上面案例的L不必相同
Kg_sympy = np.zeros((4, 4))
for i in range(4):
    for j in range(4):
        expr = Ps * sp.integrate(dNs[i] * dNs[j], (x, 0, Ls))
        Kg_sympy[i, j] = float(expr.subs({Ps: P_test, Ls: L_test}))

kg_full = local_geometric_stiffness(P_test, L_test)
bend_idx = [1, 2, 4, 5]
Kg_from_code = kg_full[np.ix_(bend_idx, bend_idx)]

assert np.allclose(Kg_from_code, Kg_sympy, rtol=1e-10), (
    f"local_geometric_stiffness()跟sympy從Hermite形狀函數獨立積分的結果不符\n"
    f"code=\n{Kg_from_code}\nsympy=\n{Kg_sympy}"
)
# 軸向自由度(索引0,3)不應受幾何剛度影響
assert np.allclose(kg_full[[0, 3]][:, [0, 3]], 0.0), "幾何剛度矩陣不應影響軸向自由度"
print("PASS: 公式跟獨立sympy積分結果逐項吻合 (非照抄, 是重新積分驗證)\n")

# ---- 案例C: 物理合理性 -- 拉力增勁度、壓力減勁度 ----
print("=== 案例C: 拉力增勁度/壓力減勁度 (P-Delta物理方向) ===")
k_elastic = member_stiffness_local(E, I, A, L)
k_tension = member_stiffness_local(E, I, A, L, P=+1000.0)
k_compression = member_stiffness_local(E, I, A, L, P=-1000.0)

# 側向平移自由度v1 (索引1) 的自身勁度項
kv1_elastic = k_elastic[1, 1]
kv1_tension = k_tension[1, 1]
kv1_compression = k_compression[1, 1]

assert kv1_tension > kv1_elastic > kv1_compression, (
    f"側向勁度應滿足 拉力({kv1_tension}) > 純彈性({kv1_elastic}) > 壓力({kv1_compression})"
)
print(f"PASS: k_v1v1 拉力={kv1_tension:.2f} > 純彈性={kv1_elastic:.2f} > 壓力={kv1_compression:.2f}\n")

# ---- 案例D: 帶release時傳入非零P必須raise, 不能默默算錯 ----
print("=== 案例D: release + 非零P 必須raise NotImplementedError ===")
raised = False
try:
    member_stiffness_local(E, I, A, L, release_j=True, P=500.0)
except NotImplementedError:
    raised = True
assert raised, "release端+非零P應該raise NotImplementedError, 而不是默默算出未驗證的結果"
print("PASS: release+P同時使用正確地raise\n")

print("PASS: 幾何剛度矩陣(P-Delta)所有案例通過")
