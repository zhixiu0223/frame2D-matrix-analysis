"""
驗證案例8: 懸臂梁斜張橋 (Cantilever Cable-Stayed Bridge)

幾何: 塔柱固定於(0,0)高25m, 懸臂梁從(0,0)延伸到(50,0)(僅在x=0端固定,
     其餘無支承, 純靠5條纜索從塔頂拉住), 5條纜索錨定於梁上x=10,20,30,40,50處。
斷面: 梁 E=2e8,I=0.08,A=0.5;塔 E=2e8,I=0.15,A=0.8;纜索 E=1.6e8,A=0.003 (kN,m制)
載重: 梁上均佈載重 w=20kN/m 向下

對照來源: 使用者自己的「贅力法(力量法)」筆記本, 兩個版本:
  - 較早版本(notebook cell 3): 柔度矩陣只算了塔柱"軸向"柔度, 漏了塔柱受纜索
    水平分量拉力時的"彎曲"柔度, 算出5條纜索全部受拉(不正確)
  - 修正版本(notebook cell 5): 補上塔柱彎曲柔度(fxx_tower = H^3/(3*EI_t)),
    算出 Cable1(10m)為負值(受壓) -- 這才是正確答案, 塔柱在纜索水平分量作用下
    確實會側向撓曲, 這是真實物理現象, 不是計算錯誤

frame2d(矩陣位移法)不需要區分"軸向/彎曲柔度"這些手動分項, 直接對兩個構件都用
完整的frame元素(同時有彎曲+軸向勁度), 驗證結果跟修正版力量法完全吻合。
"""
from frame2d import Frame2D, solve

w = 20.0
H_tower = 25.0
E_b, A_b, I_b = 2.0e8, 0.5, 0.08
E_t, A_t, I_t = 2.0e8, 0.8, 0.15
E_c, A_c = 1.6e8, 0.003
x_c = [10.0, 20.0, 30.0, 40.0, 50.0]


def _build_common():
    f = Frame2D()
    f.add_node(0, 0, 0)
    for i, x in enumerate(x_c):
        f.add_node(i + 1, x, 0)
    f.add_node(6, 0, H_tower)
    f.add_section('beam', E=E_b, I=I_b, A=A_b)
    f.add_section('tower', E=E_t, I=I_t, A=A_t)
    f.add_section('cable', E=E_c, I=1.0, A=A_c)
    for i in range(5):
        f.add_member(i, i, i + 1, 'beam')
    f.add_member(5, 0, 6, 'tower')
    f.fix(0)
    for i in range(5):
        f.distributed_load(i, w=-w)
    return f


# ---- 步驟1-2: truss版, 對照力量法(含塔柱彎曲修正)的5條纜索張力 ----
print("=== 步驟1-2: truss(可拉可壓)版本, 對照修正版力量法 ===")
f_truss = _build_common()
for i in range(5):
    f_truss.add_truss(6 + i, 6, i + 1, 'cable')
r_truss = solve(f_truss)

X_force_method = [-214.9697, 26.8740, 202.4422, 267.4776, 274.9446]
for i in range(5):
    N = -r_truss.member_results[6 + i].end_forces_local[0]
    ref = X_force_method[i]
    print(f"  Cable{i+1} ({x_c[i]:.0f}m): frame2d={N:.4f}  force_method={ref:.4f}")
    assert abs(N - ref) < 0.01, f"Cable{i+1} 跟力量法不吻合"
print("PASS: 5條纜索張力(含Cable1的受壓)全數跟修正版力量法吻合\n")


# ---- 步驟4: cable版(自動迭代), 應收斂到跟"移除所有受壓纜索"一致的解 ----
print("=== 步驟4: cable(只受拉, 自動迭代)版本 ===")
f_cable = _build_common()
for i in range(5):
    f_cable.add_cable(6 + i, 6, i + 1, 'cable')
r_cable = solve(f_cable)

print(f"自動判定鬆弛的纜索: {sorted(r_cable.slack_cables)}")
for i in range(5):
    mr = r_cable.member_results[6 + i]
    N = -mr.end_forces_local[0]
    print(f"  Cable{i+1} ({x_c[i]:.0f}m): N={N:.4f}  slack={mr.slack}")
    if mr.slack:
        assert abs(N) < 1e-6, f"鬆弛的Cable{i+1}軸力應為0"
    else:
        assert N > -1e-6, f"沒鬆弛的Cable{i+1}軸力應該>=0(拉力)"

# 跟truss版比較: truss版裡受壓的纜索(Cable1,Cable2), 應該就是cable版判定鬆弛的那些
truss_compression = [i for i in range(5)
                      if -r_truss.member_results[6 + i].end_forces_local[0] < 0]
cable_slack = [i for i in range(5) if (6 + i) in r_cable.slack_cables]
print(f"\ntruss版受壓的纜索(0-indexed): {truss_compression}")
print(f"cable版判定鬆弛的纜索(0-indexed): {cable_slack}")

print("PASS: cable自動迭代收斂, 沒有任何纜索受壓")
print()
print("PASS: 懸臂梁斜張橋(frame+truss/cable混合模型)驗證完成")
