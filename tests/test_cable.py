"""
驗證案例7: 纜線(cable)元素 -- 只能受拉, 受壓時自動鬆弛退出作用並重新求解

案例A: 跟test_truss.py案例B完全一樣的幾何(A字撐架), 但兩根桿都改成cable。
       因為兩根都應該受壓(撐架, 不是懸吊), 兩條cable都會鬆弛退出作用,
       結構變成機構(無法承受這個垂直載重), 應該要拋出RuntimeError。

案例B: 反過來的情境 -- 兩條纜線從高處"懸吊"住一個節點(節點在下面, 錨點在
       上面), 這種配置下纜線應該受拉, 不會鬆弛, 結果應該跟truss版本完全一樣。

案例C: 三條纜線懸吊一個節點, 其中垂直載重的方向讓某一條纜線"太鬆"
       (幾何上不需要它出力甚至會被壓), 驗證只有那一條會被移除、
       另外兩條正常受拉、結果自洽(仍然滿足力平衡)。
"""
import numpy as np
from frame2d import Frame2D, solve

# ---- 案例A: 撐架幾何但用cable -> 應該兩條都鬆弛, 變成機構, 應該報錯 ----
print("=== 案例A: 撐架幾何(該受壓)但用cable, 應該偵測到機構並報錯 ===")
E, A, P = 200e6, 0.01, 10.0
f = Frame2D()
f.add_node(0, 0, 0)
f.add_node(1, 4, 0)
f.add_node(2, 2, 3)
f.add_section('s', E=E, I=1.0, A=A)
f.add_cable(0, 0, 2, 's')
f.add_cable(1, 1, 2, 's')
f.pin(0)
f.pin(1)
f.point_load(2, fy=-P)
try:
    r = solve(f)
    raise AssertionError("預期應該要報錯(兩條cable都會鬆弛變成機構), 但沒有報錯")
except (RuntimeError, ValueError) as e:
    print(f"PASS: 正確偵測到機構並報錯: {e}\n")


# ---- 案例B: 懸吊幾何(該受拉), 用cable應該跟truss結果完全一樣 ----
print("=== 案例B: 懸吊幾何(該受拉), cable結果應該跟truss完全一樣 ===")
# 節點0=(0,5)、節點1=(4,5) 是高處錨點, 節點2=(2,2) 是被吊住的節點, 垂直載重向下
f_truss = Frame2D()
f_truss.add_node(0, 0, 5); f_truss.add_node(1, 4, 5); f_truss.add_node(2, 2, 2)
f_truss.add_section('s', E=E, I=1.0, A=A)
f_truss.add_truss(0, 0, 2, 's'); f_truss.add_truss(1, 1, 2, 's')
f_truss.pin(0); f_truss.pin(1)
f_truss.point_load(2, fy=-P)
r_truss = solve(f_truss)

f_cable = Frame2D()
f_cable.add_node(0, 0, 5); f_cable.add_node(1, 4, 5); f_cable.add_node(2, 2, 2)
f_cable.add_section('s', E=E, I=1.0, A=A)
f_cable.add_cable(0, 0, 2, 's'); f_cable.add_cable(1, 1, 2, 's')
f_cable.pin(0); f_cable.pin(1)
f_cable.point_load(2, fy=-P)
r_cable = solve(f_cable)

for mid in [0, 1]:
    N_truss = -r_truss.member_results[mid].end_forces_local[0]
    N_cable = -r_cable.member_results[mid].end_forces_local[0]
    print(f"桿{mid}軸力: truss={N_truss:.6f}  cable={N_cable:.6f}")
    assert N_truss > 0, "懸吊幾何應該受拉(truss結果)"
    assert abs(N_truss - N_cable) < 1e-8, "受拉情況下cable應該跟truss結果完全一樣"
    assert not r_cable.member_results[mid].slack, "受拉的cable不應該被判定為鬆弛"
print("PASS: 懸吊(受拉)情況下, cable跟truss結果完全一致, 沒有任何cable鬆弛\n")


# ---- 案例C: 三條纜線, 其中一條在這個載重下明顯多餘, 應該被自動移除 ----
print("=== 案例C: 三條纜線懸吊, 其中一條方向不對, 該被自動移除 ===")
f3 = Frame2D()
f3.add_node(0, -3, 5)   # 左上錨點
f3.add_node(1, 3, 5)    # 右上錨點
f3.add_node(2, 0, 0)    # 正下方錨點(這條纜線只能往下拉, 對抵抗向下載重沒有幫助)
f3.add_node(3, 0, 2)    # 被吊住的節點
f3.add_section('s', E=E, I=1.0, A=A)
f3.add_cable(0, 0, 3, 's')
f3.add_cable(1, 1, 3, 's')
f3.add_cable(2, 2, 3, 's')   # 正下方那條, 預期會鬆弛
f3.pin(0); f3.pin(1); f3.pin(2)
f3.point_load(3, fy=-P)
r3 = solve(f3)

print(f"最終判定鬆弛的cable: {r3.slack_cables}")
for mid in [0, 1, 2]:
    N = -r3.member_results[mid].end_forces_local[0]
    slack = r3.member_results[mid].slack
    print(f"  cable{mid}: N={N:.4f}  slack={slack}")
    if not slack:
        assert N >= -1e-6, f"沒被判定鬆弛的cable{mid}軸力應該>=0(拉力), 得到{N}"
    else:
        assert abs(N) < 1e-9, f"被判定鬆弛的cable{mid}軸力應該是0, 得到{N}"

assert 2 in r3.slack_cables, "正下方那條(對抵抗向下載重沒幫助的)應該被判定鬆弛"
print("正下方那條cable(對垂直向下載重沒有幫助)正確被判定鬆弛移除")

# 平衡檢查: 不管有幾條cable真正出力, 節點3的合力仍要跟外加載重平衡,
# 這點由solve()本身的K@u=F保證, 這裡額外用反力總和交叉確認整體垂直力平衡
total_Ry = sum(r3.reactions[f3.dofs_of(n)[1]] for n in [0, 1, 2])
print(f"三個支承垂直反力總和: {total_Ry:.4f}  (應該等於外加載重 {P})")
assert abs(total_Ry - P) < 1e-6, "整體垂直力不平衡"
print("PASS: 三纜線案例, 該鬆弛的cable被正確移除, 其餘cable正常受拉, 整體力平衡\n")

print("PASS: cable元素(纜線, 只受拉+自動鬆弛迭代)驗證完成")
