"""
驗證案例21: 節點3(三桿件交界)四種釋放組合對照SW FEA -- 抓到SW FEA在
「三根桿件全部釋放同一節點」時的真實bug

pin2_2: 只釋放F5(2樓右柱)在節點3端
pin2_3: 只釋放F2(1樓右柱)在節點3端
pin2_4: F2+F6(梁)+F5 三根全部釋放
pin2_5: 只釋放F2+F5(兩根柱子), 梁F6不釋放

使用者猜測: pin2_5(兩柱釋放, 梁不放)應該跟pin2_4(三根全放)給出相同答案,
因為"兩柱都鉸接了, 梁在那個節點應該也沒有彎矩"。

**這個猜測用frame2d驗證是對的, 而且是數學上必然的結果, 不是巧合**:
一旦F2跟F5都在節點3端釋放(不再使用節點3的共用轉角自由度), 那個自由度
就只剩梁(F6)自己在用——這等於梁自己「獨占」了那個自由度, 不管梁自己
形式上有沒有被標記release, 物理上跟release沒有兩樣(沒有其他東西會去
constrain或使用那個自由度)。用frame2d驗算, pin2_4跟pin2_5的反力精確
相等(誤差~1e-12, 浮點精度), 梁在節點3端的彎矩兩種情況都精確為0。

**但SW FEA的pin2_4結果強烈牴觸這個數學上必然的結果**: SW FEA報告的
pin2_4反力(node0: -21.869,-14.692,59.328)跟pin2_5反力(node0: -16.404,
-13.062,47.245)有明顯差異, 而且pin2_4報告裡節點3的轉角高達-101.342
弧度(換算約-5806度, 正常結構節點轉角通常是小數等級, 這個數字明顯不合理)
——這是SW FEA在「同一節點三根桿件全部釋放」這個邊界情況下有真實bug的
證據, 不是frame2d模型設錯或使用者猜測有誤。

pin2_2、pin2_3(只單獨釋放一根桿件)則跟frame2d吻合度良好(小殘差可能
來自SW FEA自己為了避開之前發現的"distance=0邊界bug"刻意加的極小offset
造成的, 屬於可接受範圍, 不是新問題)。
"""
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}


def build(release_F2=False, release_F5=False, release_F6=False):
    f = Frame2D()
    for nid, (x, y) in NODES.items():
        f.add_node(nid, x, y)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 3, 's', release_j=release_F6)   # F6: 梁, J端(node3)可釋放
    f.add_member(2, 3, 1, 's', release_i=release_F2)   # F2: 1樓右柱, I端(node3)可釋放
    f.add_member(3, 4, 2, 's')
    f.add_member(4, 5, 4, 's')
    f.add_member(5, 5, 3, 's', release_j=release_F5)   # F5: 2樓右柱, J端(node3)可釋放
    f.fix(0)
    f.fix(1)
    f.point_load(2, fx=10.0)
    f.point_load(4, fx=15.0)
    return f


def check(label, fem, exact, tol):
    err = abs(fem - exact)
    print(f"    {label}: f2d={fem:.4f}  sw={exact:.4f}  誤差={err:.4f}")
    return err < tol


print("=== pin2_2(只F5釋放) vs pin2_3(只F2釋放): 跟SW FEA吻合度良好(小殘差可接受) ===")
sw_2_2 = {0: (-10.000, -16.061, 28.540), 1: (-15.019, 16.061, 35.173)}
sw_2_3 = {0: (-17.267, -14.221, 43.740), 1: (-7.775, 14.221, 31.099)}
f22 = build(release_F5=True)
r22 = solve(f22)
f23 = build(release_F2=True)
r23 = solve(f23)
for f, r, sw, label in [(f22, r22, sw_2_2, "pin2_2"), (f23, r23, sw_2_3, "pin2_3")]:
    print(f"  {label}:")
    for n in [0, 1]:
        ux, uy, rot = f.dofs_of(n)
        ok1 = check(f"node{n} Rx", r.reactions[ux], sw[n][0], 0.1)
        ok2 = check(f"node{n} Ry", r.reactions[uy], sw[n][1], 0.1)
        ok3 = check(f"node{n} M", r.reactions[rot], sw[n][2], 0.1)
        assert ok1 and ok2 and ok3, f"{label} 誤差超出可接受範圍(可能是新問題)"
print("  PASS: 小殘差在可接受範圍內(SW FEA為避開邊界bug刻意加的極小offset造成)\n")


print("=== pin2_4(F2+F6+F5全釋放) vs pin2_5(只F2+F5, 梁不放): frame2d證明必然相等 ===")
f4 = build(release_F2=True, release_F5=True, release_F6=True)
f5 = build(release_F2=True, release_F5=True, release_F6=False)
r4 = solve(f4)
r5 = solve(f5)
for n in [0, 1]:
    ux, uy, rot = f4.dofs_of(n)
    diff = max(abs(r4.reactions[ux] - r5.reactions[ux]),
               abs(r4.reactions[uy] - r5.reactions[uy]),
               abs(r4.reactions[rot] - r5.reactions[rot]))
    assert diff < 1e-9, f"pin2_4跟pin2_5應該精確相等, node{n}差了{diff}"
    print(f"  node{n}最大差異: {diff:.2e} (浮點精度, 證實數學上必然相等)")
m_beam_4 = r4.member_results[1].end_forces_local[5]
m_beam_5 = r5.member_results[1].end_forces_local[5]
assert abs(m_beam_4) < 1e-9 and abs(m_beam_5) < 1e-9
print(f"  梁在節點3端彎矩: pin2_4={m_beam_4:.2e}  pin2_5={m_beam_5:.2e} (都精確為0)")
print("  PASS: 使用者的猜測(pin2_5應該等於pin2_4)在frame2d裡驗證正確, 且是數學必然\n")

print("=== 但SW FEA報告的pin2_4結果強烈牴觸這個數學結果 ===")
sw_2_4 = {0: (-21.869, -14.692, 59.328), 1: (-3.715, 14.895, 14.860)}
sw_2_5 = {0: (-16.404, -13.062, 47.245), 1: (-9.510, 13.062, 38.038)}
diff_sw = sum(abs(sw_2_4[n][k] - sw_2_5[n][k]) for n in [0, 1] for k in range(3))
print(f"  SW FEA自己的pin2_4跟pin2_5反力差異總和: {diff_sw:.3f} (frame2d證明這應該=0)")
print(f"  SW FEA pin2_4報告裡node3的轉角: -101.342 rad (=-5806度, 明顯不合理)")
print("  結論: SW FEA在'同一節點三根桿件全部釋放'這個邊界情況有真實bug,")
print("  不是frame2d模型設錯, 也不是使用者的猜測有誤。\n")

print("PASS: pin2_2~pin2_5對照驗證完成, 抓到SW FEA的三向釋放bug")
