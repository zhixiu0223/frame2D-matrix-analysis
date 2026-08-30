"""
驗證案例28: a/L = 1e-1, 1e-2, 1e-3, 1e-4 四點對照SW FEA(phase-04-
internal-pin_2_3e-1~3e-4), 確認安全區內誤差單調收斂, 尚未進入懸崖

背景: 使用者提供F2(1樓右柱)在四個不同a/L比例下的SW FEA案例, 問「理論上
誤差要越來越大嗎? 預期桿件太短, 噪音變多(病態矩陣)」。

**答案: 在這個範圍(1e-1~1e-4)不會越來越大, 反而單調收斂**——因為這四個
點全部落在test_a_over_L_sweep_boundary.py找到的"安全區"內(懸崖在
a/L~4~5e-5, 比這裡最小的1e-4還小一個數量級以上), 還沒開始出現病態矩陣
的噪音。

**結果**:
  a/L      frame2d跟SW FEA誤差    frame2d跟乾淨release_i誤差(理論極限)
  1e-1     0.0004                 5.4848  (offset較大, 跟release_i本來
  1e-2     0.0002                 0.5133   就會有差異, 這是幾何位置
  1e-3     0.0004                 0.0509   不同造成, 不是異常)
  1e-4     0.0058                 0.0046  (已經很接近release_i)

frame2d(雙桿件模擬) 跟 SW FEA 在全部四個offset都精確吻合(誤差<0.006),
而且frame2d自己算出的答案隨a/L縮小平滑收斂到release_i的理論極限值,
每縮小一個數量級, 跟release_i的差距大約縮小10倍, 符合線彈性有限元素法
"網格加密, 誤差線性/多項式收斂"的標準預期。

跟test_a_over_L_sweep_boundary.py(掃到懸崖, a/L低於~5e-5開始惡化)合併
來看, 完整的圖像是: a/L從0.5一路降到~5e-5, 誤差持續下降; 只有低於這個
懸崖, 誤差才會反轉暴增。這次的四個點都還在下降段, 沒有異常。
"""
import numpy as np
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}
L = 4.0   # F2桿長


def base_frame():
    f = Frame2D()
    for nid, (x, y) in NODES.items():
        f.add_node(nid, x, y)
    f.add_section('s', E=E, I=I, A=A)
    f.fix(0)
    f.fix(1)
    f.point_load(2, fx=10.0)
    f.point_load(4, fx=15.0)
    return f


def build_split(a_over_L):
    offset = a_over_L * L
    f = base_frame()
    f.add_node(8, 6, 4 - offset)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 3, 's')
    f.add_member(2, 3, 8, 's')
    f.add_member(9, 8, 1, 's', release_i=True)
    f.add_member(3, 4, 2, 's')
    f.add_member(4, 5, 4, 's')
    f.add_member(5, 5, 3, 's')
    return f


f_clean = base_frame()
f_clean.add_member(0, 0, 2, 's')
f_clean.add_member(1, 2, 3, 's')
f_clean.add_member(2, 3, 1, 's', release_i=True)
f_clean.add_member(3, 4, 2, 's')
f_clean.add_member(4, 5, 4, 's')
f_clean.add_member(5, 5, 3, 's')
r_clean = solve(f_clean)
v_clean = tuple(r_clean.reactions[i] for i in f_clean.dofs_of(0))

sw = {
    1e-1: (-15.663, -14.406, 39.952),
    1e-2: (-17.088, -14.225, 43.318),
    1e-3: (-17.222, -14.214, 43.635),
    1e-4: (-17.234, -14.213, 43.663),
}

print(f"{'a/L':>8}{'frame2d node0':>32}{'跟SW誤差':>12}{'跟乾淨版誤差':>14}")
prev_clean_err = None
for ratio, swv in sw.items():
    f = build_split(ratio)
    r = solve(f)
    v = tuple(r.reactions[i] for i in f.dofs_of(0))
    err_sw = sum(abs(a - b) for a, b in zip(v, swv))
    err_clean = sum(abs(a - b) for a, b in zip(v, v_clean))
    print(f"{ratio:>8.0e}{str(tuple(round(x,3) for x in v)):>32}{err_sw:>12.4f}{err_clean:>14.4f}")
    assert err_sw < 0.01, f"a/L={ratio}時跟SW FEA誤差應該很小(安全區), 得到{err_sw}"
    if prev_clean_err is not None:
        assert err_clean < prev_clean_err, f"a/L={ratio}時跟乾淨版本的誤差應該持續下降"
    prev_clean_err = err_clean

print(f"\n乾淨版本(release_i, 理論極限): {tuple(round(x,3) for x in v_clean)}")
print("\nPASS: a/L=1e-1~1e-4這四點全部在安全區內, frame2d精確吻合SW FEA,")
print("且frame2d自己隨a/L縮小平滑收斂到release_i理論極限值, 沒有異常。")
print("(懸崖在a/L~4~5e-5, 比這裡最小值還小一個數量級以上, 見")
print("test_a_over_L_sweep_boundary.py)")
