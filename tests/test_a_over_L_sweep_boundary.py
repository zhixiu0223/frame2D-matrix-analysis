"""
驗證案例27: a/L(鉸接分割長度/桿件全長比例)敏感度掃描, 精確定位雙桿件
模擬鉸接的數值安全邊界

背景: 使用者引用chatGPT的建議, 認為不該只測單一offset(0.0001m), 應該
做完整的a/L比例掃描, 才能精確定義frame2d的使用邊界。

**掃描結果(F5單獨鉸接案例, L=4m, 對照乾淨release_j)**:

  a/L        offset(m)    誤差(vs乾淨release_j)
  5e-1       2.000        8.374056   <- 這裡已經不是"逼近端點", 是真的
  2e-1       0.800        2.903053      移到桿件中段, 跟release_j本來就
  1e-1       0.400        1.317740      該有差異, 不算異常
  5e-2       0.200        0.625821
  1e-2       0.040        0.120026
  1e-3       0.004        0.011887   <- 誤差單調下降
  5e-4       0.002        0.005919
  1e-4       0.0004       0.001516
  5e-5       0.0002       0.000915   <- 誤差最低點(安全區底線)
  4e-5       0.00016      0.060167   <- 懸崖: 誤差瞬間放大66倍!
  2.5e-5     0.0001       0.060338   <- SW FEA實際使用的比例, 已經在懸崖裡
  2e-5       0.00008      0.543705
  1e-5       0.00004      1.438255
  1e-6       0.000004     88.010895  <- 完全失控

**結論: 這不是「offset越小越接近release_j」的單調關係, 是一個U型曲線**
- a/L從0.5降到5e-5: 誤差單調下降(offset越小越接近release_j, 符合直覺)
- a/L低於~4e-5: 誤差反轉暴增(病態矩陣開始主導, 不符合直覺, 但這正是
  chatGPT分析裡"L³在分母"造成的病態矩陣機制)
- 懸崖位置: 在a/L=5e-5(誤差0.0009)到a/L=4e-5(誤差0.060)之間(x86實測),
  誤差瞬間放大66倍以上。**這個位置是浮點雜訊決定的, 隨平台移動**: 在Termux
  (arm64, Python 3.14)上懸崖提早到a/L=1e-4~5e-5之間。測試因此只斷言結論
  (U型曲線、懸崖落在1e-5~1e-4), 不斷言恰好在哪一格

**這解釋了SW FEA殘差的精確成因**: SW FEA實際使用的offset(0.0001m,
對F5這根4m桿件而言a/L=2.5e-5)剛好落在懸崖裡面, 不是安全區的邊緣,
是已經進入惡化區間。

**使用邊界建議(寫入README.md)**:
  - a/L >= 1e-2(offset至少是桿長的1%): 完全安全, 誤差<0.15
  - a/L在1e-3~1e-4之間: 安全, 誤差<0.02
  - a/L < 5e-5: 不建議, 進入病態矩陣區間, 誤差開始不可預期地放大
  - 結論不變: 鉸接在桿件端點時, 直接用release_i/release_j(不需要
    offset, 不受這個問題影響); 只有鉸接真的在桿件中段時才需要節點
    分割, 且offset應該用有意義的長度(建議>=桿長的1%), 不要用極小值
    去逼近端點。
"""
from frame2d import Frame2D, solve

E, I, A = 200e6, 8e-5, 0.01
NODES = {0: (0, 0), 1: (6, 0), 2: (0, 4), 3: (6, 4), 4: (0, 8), 5: (6, 8)}
L = 4.0   # F5桿長


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


def build_split(offset):
    f = base_frame()
    f.add_node(7, 6, 4 + offset)
    f.add_member(0, 0, 2, 's')
    f.add_member(1, 2, 3, 's')
    f.add_member(2, 3, 1, 's')
    f.add_member(3, 4, 2, 's')
    f.add_member(4, 5, 4, 's')
    f.add_member(5, 5, 7, 's', release_j=True)
    f.add_member(10, 7, 3, 's')
    return f


f_clean = base_frame()
f_clean.add_member(0, 0, 2, 's')
f_clean.add_member(1, 2, 3, 's')
f_clean.add_member(2, 3, 1, 's')
f_clean.add_member(3, 4, 2, 's')
f_clean.add_member(4, 5, 4, 's')
f_clean.add_member(5, 5, 3, 's', release_j=True)
r_clean = solve(f_clean)
v_clean = tuple(r_clean.reactions[i] for i in f_clean.dofs_of(0))

RATIOS = [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 4e-5, 3.5e-5, 3e-5, 2.5e-5, 2e-5, 1e-5]


def analyze_sweep(sweep):
    """sweep: [(a/L, 誤差), ...], a/L由大到小。回傳 (安全區最佳誤差, 懸崖a/L, 懸崖誤差)。

    懸崖的**位置**是病態矩陣的數值雜訊決定的, 會隨平台(CPU架構、BLAS/LAPACK版本、
    Python/numpy版本)移動: x86 上在 5e-5→4e-5 之間, 在 Termux(arm64, Python 3.14)上
    已經提早到 1e-4→5e-5。所以這裡不再斷言「懸崖恰好在哪一格」, 只斷言結論本身:
      (1) a/L >= 2e-4 的安全區內, 誤差隨 offset 變小而持續下降
      (2) 往下走一定會遇到懸崖(誤差比前面最好的值放大10倍以上), 而且落在 1e-5~1e-4 之間
      (3) 懸崖前能達到的最佳誤差夠小(<0.005)
    """
    prev_err, best = None, None
    for ratio, err in sweep:
        if ratio >= 2e-4:
            if prev_err is not None:
                assert err <= prev_err * 1.5, f"安全區(a/L>=2e-4)內誤差應該持續下降: a/L={ratio:g}誤差{err}"
            prev_err = err
        if best is not None and err > 10 * best:
            return best, ratio, err
        best = err if best is None else min(best, err)
    return best, None, None


sweep = []
print(f"{'a/L':>12}{'offset(m)':>14}{'誤差':>16}")
for ratio in RATIOS:
    off = ratio * L
    f = build_split(off)
    r = solve(f)
    v = tuple(r.reactions[i] for i in f.dofs_of(0))
    err = sum(abs(a - b) for a, b in zip(v, v_clean))
    sweep.append((ratio, err))
    marker = "  <- SW FEA實際使用比例" if abs(ratio - 2.5e-5) < 1e-7 else ""
    print(f"{ratio:>12.2e}{off:>14.6f}{err:>16.6f}{marker}")

best, cliff_ratio, cliff_err = analyze_sweep(sweep)
assert cliff_ratio is not None, "應該要偵測到誤差暴增的懸崖"
assert 1e-5 <= cliff_ratio <= 1e-4, f"懸崖應該落在 a/L = 1e-5 ~ 1e-4 之間, 實際在 {cliff_ratio:g}"
assert best < 0.005, f"懸崖前的最佳誤差應該 < 0.005, 實際 {best}"
print(f"\n確認: 安全區最佳誤差 {best:.6f} -> 懸崖出現在 a/L={cliff_ratio:g} (誤差 {cliff_err:.6f}), "
      f"放大 {cliff_err / best:.1f} 倍。懸崖位置隨平台而異, 但「U型曲線、a/L 不可太小」的結論不變")

# 用另一個平台(Termux arm64 / Python 3.14)實測到的數值驗證這套判準對「懸崖提早出現」也成立
termux_sweep = [(1e-2, 0.120026), (5e-3, 0.059698), (2e-3, 0.023804), (1e-3, 0.011885),
                (5e-4, 0.005919), (2e-4, 0.002512), (1e-4, 0.000362), (5e-5, 0.043414)]
b2, c2, e2 = analyze_sweep(termux_sweep)
assert c2 == 5e-5 and b2 < 0.005, f"Termux實測數值應該判成: 懸崖在5e-5, 實際 {c2}, {b2}"
print(f"另一平台實測數值(Termux): 最佳誤差 {b2:.6f} -> 懸崖在 a/L={c2:g} (誤差 {e2:.6f}), 同一套判準成立")

print("\nPASS: a/L敏感度掃描完成, 安全邊界約在 a/L~1e-4 以上(懸崖的精確位置隨平台變動)")
print("(SW FEA實際使用的a/L=2.5e-5落在懸崖惡化區間內, 這解釋了之前發現的殘差成因)")
