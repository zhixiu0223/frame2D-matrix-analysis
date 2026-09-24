"""
驗證案例: frame2d.ground_motion_io -- 讀取 PEER NGA .AT2 地震紀錄格式, 只依賴 numpy

層1 標準格式(有逗號、每行5筆): 解析出的 dt/加速度陣列逐點跟原始資料相符。
層2 格式容錯: 沒有逗號分隔、每行8筆、數值間距不固定, 都要能正確解析(PEER NGA 不同資料庫/
    工具匯出的檔案在這些地方常有差異, 見模組說明)。
層3 peer_nga_to_points(): 轉成 (t, ag) 點列表, t=i·dt 逐點驗證; max_points 抽稀時頭尾保留、
    點數不超過上限、原始資料點數不夠時不抽稀。
層4 明確拒絕: 找不到NPTS/DT那一行、NPTS<=0、DT<=0、數值筆數不夠。
"""
import numpy as np

from frame2d.ground_motion_io import parse_peer_nga, peer_nga_to_points


def make_at2(npts, dt, values, per_line=5, comma=True, sec_suffix=True):
    sep = ',' if comma else ''
    suffix = ' SEC' if sec_suffix else ''
    header = f"EXAMPLE RECORD\nSTATION INFO\nACCELERATION TIME HISTORY IN UNITS OF G\nNPTS={npts}{sep} DT={dt}{suffix}\n"
    lines = []
    for i in range(0, len(values), per_line):
        lines.append(" ".join(f"{v:.7f}" for v in values[i:i + per_line]))
    return header + "\n".join(lines) + "\n"


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------
print("=== 層1: 標準格式(逗號分隔, 每行5筆) ===")
dt = 0.01
n = 37
vals = (0.08 * np.sin(2 * np.pi * 1.3 * np.arange(n) * dt) * np.exp(-0.1 * np.arange(n) * dt)).round(7)
text = make_at2(n, dt, vals, per_line=5)
dt_got, acc_got = parse_peer_nga(text)
print(f"  dt = {dt_got}(給定 {dt}), 點數 = {len(acc_got)}(給定 {n})")
assert dt_got == dt and len(acc_got) == n
assert np.allclose(acc_got, vals), "逐點應該跟原始資料相符"

# 檔案末尾如果有多出來的數值(例如某些工具在資料後面留了一個換行或雜訊), 解析結果只能取前
# NPTS 個, 不能多拿(這是唯一能讓「多取1個」這類off-by-one錯誤露餡的情況, 一般乾淨的檔案
# 剛好只有NPTS個數值, 多拿1個時因為根本沒有第NPTS+1個token可拿, 不會露出任何差異)
text_extra = text.rstrip('\n') + f" {999.0}\n"
dt_e, acc_e = parse_peer_nga(text_extra)
assert len(acc_e) == n, f"檔案末尾多一個數值時, 還是應該只取前{n}個, 實際取了{len(acc_e)}個"
assert np.allclose(acc_e, vals), "多出來的那個數值不應該混進結果裡"

# ------------------------------------------------------------------
print("=== 層2: 格式容錯 ===")
for per_line, comma, sec in [(8, True, True), (5, False, True), (6, True, False), (1, False, False)]:
    text_v = make_at2(n, dt, vals, per_line=per_line, comma=comma, sec_suffix=sec)
    dt_v, acc_v = parse_peer_nga(text_v)
    ok = dt_v == dt and len(acc_v) == n and np.allclose(acc_v, vals)
    print(f"  每行{per_line}筆, 逗號={comma}, 'SEC'字樣={sec}: {'OK' if ok else 'FAIL'}")
    assert ok

# ------------------------------------------------------------------
print("=== 層3: peer_nga_to_points() ===")
text = make_at2(n, dt, vals)
pts = peer_nga_to_points(text)
assert len(pts) == n
for i in (0, 1, n // 2, n - 1):
    assert abs(pts[i][0] - i * dt) < 1e-12, f"t[{i}]"
    assert abs(pts[i][1] - vals[i]) < 1e-9, f"ag[{i}]"
print(f"  {n} 個點, 逐點 t=i·dt 且加速度值相符")

pts_ds = peer_nga_to_points(text, max_points=10)
assert len(pts_ds) <= 10
assert abs(pts_ds[0][0] - 0.0) < 1e-12 and abs(pts_ds[-1][0] - (n - 1) * dt) < 1e-12, "抽稀後頭尾要保留"
print(f"  抽稀到最多10點: 實際 {len(pts_ds)} 點, 頭尾都保留")

pts_nods = peer_nga_to_points(text, max_points=1000)   # 上限比實際點數大, 不該抽稀
assert len(pts_nods) == n, "上限比實際點數多時不應該抽稀"
print(f"  max_points比實際點數大時不抽稀: {len(pts_nods)} 點(原始 {n} 點)")

# ------------------------------------------------------------------
print("=== 層4: 明確拒絕 ===")


def expect_err(label, fn):
    try:
        fn()
    except ValueError as e:
        print(f"  {label}: ValueError OK ({str(e)[:40]}...)")
        return
    raise AssertionError(f"{label}: 應該報錯")


expect_err("找不到NPTS/DT那一行", lambda: parse_peer_nga("just some random text\nwith no header info\n"))
expect_err("NPTS<=0", lambda: parse_peer_nga("h\nh\nh\nNPTS=0, DT=0.01 SEC\n1 2 3\n"))
expect_err("DT<=0", lambda: parse_peer_nga("h\nh\nh\nNPTS=3, DT=0.0 SEC\n1 2 3\n"))
expect_err("數值筆數不夠", lambda: parse_peer_nga("h\nh\nh\nNPTS=10, DT=0.01 SEC\n1 2 3\n"))
expect_err("NPTS/DT那一行之後有非數字內容", lambda: parse_peer_nga("h\nh\nh\nNPTS=3, DT=0.01 SEC\n1 abc 3\n"))

print("\n全部通過: PEER NGA .AT2 標準格式、格式容錯、peer_nga_to_points抽稀都正確。")
