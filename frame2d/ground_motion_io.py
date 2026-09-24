"""
讀取真實地震紀錄檔案 (動力分析路線的延伸: 真實地震紀錄支援)。

目前支援 PEER NGA(太平洋地震工程研究中心強震資料庫)的 .AT2 文字格式, 這是結構工程界最
常用的公開強震紀錄格式之一(NGA-West1/West2/Sub 資料庫都用這個格式)。解析出來的
(dt, 加速度陣列)或 (t, ag) 點列表可以直接餵給 `seismic.custom_ground_motion()`, 接上
D5~D8 已經驗證過的非線性地震反應分析管線——**這個檔案只負責讀檔解析, 不含任何新的求解
邏輯**。

.AT2 格式(容許一定的格式差異, 見 `parse_peer_nga()` 的說明)
--------------------------------------------------------------
    第1行: 紀錄說明文字(地震名稱等, 自由格式)
    第2行: 資料來源說明
    第3行: 單位說明(通常是"ACCELERATION TIME HISTORY IN UNITS OF G"這類文字)
    第4行: 關鍵參數, 例如 "NPTS=  4000, DT=   .0100 SEC"(NPTS=點數, DT=取樣間隔秒)
    第5行以後: 加速度數值(單位: g), 每行固定筆數(常見5或8筆), 定寬或空白分隔

不同資料庫/工具匯出的檔案在「每行幾筆」「數值格式」上常有差異, 這裡的解析器**不假設固定
的每行筆數或欄寬**——只用正規表示式找出 NPTS/DT 那一行, 之後把剩下所有行的文字直接用空白
切開轉成浮點數, 取前 NPTS 個, 對格式差異有一定的容錯能力。
"""
import re

import numpy as np


def parse_peer_nga(text: str):
    """解析 PEER NGA .AT2 格式的文字內容, 回傳 (dt, accel_g)。

    text: 檔案的完整文字內容(自己讀檔案、貼上內容都可以, 這個函式只處理字串)。
    dt: 取樣間隔(s)。accel_g: numpy 陣列, 加速度數值(單位: g, 跟檔案裡的單位一致, 不做轉換
        ——PEER NGA 慣例上就是 g, 如果拿到的檔案單位不是g, 請先自己換算)。

    明確拒絕: 找不到 NPTS/DT 那一行、NPTS或DT不是正數、實際抓到的數值筆數比 NPTS 少。
    """
    m = re.search(r'NPTS\s*=\s*([\d.]+)\s*,?\s*DT\s*=\s*([\d.eE+-]+)', text, flags=re.IGNORECASE)
    if not m:
        raise ValueError("找不到 NPTS/DT 那一行(格式應該類似 'NPTS=4000, DT=0.01 SEC'), "
                         "這可能不是 PEER NGA .AT2 格式的檔案。")
    npts = int(float(m.group(1)))
    dt = float(m.group(2))
    if npts <= 0:
        raise ValueError(f"NPTS必須是正數, 解析到{npts}")
    if not (dt > 0 and np.isfinite(dt)):
        raise ValueError(f"DT必須是正數, 解析到{dt}")

    rest = text[m.end():]
    nl = rest.find('\n')                # 跳過NPTS/DT那一行剩下的文字(例如尾巴的"SEC"字樣)
    if nl != -1:
        rest = rest[nl + 1:]
    tokens = rest.split()
    try:
        values = [float(tok) for tok in tokens[:npts]]
    except ValueError as e:
        raise ValueError(f"NPTS/DT那一行之後的數值解析失敗(可能格式不對): {e}")
    if len(values) < npts:
        raise ValueError(f"NPTS={npts}, 但只解析到 {len(values)} 個數值, 檔案內容可能不完整。")

    return dt, np.array(values)


def peer_nga_to_points(text: str, max_points: int = None):
    """`parse_peer_nga()` 的結果轉成 `seismic.custom_ground_motion()` 要的 [(t, ag_g), ...]
    點列表。

    max_points: 選用, 真實強震紀錄常有上萬個取樣點, 直接全部塞進
        `custom_ground_motion()`(分段線性內插)雖然數學上沒問題, 但點數太多時(尤其是網頁
        傳輸)不必要地肥大——這裡如果給了 max_points 且原始點數超過它, 會**等間隔抽稀**
        (不是每隔N點簡單跳著取, 是均勻覆蓋整個時間範圍, 頭尾一定保留), 抽稀後波形仍然近似
        原始加速度歷程的包絡, 但不會逐點精確——需要逐點精確時不要設這個參數。
    """
    dt, accel_g = parse_peer_nga(text)
    n = len(accel_g)
    t = np.arange(n) * dt
    if max_points is not None and n > max_points:
        idx = np.unique(np.round(np.linspace(0, n - 1, max_points)).astype(int))
        t, accel_g = t[idx], accel_g[idx]
    return list(zip(t.tolist(), accel_g.tolist()))
