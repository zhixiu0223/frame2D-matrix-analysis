"""
把 postprocess.py 已驗證過的內力分佈(N/V/M)跟變形形狀,整理成前端 SVG
畫圖需要的取樣點資料。不重新推導任何公式,只是呼叫既有函式 + 做自動縮放。

刻意跟 webapi_stdlib/diagrams.py 內容重複(不用 import 共用),因為
webapi/ 跟 webapi_stdlib/ 是兩個獨立的部署路徑(前者給 Colab/GCP 用
FastAPI,後者給 Termux 用純標準函式庫),不希望其中一個的存在與否
影響另一個能不能跑。
"""
import math

from frame2d.postprocess import member_internal_forces, member_deformed_shape


def build_diagrams_and_deformed(f, result, n_sample=21, target_fraction=0.15):
    """回傳 (diagrams, deformed, deform_scale)。

    diagrams: {member_id: {"x":[...], "N":[...], "V":[...], "M":[...]}}
        沿桿件局部座標取樣 n_sample 個點(重用 member_internal_forces,
        跟數值驗證用的是同一套函式,正負號慣例保證一致)。

    deformed: {member_id: {"X":[...], "Y":[...]}}
        桿件變形後的全域座標,縮放倍率是自動算出來的(讓最大變形量
        大約是結構整體外接矩形對角線的 target_fraction 倍,好讓變形
        曲線在畫面上看得出來,不是像真實比例那樣幾乎貼齊原線)。

    deform_scale: 實際套用的放大倍率(顯示在畫面上讓使用者知道
        「這是放大過的」,不是真實變形量)。
    """
    xs = [n.x for n in f.nodes.values()]
    ys = [n.y for n in f.nodes.values()]
    diag = math.hypot(max(xs) - min(xs), max(ys) - min(ys)) if xs else 1.0
    if diag < 1e-9:
        diag = 1.0

    raw_max = result.max_displacement()
    raw_max_val = raw_max.value if raw_max is not None else 0.0
    deform_scale = (diag * target_fraction / raw_max_val) if raw_max_val > 1e-12 else 1.0

    diagrams = {}
    deformed = {}
    for mid in f.members:
        x, N, V, M = member_internal_forces(f, result, mid, n=n_sample)
        diagrams[mid] = {
            "x": [float(v) for v in x],
            "N": [float(v) for v in N],
            "V": [float(v) for v in V],
            "M": [float(v) for v in M],
        }
        X, Y = member_deformed_shape(f, result, mid, scale=deform_scale, n=n_sample)
        deformed[mid] = {"X": [float(v) for v in X], "Y": [float(v) for v in Y]}

    return diagrams, deformed, deform_scale
