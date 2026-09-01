"""
把 postprocess.py 已驗證過的內力分佈(N/V/M)跟變形形狀,整理成前端 SVG
畫圖需要的取樣點資料。不重新推導任何公式,只是呼叫既有函式 + 做自動縮放。

刻意跟 webapi/diagrams.py 內容重複(見該檔案開頭說明)。
"""
import math

from frame2d.postprocess import member_internal_forces, member_deformed_shape


def build_diagrams_and_deformed(f, result, n_sample=21, target_fraction=0.15):
    """回傳 (diagrams, deformed, deform_scale)。定義同 webapi/diagrams.py。"""
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
