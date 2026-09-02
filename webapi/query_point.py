"""
查詢桿件內部「任意位置」(不是節點, 是桿件中間某一點)的內力(N/V/M)
跟位移量。不重新推導任何公式, 直接重用已驗證過的
postprocess.member_internal_forces() 跟 postprocess.member_deformed_shape(),
用夠密的取樣(n_sample)後在目標位置內插:
  - N, V 本身沿桿長是分段線性, 內插值是精確值
  - M 在均佈載重段內是二次曲線, 內插值有極小誤差(取樣越密誤差越小,
    n_sample=257 時對一般跨度已經遠低於工程有效位數)
  - 撓度(member_deformed_shape)本身就是用M(x)對EI做兩次數值積分算出來的
    精確解, 內插只是把它對齊到查詢的x, 不會額外損失精度

這裡不修改 postprocess.py, 純粹是查詢介面。
"""
import math

import numpy as np

from frame2d.postprocess import member_internal_forces, member_deformed_shape


def query_point(f, result, member_id: int, mode: str, value: float, n_sample: int = 257) -> dict:
    m = f.members[member_id]
    ni, nj = f.nodes[m.node_i], f.nodes[m.node_j]
    L = math.hypot(nj.x - ni.x, nj.y - ni.y)

    x_query = value * L if mode == "relative" else value
    if not (0.0 <= x_query <= L):
        raise ValueError(f"位置超出桿件範圍: 0 <= x <= {L:.4f} m, 收到 {x_query:.4f} m")

    x, N, V, M = member_internal_forces(f, result, member_id, n=n_sample)
    N_q = float(np.interp(x_query, x, N))
    V_q = float(np.interp(x_query, x, V))
    M_q = float(np.interp(x_query, x, M))

    X, Y = member_deformed_shape(f, result, member_id, scale=1.0, n=n_sample)
    Xd = float(np.interp(x_query, x, X))
    Yd = float(np.interp(x_query, x, Y))

    angle = math.atan2(nj.y - ni.y, nj.x - ni.x)
    ox = ni.x + x_query * math.cos(angle)
    oy = ni.y + x_query * math.sin(angle)
    dx, dy = Xd - ox, Yd - oy

    return {
        "member_id": member_id,
        "L": L,
        "x": x_query,
        "world_x": ox, "world_y": oy,
        "N": N_q, "V": V_q, "M": M_q,
        "dx": dx, "dy": dy,
        "disp": math.hypot(dx, dy),
    }
