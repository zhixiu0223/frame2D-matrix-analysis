"""
門型鋼架的遲滯迴圈示範 (D7: 循環塑鉸 + 準靜態反覆載重)

單層單跨門型鋼架, 6 個塑鉸(柱底、柱頂、梁兩端), 用「位移逐步增大的正反向循環」推屋頂節點,
畫出:
  左圖  底剪力 - 屋頂位移 遲滯迴圈(整體結構的遲滯行為)
  右圖  柱底塑鉸的 彎矩 - 塑性轉角 迴圈(單一構件的遲滯行為)
並印出: 降伏順序、每一圈的耗能(迴圈面積)與等效黏性阻尼比。

執行:   PYTHONPATH=. python examples/cyclic_portal_demo.py
輸出:   examples/output_plots/cyclic_portal_hysteresis.png

單位: SI (N, m); 圖上換成 kN、mm。模型參數(可自己改):
  跨度 6 m、柱高 4 m; E = 200 GPa, I = 8e7 mm⁴, A = 5000 mm²
  塑性彎矩 Mp = 100 kN·m; 降伏後硬化剛度 R = 0.1·EI/L (每根桿件依自己的長度)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState, make_protocol, run_cyclic
from frame2d.dofmanager import initial_hinge_states

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0
MP = 100e3                                   # N·m


def build():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    for mid, (i, j, L) in enumerate([(0, 1, H), (1, 2, W), (3, 2, H)]):
        R = 0.1 * E * I / L
        f.add_member(mid, node_i=i, node_j=j, section='H', Mp_i=MP, Mp_j=MP,
                     R_post_yield_i=R, R_post_yield_j=R)
    f.fix(0).fix(3)
    return f


def _pick_labels():
    """有中文字型(CJK)就用中文圖標, 沒有(例如 Termux 沒裝字型)就用英文, 避免圖上出現方框。"""
    cjk = [n for n in {f.name for f in font_manager.fontManager.ttflist}
           if any(k in n for k in ("CJK", "JhengHei", "PingFang", "Heiti", "Noto Sans TC", "Noto Serif TC"))]
    if cjk:
        plt.rcParams["font.sans-serif"] = sorted(cjk) + ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return {"x0": "屋頂位移 (mm)", "y0": "底剪力 (kN)", "t0": "整體結構遲滯迴圈",
                "x1": "累積塑性轉角 (mrad)", "y1": "塑鉸彎矩 (kN·m)", "t1": "左柱柱底塑鉸 M - θp"}
    return {"x0": "Roof displacement (mm)", "y0": "Base shear (kN)", "t0": "Global hysteresis loop",
            "x1": "Cumulative plastic rotation (mrad)", "y1": "Hinge moment (kN·m)",
            "t1": "Left column base hinge: M - θp"}


def main():
    lab = _pick_labels()
    f = build()
    hinge_states = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
    control_dof = f.dofs_of(1)[0]                                   # 屋頂左節點 ux
    base_dofs = [f.dofs_of(0)[0], f.dofs_of(3)[0]]                  # 兩個柱底的水平反力
    amplitudes = [0.015, 0.030, 0.045, 0.060]                       # m (位移比 0.4% ~ 1.5%)
    protocol = make_protocol(amplitudes, n_cycles=2)
    res = run_cyclic(f, hinge_states, [control_dof], [1.0], protocol, d_nominal=0.002,
                     base_reaction_dofs=base_dofs)

    print(f"共 {len(res.u)} 個記錄點, {sum(e['kind'] == 'yield' for e in res.events)} 次降伏事件, "
          f"{sum(e['kind'] == 'unload' for e in res.events)} 次卸載事件\n")

    print("降伏順序(第一次降伏):")
    seen = set()
    names = {0: "左柱", 1: "梁", 2: "右柱"}
    ends = {0: "i端", 1: "j端"}
    for ev in res.events:
        key = (ev['member'], ev['end'])
        if ev['kind'] == 'yield' and key not in seen:
            seen.add(key)
            print(f"  {len(seen)}. {names[ev['member']]} {ends[ev['end']]}: "
                  f"屋頂位移 {ev['u'] * 1000:6.2f} mm, 底剪力 {ev['F'] / 1000:7.2f} kN")

    # 每一圈(從 +a 到下一個 +a)的耗能與等效黏性阻尼比
    print("\n每個幅值最後一圈的耗能與等效黏性阻尼比 ξ_eq = E_D / (4π·E_S), E_S = ½·F_max·u_max:")
    print(f"  {'幅值(mm)':>9} {'F_max(kN)':>10} {'耗能 E_D(kJ)':>13} {'ξ_eq':>7}")
    peaks = [i for i in range(1, len(res.u) - 1)
             if res.u[i] >= res.u[i - 1] and res.u[i] >= res.u[i + 1] and res.u[i] > 1e-9]
    for a in amplitudes:
        idx = [i for i in peaks if abs(res.u[i] - a) < 1e-9]
        i0, i1 = idx[-2], idx[-1]                                    # 該幅值最後一圈
        Ed = res.loop_work(i0, i1)
        Fmax = max(res.F[i0:i1 + 1])
        xi = Ed / (4 * np.pi * 0.5 * Fmax * a)
        print(f"  {a * 1000:9.1f} {Fmax / 1000:10.2f} {Ed / 1000:13.3f} {xi:7.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ax[0].plot(res.u * 1000, res.F / 1000, lw=1.4, color="#7c3aed")
    ax[0].axhline(0, color="#999", lw=0.6)
    ax[0].axvline(0, color="#999", lw=0.6)
    ax[0].set_xlabel(lab["x0"])
    ax[0].set_ylabel(lab["y0"])
    ax[0].set_title(lab["t0"])
    ax[0].grid(alpha=0.3)
    key = (0, 0)                                                      # 左柱柱底塑鉸
    ax[1].plot(res.hinge_theta_p[key] * 1000, res.hinge_M[key] / 1000, lw=1.4, color="#ea580c")
    ax[1].axhline(0, color="#999", lw=0.6)
    ax[1].axvline(0, color="#999", lw=0.6)
    ax[1].set_xlabel(lab["x1"])
    ax[1].set_ylabel(lab["y1"])
    ax[1].set_title(lab["t1"])
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "cyclic_portal_hysteresis.png")
    fig.savefig(out, dpi=130)
    print(f"\n圖已存到 {out}")


if __name__ == "__main__":
    main()
