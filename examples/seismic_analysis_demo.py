"""
一站式非線性地震反應分析示範 (D8): 用 `frame2d.seismic.seismic_analysis()` 一次呼叫做完
「模態→Rayleigh阻尼→地震輸入→非線性時程」, 畫出屋頂位移、能量平衡(輸入能量怎麼分配到
動能/阻尼耗能/桿件內力作功)、絕對加速度。

執行: PYTHONPATH=. python examples/seismic_analysis_demo.py
輸出: examples/output_plots/seismic_analysis_demo.png

**這裡的地震歷程是虛構的簡化脈衝(衰減正弦波), 不是真實地震紀錄** —— 只是用來展示
`seismic_analysis()` 怎麼用; 真實地震分析要接真實的地震歷程。這個示範刻意用跟
D6(`nonlinear_seismic_portal_demo.py`)同一個模型跟同一段地震歷程, 只是這次不用自己手動
串 eigen()→rayleigh_damping_matrix()→ground_motion_force()→nonlinear_newmark_integrate()
四步, 一行呼叫就好, 而且多了能量平衡的計算與檢核。
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from frame2d import Frame2D
from frame2d.seismic import seismic_analysis

E, I, A = 200e9, 8e-5, 5e-3
H, W = 4.0, 6.0
G = 9.80665
MP = 60e3
R_RATIO = 0.1
AMP_G = 0.5


def build():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A)
    for mid, (i, j, L) in enumerate([(0, 1, H), (1, 2, W), (3, 2, H)]):
        R = R_RATIO * E * I / L
        f.add_member(mid, node_i=i, node_j=j, section='H', Mp_i=MP, Mp_j=MP, R_post_yield_i=R, R_post_yield_j=R)
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


def synthetic_pulse(t):
    return AMP_G * G * np.sin(2 * np.pi * 1.5 * t) * np.exp(-0.25 * t)


def _pick_labels():
    cjk = [f.name for f in font_manager.fontManager.ttflist
           if any(k in f.name for k in ("CJK", "JhengHei", "PingFang", "Heiti", "Noto Sans TC", "Noto Serif TC"))]
    if cjk:
        plt.rcParams["font.sans-serif"] = sorted(set(cjk)) + ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return {"t0": "屋頂相對位移與絕對加速度", "t1": "能量平衡: 輸入能量怎麼分配",
                "x0": "時間 (s)", "y0a": "相對位移 (mm)", "y0b": "絕對加速度 (g)",
                "x1": "時間 (s)", "y1": "累積能量 (kN·m)",
                "l_ke": "動能", "l_wd": "阻尼耗能", "l_wi": "桿件內力作功(彈性+塑性)", "l_we": "外力(等效地震力)作功"}
    return {"t0": "Roof relative displacement & absolute acceleration", "t1": "Energy balance",
            "x0": "Time (s)", "y0a": "Rel. disp (mm)", "y0b": "Abs. accel (g)",
            "x1": "Time (s)", "y1": "Cumulative energy (kN·m)",
            "l_ke": "Kinetic", "l_wd": "Damping dissipated", "l_wi": "Member internal work", "l_we": "External (effective EQ force) work"}


def main():
    f = build()
    res = seismic_analysis(f, synthetic_pulse, direction='x', zeta=0.05)
    tip = f.dofs_of(1)[0]

    print(f"第1模態週期 T1 = {res.modal.period[0]:.4f} s")
    print(f"Rayleigh係數: alpha={res.alpha:.5f}, beta={res.beta:.6f}")
    print(f"尖峰相對位移 = {res.peak_displacement(1, 'x') * 1000:.2f} mm")
    print(f"尖峰層間位移角 = {res.peak_drift_ratio(1, 0, H, 'x') * 100:.3f} %")
    n_yield = sum(1 for e in res.events if e['kind'] == 'yield')
    print(f"降伏事件數 = {n_yield}")

    eb = res.energy_balance()
    resid = eb['Wext'] - (eb['KE'] + eb['Wdamp'] + eb['Wint'])
    print(f"能量平衡殘差(應接近0) = {np.max(np.abs(resid)):.3e}(外力作功尺度 {np.max(np.abs(eb['Wext'])):.1f} N·m)")
    print(f"最終: 動能={eb['KE'][-1]:.1f}, 阻尼耗能={eb['Wdamp'][-1]:.1f}, 桿件內力作功={eb['Wint'][-1]:.1f} N·m")

    a_abs = res.absolute_acceleration()

    lab = _pick_labels()
    fig, axs = plt.subplots(2, 1, figsize=(8, 7.5))
    ax0 = axs[0]
    ax0.plot(res.t, res.u[:, tip] * 1000, color="#0891b2", lw=1.1, label=lab["y0a"])
    ax0.axhline(0, color="#999", lw=0.6)
    ax0.set_xlabel(lab["x0"])
    ax0.set_ylabel(lab["y0a"], color="#0891b2")
    ax0.set_title(lab["t0"])
    ax0b = ax0.twinx()
    ax0b.plot(res.t, a_abs[:, tip] / G, color="#dc2626", lw=0.8, alpha=0.7, label=lab["y0b"])
    ax0b.set_ylabel(lab["y0b"], color="#dc2626")
    ax0.grid(alpha=0.3)

    ax1 = axs[1]
    ax1.stackplot(res.t, eb['KE'] / 1000, eb['Wdamp'] / 1000, eb['Wint'] / 1000,
                  labels=[lab["l_ke"], lab["l_wd"], lab["l_wi"]],
                  colors=["#0891b2", "#64748b", "#ea580c"], alpha=0.85)
    ax1.plot(res.t, eb['Wext'] / 1000, color="#111827", lw=1.3, ls="--", label=lab["l_we"])
    ax1.set_xlabel(lab["x1"])
    ax1.set_ylabel(lab["y1"])
    ax1.set_title(lab["t1"])
    ax1.legend(loc="upper left", fontsize=8)
    ax1.grid(alpha=0.3)

    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "seismic_analysis_demo.png")
    fig.savefig(out, dpi=130)
    print(f"\n圖已存到 {out}")


if __name__ == "__main__":
    main()
