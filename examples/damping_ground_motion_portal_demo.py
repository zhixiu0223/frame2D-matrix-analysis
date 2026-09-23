"""
Rayleigh 阻尼 + 地震輸入示範 (D5): 同一個門型鋼架, 這次加上阻尼、用一段簡化的地震脈衝
(衰減正弦波, 不是真實地震歷程, 只是示範用)當地面加速度輸入, 算屋頂的相對位移與絕對加速度。

執行: PYTHONPATH=. python examples/damping_ground_motion_portal_demo.py
輸出: examples/output_plots/damping_ground_motion_portal.png

**這裡的地震歷程是虛構的簡化脈衝(衰減正弦波), 不是真實地震紀錄** —— 只是用來展示
Rayleigh 阻尼 + 地面加速度輸入怎麼串起來用; 真實地震分析要接真實的地震歷程(PEER NGA、
921 集集地震紀錄等)進 `ground_motion_force()`。
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from frame2d import Frame2D
from frame2d.damping import rayleigh_damping_matrix, rayleigh_damping_ratio
from frame2d.excitation import absolute_acceleration, ground_motion_force
from frame2d.modal import eigen
from frame2d.newmark import newmark_integrate

E, I, A, RHO = 200e9, 8e-5, 5e-3, 7850.0
H, W = 4.0, 6.0
G = 9.80665


def build():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A, rho=RHO)
    f.add_member(0, 0, 1, 'H').add_member(1, 1, 2, 'H').add_member(2, 3, 2, 'H')
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


def synthetic_pulse(pga_g=0.3, freq_hz=2.0, decay=0.3):
    """虛構的衰減正弦波地震脈衝, ag(t) = pga·g·sin(2πf t)·e^(-decay·t)。只是示範用。"""
    amp = pga_g * G
    return lambda t: amp * np.sin(2 * np.pi * freq_hz * t) * np.exp(-decay * t)


def _pick_labels():
    cjk = [f.name for f in font_manager.fontManager.ttflist
           if any(k in f.name for k in ("CJK", "JhengHei", "PingFang", "Heiti", "Noto Sans TC", "Noto Serif TC"))]
    if cjk:
        plt.rcParams["font.sans-serif"] = sorted(set(cjk)) + ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return {"t0": "地面加速度輸入(虛構脈衝)", "t1": "屋頂相對位移", "t2": "屋頂絕對加速度",
                "x": "時間 (s)", "y0": "地面加速度 (g)", "y1": "相對位移 (mm)", "y2": "絕對加速度 (g)"}
    return {"t0": "Ground acceleration (synthetic pulse)", "t1": "Roof relative displacement",
            "t2": "Roof absolute acceleration", "x": "Time (s)", "y0": "Ground accel (g)",
            "y1": "Rel. disp (mm)", "y2": "Abs. accel (g)"}


def main():
    f = build()
    md = eigen(f, n_modes=4, mass='consistent')
    print("模態週期:", [f"{t:.4f}" for t in md.period])

    zeta = 0.05
    C, alpha, beta = rayleigh_damping_matrix(f, md.omega[0], md.omega[2], zeta, mass_kind='consistent')
    print(f"Rayleigh係數: alpha={alpha:.5f}, beta={beta:.6f}")
    for i in range(md.n_modes):
        print(f"  模態{i + 1} T={md.period[i]:.4f}s 實際阻尼比={rayleigh_damping_ratio(alpha, beta, md.omega[i]):.4f}")

    ag = synthetic_pulse()
    force = ground_motion_force(f, 'x', ag, mass_kind='consistent')
    dt = md.period[0] / 100
    nsteps = int(round(10 / dt))
    res = newmark_integrate(f, dt, nsteps, force=force, damping_matrix=C, mass_kind='consistent')
    tip = f.dofs_of(1)[0]
    a_abs = absolute_acceleration(res, 'x', ag)

    peak_disp = np.max(np.abs(res.u[:, tip])) * 1000
    peak_abs_accel = np.max(np.abs(a_abs[:, tip])) / G
    print(f"屋頂尖峰相對位移 = {peak_disp:.2f} mm, 尖峰絕對加速度 = {peak_abs_accel:.3f} g")

    lab = _pick_labels()
    fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    ag_series = np.array([ag(tt) for tt in res.t]) / G
    axs[0].plot(res.t, ag_series, color="#64748b", lw=1.0)
    axs[0].set_title(lab["t0"])
    axs[0].set_ylabel(lab["y0"])
    axs[0].grid(alpha=0.3)
    axs[1].plot(res.t, res.u[:, tip] * 1000, color="#0891b2", lw=1.1)
    axs[1].axhline(0, color="#999", lw=0.6)
    axs[1].set_title(lab["t1"])
    axs[1].set_ylabel(lab["y1"])
    axs[1].grid(alpha=0.3)
    axs[2].plot(res.t, a_abs[:, tip] / G, color="#dc2626", lw=1.1)
    axs[2].axhline(0, color="#999", lw=0.6)
    axs[2].set_title(lab["t2"])
    axs[2].set_xlabel(lab["x"])
    axs[2].set_ylabel(lab["y2"])
    axs[2].grid(alpha=0.3)
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "damping_ground_motion_portal.png")
    fig.savefig(out, dpi=130)
    print(f"\n圖已存到 {out}")


if __name__ == "__main__":
    main()
