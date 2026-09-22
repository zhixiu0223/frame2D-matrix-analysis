"""
線性時程分析示範 (D4): 同一個門型鋼架, 這次用 Newmark 積分算受諧和力/脈衝力時的位移時間歷程,
並跟 D3 的模態疊加互相對照(同一個模型、同一個結構, 只是分析方法不同)。

執行: PYTHONPATH=. python examples/newmark_portal_demo.py
輸出: examples/output_plots/newmark_portal_history.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from frame2d import Frame2D
from frame2d.excitation import force_series_from_pattern, harmonic, pulse
from frame2d.modal import eigen
from frame2d.newmark import newmark_integrate

E, I, A, RHO = 200e9, 8e-5, 5e-3, 7850.0
H, W = 4.0, 6.0


def build():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A, rho=RHO)
    f.add_member(0, 0, 1, 'H').add_member(1, 1, 2, 'H').add_member(2, 3, 2, 'H')
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)
    return f


def _pick_labels():
    cjk = [n for n in {f.name for f in font_manager.fontManager.ttflist}
           if any(k in n for k in ("CJK", "JhengHei", "PingFang", "Heiti", "Noto Sans TC", "Noto Serif TC"))]
    if cjk:
        plt.rcParams["font.sans-serif"] = sorted(cjk) + ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return {"t0": "諧和力(接近共振): 屋頂位移時間歷程", "t1": "脈衝力: 屋頂位移時間歷程",
                "x": "時間 (s)", "y": "屋頂(N1)水平位移 (mm)"}
    return {"t0": "Harmonic force (near resonance): roof displacement",
            "t1": "Pulse force: roof displacement", "x": "Time (s)", "y": "Roof (N1) horizontal disp (mm)"}


def main():
    f = build()
    md = eigen(f, n_modes=1, mass='consistent')
    T1, w1 = md.period[0], md.omega[0]
    print(f"第1模態: T1 = {T1:.4f} s, ω1 = {w1:.4f} rad/s")

    tip_dof_pattern = np.zeros(md.n_dof)
    tip_dof_pattern[f.dofs_of(1)[0]] = 1.0

    # ---- 案例A: 諧和力, 頻率比0.9(接近共振) ----
    F0 = 20000.0
    omega_f = 0.9 * w1
    force_h = force_series_from_pattern(tip_dof_pattern, harmonic(F0, omega_f))
    dt = T1 / 100
    n_cycles = 15
    nsteps = int(round(n_cycles * (2 * np.pi / omega_f) / dt))
    res_h = newmark_integrate(f, dt, nsteps, force=force_h, mass_kind='consistent')
    tip = f.dofs_of(1)[0]
    print(f"諧和力(ω=0.9ω1, F0={F0}N): 尖峰位移 = {np.max(np.abs(res_h.u[:, tip])) * 1000:.2f} mm")

    # ---- 案例B: 矩形脈衝力 ----
    F0b, t_start, duration = 40000.0, 0.0, T1 / 4
    force_p = force_series_from_pattern(tip_dof_pattern, pulse(F0b, t_start, duration))
    dt2 = T1 / 200
    nsteps2 = int(round(6 * T1 / dt2))
    res_p = newmark_integrate(build(), dt2, nsteps2, force=force_p, mass_kind='consistent')
    print(f"矩形脈衝(F0={F0b}N, 持續{duration * 1000:.1f}ms): 尖峰位移 = {np.max(np.abs(res_p.u[:, tip])) * 1000:.2f} mm")

    m0 = res_h.member_force_history(0)
    print(f"諧和力案例: 左柱柱底彎矩尖峰 = {np.max(np.abs(m0[:, 2])) / 1000:.2f} kN·m")

    lab = _pick_labels()
    fig, ax = plt.subplots(2, 1, figsize=(8, 7.5))
    ax[0].plot(res_h.t, res_h.u[:, tip] * 1000, color="#0891b2", lw=1.2)
    ax[0].axhline(0, color="#999", lw=0.6)
    ax[0].set_title(lab["t0"])
    ax[0].set_xlabel(lab["x"])
    ax[0].set_ylabel(lab["y"])
    ax[0].grid(alpha=0.3)
    ax[1].plot(res_p.t, res_p.u[:, tip] * 1000, color="#dc2626", lw=1.2)
    ax[1].axhline(0, color="#999", lw=0.6)
    ax[1].axvspan(t_start, t_start + duration, color="#dc2626", alpha=0.1)
    ax[1].set_title(lab["t1"])
    ax[1].set_xlabel(lab["x"])
    ax[1].set_ylabel(lab["y"])
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "newmark_portal_history.png")
    fig.savefig(out, dpi=130)
    print(f"\n圖已存到 {out}")


if __name__ == "__main__":
    main()
