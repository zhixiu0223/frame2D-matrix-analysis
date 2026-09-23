"""
非線性地震反應示範 (D6): 同一個門型鋼架, 這次塑鉸容量設得比彈性需求低, 讓它在地震歷程中
真的降伏, 畫出力-位移遲滯迴圈與屋頂位移時間歷程, 跟 D5(純彈性)的結果對照。

執行: PYTHONPATH=. python examples/nonlinear_seismic_portal_demo.py
輸出: examples/output_plots/nonlinear_seismic_portal.png

**這裡的地震歷程是虛構的簡化脈衝(衰減正弦波), 不是真實地震紀錄** —— 只是用來展示
非線性時程分析怎麼用; 真實地震分析要接真實的地震歷程。
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState
from frame2d.damping import rayleigh_damping_matrix
from frame2d.dofmanager import initial_hinge_states
from frame2d.excitation import ground_motion_force
from frame2d.modal import eigen
from frame2d.newmark import newmark_integrate
from frame2d.nonlinear_newmark import nonlinear_newmark_integrate

E, I, A, RHO = 200e9, 8e-5, 5e-3, 7850.0
H, W = 4.0, 6.0
G = 9.80665
MP = 60e3                                    # N·m, 刻意設得比彈性需求(~91.5 kN·m)低, 才會降伏
R_RATIO = 0.1                                 # 降伏後硬化剛度 = R_RATIO * EI/L
AMP_G = 0.5                                   # 虛構地震脈衝的峰值加速度(g的倍數)


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
    """虛構的衰減正弦波地震脈衝, 只是示範用(跟 D5 demo 同一條)。"""
    return AMP_G * G * np.sin(2 * np.pi * 1.5 * t) * np.exp(-0.25 * t)


def _pick_labels():
    cjk = [f.name for f in font_manager.fontManager.ttflist
           if any(k in f.name for k in ("CJK", "JhengHei", "PingFang", "Heiti", "Noto Sans TC", "Noto Serif TC"))]
    if cjk:
        plt.rcParams["font.sans-serif"] = sorted(set(cjk)) + ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return {"t0": "屋頂位移時間歷程: 彈性 vs 非線性(有塑鉸)", "t1": "整體結構遲滯迴圈(非線性)",
                "t2": "左柱柱底塑鉸 M-θp", "x0": "時間 (s)", "y0": "屋頂位移 (mm)",
                "x1": "屋頂位移 (mm)", "y1": "柱底總彎矩 (kN·m)", "x2": "累積塑性轉角 (mrad)",
                "y2": "塑鉸彎矩 (kN·m)", "leg_e": "彈性(D5)", "leg_n": "非線性(D6)"}
    return {"t0": "Roof displacement: elastic vs nonlinear (hinges)", "t1": "Global hysteresis loop (nonlinear)",
            "t2": "Left column base hinge M-θp", "x0": "Time (s)", "y0": "Roof disp (mm)",
            "x1": "Roof disp (mm)", "y1": "Column base moment (kN·m)", "x2": "Cumulative plastic rotation (mrad)",
            "y2": "Hinge moment (kN·m)", "leg_e": "Elastic (D5)", "leg_n": "Nonlinear (D6)"}


def main():
    f = build()
    w1 = eigen(f, n_modes=1, mass='lumped').omega[0]
    w3 = eigen(f, n_modes=3, mass='lumped').omega[2]
    zeta = 0.05
    C, alpha, beta = rayleigh_damping_matrix(f, w1, w3, zeta, mass_kind='lumped')
    T1 = 2 * np.pi / w1
    print(f"第1模態週期 T1 = {T1:.4f} s, 塑鉸容量 Mp = {MP / 1000:.1f} kN·m")

    force_e = ground_motion_force(build(), 'x', synthetic_pulse, mass_kind='lumped')
    force_n = ground_motion_force(f, 'x', synthetic_pulse, mass_kind='lumped')
    dt = T1 / 80
    nsteps = int(round(10 / dt))

    res_elastic = newmark_integrate(build(), dt, nsteps, force=force_e, damping_matrix=C, mass_kind='lumped')
    hs = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
    res_nl = nonlinear_newmark_integrate(f, hs, dt, nsteps, force=force_n, damping_matrix=C, mass_kind='lumped')

    tip = f.dofs_of(1)[0]
    peak_e = np.max(np.abs(res_elastic.u[:, tip])) * 1000
    peak_n = np.max(np.abs(res_nl.u[:, tip])) * 1000
    n_yield = sum(1 for e in res_nl.events if e['kind'] == 'yield')
    print(f"彈性(D5)尖峰位移 = {peak_e:.2f} mm; 非線性(D6)尖峰位移 = {peak_n:.2f} mm, "
          f"{n_yield} 次降伏事件")
    key = (0, 0)
    residual = res_nl.hinge_theta_p[key][-1]
    print(f"左柱柱底塑鉸殘餘塑性轉角 = {residual * 1000:.4f} mrad")

    lab = _pick_labels()
    fig = plt.figure(figsize=(9, 8.5))
    ax0 = fig.add_subplot(2, 1, 1)
    ax0.plot(res_elastic.t, res_elastic.u[:, tip] * 1000, color="#0891b2", lw=1.0, label=lab["leg_e"])
    ax0.plot(res_nl.t, res_nl.u[:, tip] * 1000, color="#dc2626", lw=1.2, label=lab["leg_n"])
    ax0.axhline(0, color="#999", lw=0.6)
    ax0.set_title(lab["t0"])
    ax0.set_xlabel(lab["x0"])
    ax0.set_ylabel(lab["y0"])
    ax0.legend()
    ax0.grid(alpha=0.3)

    ax1 = fig.add_subplot(2, 2, 3)
    base_M = res_nl.member_force_history[0][:, 2] + res_nl.member_force_history[2][:, 2]
    ax1.plot(res_nl.u[:, tip] * 1000, base_M / 1000, color="#7c3aed", lw=1.0)
    ax1.axhline(0, color="#999", lw=0.6)
    ax1.axvline(0, color="#999", lw=0.6)
    ax1.set_title(lab["t1"])
    ax1.set_xlabel(lab["x1"])
    ax1.set_ylabel(lab["y1"])
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(2, 2, 4)
    ax2.plot(res_nl.hinge_theta_p[key] * 1000, res_nl.hinge_M[key] / 1000, color="#ea580c", lw=1.0)
    ax2.axhline(0, color="#999", lw=0.6)
    ax2.axvline(0, color="#999", lw=0.6)
    ax2.set_title(lab["t2"])
    ax2.set_xlabel(lab["x2"])
    ax2.set_ylabel(lab["y2"])
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "nonlinear_seismic_portal.png")
    fig.savefig(out, dpi=130)
    print(f"\n圖已存到 {out}")


if __name__ == "__main__":
    main()
