"""
反應譜分析示範 (D3): 台灣規範風格的設計反應譜(簡化版, 不是精確查表)

跟 `cyclic_portal_demo.py`/上一輪模態 demo 同一個門型鋼架, 這次做反應譜分析: 算基底剪力、
屋頂位移估計、累積有效質量比, 並畫出設計反應譜與各模態落點。

**重要**: `taiwan_spectrum()` 只是示範用的簡化公式(固定形狀係數), **不是**「建築物耐震設計
規範及解說」的精確查表版本 —— 真正的 SDS/SD1(工址設計水平加速度反應譜)要依工址位置、
地盤分類查表得到。你在 taiwan-seismic-code-calc 那個 repo 已經有等值靜力法的完整實作,
之後可以把那份查表邏輯包成一個 `spectrum(T)` callable, 直接餵給這裡的 `response_spectrum()`
——這正是 ROADMAP 裡說的「規範譜做成 adapter, 不進核心」。

執行: PYTHONPATH=. python examples/response_spectrum_portal_demo.py
輸出: examples/output_plots/response_spectrum_portal.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from frame2d import Frame2D
from frame2d.modal import eigen
from frame2d.spectrum import response_spectrum

E, I, A, RHO = 200e9, 8e-5, 5e-3, 7850.0
H, W = 4.0, 6.0
G = 9.81                                     # m/s²


def build():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, H).add_node(2, W, H).add_node(3, W, 0)
    f.add_section('H', E=E, I=I, A=A, rho=RHO)
    f.add_member(0, 0, 1, 'H').add_member(1, 1, 2, 'H').add_member(2, 3, 2, 'H')
    f.fix(0).fix(3)
    f.add_mass(1, mx=5000.0, my=5000.0).add_mass(2, mx=5000.0, my=5000.0)   # 每層 5 t
    return f


def taiwan_spectrum(SDS=0.6, SD1=0.35, T0_ratio=1.0, TL=6.0):
    """簡化的設計水平加速度反應譜 S_aD(T)(單位 m/s², 已經乘上重力加速度 g), 形狀依「建築物耐震
    設計規範及解說」等值靜力法常見的四段式: 短週期上升段 / 平台段(=SDS) / 中長週期下降段(∝1/T,
    =SD1/T) / 長週期段(∝1/T², 只在 T>TL 才出現)。T0 = 0.2·SD1/SDS(平台段的下界週期), 這裡的
    T0_ratio 是額外的示範用縮放, 預設 1.0(不縮放)。"""
    T0 = T0_ratio * SD1 / SDS

    def spectrum(T):
        if T <= 0.2 * T0:
            return (0.4 + 0.6 * T / (0.2 * T0)) * SDS * G
        if T <= T0:
            return SDS * G
        if T <= TL:
            return SD1 * G / T
        return SD1 * TL * G / T**2
    return spectrum


def _pick_labels():
    cjk = [n for n in {f.name for f in font_manager.fontManager.ttflist}
           if any(k in n for k in ("CJK", "JhengHei", "PingFang", "Heiti", "Noto Sans TC", "Noto Serif TC"))]
    if cjk:
        plt.rcParams["font.sans-serif"] = sorted(cjk) + ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return {"t": "設計反應譜與各模態落點", "x": "週期 T (s)", "y": "S_a (g)",
                "leg_spec": "設計反應譜", "leg_mode": "模態"}
    return {"t": "Design spectrum and modal periods", "x": "Period T (s)", "y": "S_a (g)",
            "leg_spec": "Design spectrum", "leg_mode": "Mode"}


def main():
    f = build()
    md = eigen(f, n_modes=4, mass='consistent')
    spectrum = taiwan_spectrum()

    print("模態週期與 x 向累積有效質量比(相對總質量):")
    for i in range(md.n_modes):
        print(f"  模態{i + 1}: T={md.period[i]:.5f}s, S_a={spectrum(md.period[i]) / G:.3f}g, "
              f"累積={md.cum_ratio_total['x'][i]:.4f}")

    for combine in ('SRSS', 'CQC'):
        res = response_spectrum(md, spectrum, direction='x', damping=0.05, combine=combine)
        tip = f.dofs_of(1)[0]
        print(f"\n[{combine}] 累積質量比(用到的模態數={md.n_modes}) = {res.cum_ratio_total:.4f}"
              f"{' (未達規範常見門檻90%, 建議增加模態數)' if res.cum_ratio_total < 0.90 else ''}")
        print(f"[{combine}] 基底總剪力 = {res.base_shear / 1000:.2f} kN")
        print(f"[{combine}] 屋頂(N1)水平位移估計 = {res.displacements[tip] * 1000:.3f} mm")
        col_i = f.members[0]
        Mi = res.member_forces[0][2] / 1000
        print(f"[{combine}] 左柱柱底彎矩估計 = {Mi:.2f} kN·m")

    lab = _pick_labels()
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    Ts = np.linspace(0.001, 2.0, 400)
    ax.plot(Ts, [spectrum(T) / G for T in Ts], color="#7c3aed", lw=1.6, label=lab["leg_spec"])
    ax.scatter(md.period, [spectrum(T) / G for T in md.period], color="#ea580c", zorder=5, label=lab["leg_mode"])
    for i in range(md.n_modes):
        ax.annotate(f"{i + 1}", (md.period[i], spectrum(md.period[i]) / G),
                    textcoords="offset points", xytext=(5, 5), fontsize=9, color="#ea580c")
    ax.set_xlabel(lab["x"])
    ax.set_ylabel(lab["y"])
    ax.set_title(lab["t"])
    ax.set_xlim(0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "response_spectrum_portal.png")
    fig.savefig(out, dpi=130)
    print(f"\n圖已存到 {out}")


if __name__ == "__main__":
    main()
