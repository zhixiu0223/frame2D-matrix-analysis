"""
交叉驗證: frame2d.spectrum.response_spectrum() vs OpenSeesPy `responseSpectrumAnalysis` (動力 D3)

這是「選用」測試: 沒有安裝 openseespy 時印出 SKIPPED 並正常結束。安裝: pip install openseespy

**範圍限制(誠實說明)**: 只比對**單一模態**的反應譜結果(位移、反力), 不比對多模態的 SRSS/CQC
組合。原因: OpenSeesPy 的 `ops.responseSpectrumAnalysis(tsTag, dir, '-mode', m1, m2, ...)` 在
探索階段測試了幾種呼叫方式(重複 `-mode` 參數對、單一 `-mode` 帶多個值、完全不給 `-mode`),
得到的位移都跟「只用最後(或某一個)單一模態」的結果一模一樣, 不像是把多個模態做了 SRSS/CQC
組合。這套件的說明文件不在這個環境裡, 我沒有把握找到正確的多模態語法, 與其比對一個我看不懂
語意、可能本來就沒在組合的數字, 不如誠實地只比對「這套件明確支援、行為可預期」的單一模態情況
—— 單一模態下兩者位移、反力都吻合到機器精度(1e-13), 這已經對 `response_spectrum()` 的核心公式
(K D_i = f_i 的等效靜力關係、S_d = S_a/ω²、反力慣例)提供了強力的外部驗證。多模態的 SRSS/CQC
組合改由 `tests/test_spectrum.py` 層2(兩質量懸臂, 獨立柔度矩陣解析解)驗證。

另外(跟 D2 的 `test_modal_vs_openseespy.py` 記錄過的同一件事): 一致質量下 OpenSees
`modalProperties`/`responseSpectrumAnalysis` 用的參與係數跟這裡的定義不完全一致(第1模態差
0.3%左右), 所以這裡的比對用**集中質量**, 兩邊逐項相同到 1e-13。
"""
import numpy as np

try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

from frame2d import Frame2D
from frame2d.modal import eigen
from frame2d.spectrum import response_spectrum

E, I, A, RHO, L = 200e6, 8e-5, 1e-2, 7.85, 3.0
N = 8


def build_frame2d():
    f = Frame2D()
    for i in range(N + 1):
        f.add_node(i, L * i / N, 0)
    f.add_section('s', E=E, I=I, A=A, rho=RHO)
    for i in range(N):
        f.add_member(i, i, i + 1, 's')
    f.fix(0)
    for i in range(1, N + 1):
        f.support(i, ux=0.0)
    return f


def build_opensees():
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    for i in range(N + 1):
        ops.node(i, L * i / N, 0.0)
    ops.fix(0, 1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i, 1, 0, 0)
    ops.geomTransf('Linear', 1)
    for i in range(N):
        ops.element('elasticBeamColumn', i, i, i + 1, A, E, I, 1, '-mass', RHO * A)


def main():
    if ops is None:
        print("SKIPPED: 沒有安裝 openseespy (pip install openseespy), 略過與 OpenSeesPy 的交叉驗證")
        return

    Sa0 = 5.0
    for mode in (1, 2, 3):
        f = build_frame2d()
        md = eigen(f, n_modes=mode, mass='lumped')
        res = response_spectrum(md, lambda T: Sa0, direction='y', n_modes=mode)
        # 只取第 mode 個模態自己的訊號值(不組合), 對應 OpenSeesPy 只算那一個模態的情況
        tip = f.dofs_of(N)[1]
        u_f2d = res.modal_displacements[mode - 1, tip]
        base_uy_f2d = res.modal_reactions[f.dofs_of(0)[1]][mode - 1]
        base_rot_f2d = res.modal_reactions[f.dofs_of(0)[2]][mode - 1]

        build_opensees()
        ops.eigen(mode)
        ops.modalProperties('-print')
        ops.timeSeries('Constant', 2, '-factor', Sa0)
        ops.responseSpectrumAnalysis(2, 2, '-mode', mode)
        u_os = ops.nodeDisp(N, 2)
        ops.reactions()
        base_uy_os = ops.nodeReaction(0, 2)
        base_rot_os = ops.nodeReaction(0, 3)

        d_disp = abs(u_f2d / u_os - 1.0) if u_os else abs(u_f2d)
        # 反力慣例相反號(見 test_spectrum.py 的說明: frame2d 用 R=Ku-F, 這裡兩邊都只比大小)
        d_uy = abs(abs(base_uy_f2d) / abs(base_uy_os) - 1.0) if base_uy_os else abs(base_uy_f2d)
        d_rot = abs(abs(base_rot_f2d) / abs(base_rot_os) - 1.0) if base_rot_os else abs(base_rot_f2d)
        print(f"[第 {mode} 個模態單獨] 位移相對差 {d_disp:.2e}, 底部反力(uy) 相對差 {d_uy:.2e}, "
              f"(rot) 相對差 {d_rot:.2e}")
        assert d_disp < 1e-9 and d_uy < 1e-9 and d_rot < 1e-9, f"第 {mode} 個模態與 OpenSeesPy 不一致"

    print("\n全部通過: 單一模態的反應譜位移與反力跟 OpenSeesPy responseSpectrumAnalysis 一致到機器精度。")


if __name__ == "__main__":
    main()
