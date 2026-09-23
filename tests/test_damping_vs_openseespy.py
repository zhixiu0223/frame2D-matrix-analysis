"""
交叉驗證: frame2d.damping (Rayleigh阻尼) + frame2d.excitation (地震輸入) vs OpenSeesPy
(動力 D5)

這是「選用」測試: 沒有安裝 openseespy 時印出 SKIPPED 並正常結束。安裝: pip install openseespy

**範圍限制(誠實說明)**: 只比對 **SDOF** 案例, 不比對 MDOF。原因: `ops.rayleigh(alphaM, betaK,
betaKinit, betaKcomm)` 的 βK 項是對**每個元素自己的(未凝縮)勁度矩陣**做比例阻尼, 不是對「消去
無質量自由度後的凝縮系統 K_eff」做比例阻尼。這兩者在 SDOF(只有一個動態自由度、沒有無質量自由
度需要凝縮)時完全相同, 但在 MDOF(例如集中質量下的懸臂梁, 轉角自由度沒有質量)時會有實質差異
——實測一個 6 元素懸臂梁案例, 兩邊差了 50%(不是數值誤差, 是方法論本身不同)。這正是
newmark.py 模組說明裡提到的「Rayleigh 阻尼應該建立在凝縮後的 K_eff 上, 不是原始 K 上」這件事的
反面案例: frame2d 刻意選擇在凝縮後的系統上做 Rayleigh 阻尼(見 damping.py 模組說明), OpenSees
的內建 `rayleigh()` 沒有這個選項。所以 MDOF 不放進機器精度比對, 改用 test_damping.py 層4 的
獨立模態投影檢核(用 modal.eigen() 算出的模態形狀驗證每個模態的阻尼比, 不依賴 OpenSees)。

在 SDOF 範圍內(從靜止開始、諧和地震輸入), frame2d 跟 OpenSeesPy 對到機器精度, 這對
rayleigh_coefficients()、ground_motion_force()、以及兩者透過 newmark_integrate() 串起來的
整條路徑都是很強的外部驗證。
"""
import numpy as np

try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

from frame2d import Frame2D
from frame2d.damping import rayleigh_damping_matrix
from frame2d.excitation import ground_motion_force, harmonic
from frame2d.modal import eigen
from frame2d.newmark import newmark_integrate

E, I, A, L = 200e6, 8e-5, 1e-2, 3.0
M_TIP = 2.0


def sdof_cantilever():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 's')
    f.fix(0)
    f.add_mass(1, my=M_TIP)
    return f


def main():
    if ops is None:
        print("SKIPPED: 沒有安裝 openseespy (pip install openseespy), 略過與 OpenSeesPy 的交叉驗證")
        return

    f = sdof_cantilever()
    w1 = eigen(f, n_modes=1).omega[0]
    zeta = 0.05
    A0, omega_f = 3.0, 0.65 * w1
    dt = (2 * np.pi / w1) / 173
    nsteps = 900

    f2 = sdof_cantilever()
    C, alpha, beta = rayleigh_damping_matrix(f2, w1, w1 * 1.3, zeta)
    force = ground_motion_force(f2, 'y', harmonic(A0, omega_f))
    res = newmark_integrate(f2, dt, nsteps, force=force, damping_matrix=C)
    tip = f2.dofs_of(1)[1]

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(0, 0, 0)
    ops.node(1, L, 0)
    ops.fix(0, 1, 1, 1)
    ops.mass(1, 0.0, M_TIP, 0.0)
    ops.geomTransf('Linear', 1)
    ops.element('elasticBeamColumn', 1, 0, 1, A, E, I, 1)
    ops.rayleigh(alpha, beta, 0.0, 0.0)
    ops.timeSeries('Trig', 1, 0.0, 1e9, 2 * np.pi / omega_f, '-factor', A0)
    ops.pattern('UniformExcitation', 1, 2, '-accel', 1)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-14, 50)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    us = [0.0]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            raise RuntimeError("OpenSees 這一步不收斂")
        us.append(ops.nodeDisp(1, 2))
    us = np.array(us)
    d = float(np.max(np.abs(res.u[:, tip] - us)))
    scale = float(np.max(np.abs(us)))
    print(f"[SDOF, 有阻尼(ζ={zeta})的諧和地震輸入] {nsteps} 步, 最大差 = {d:.3e}(相對振幅 {d / scale:.3e})")
    assert d < 1e-8 * scale, "與 OpenSeesPy 不一致"

    print("\n全部通過: SDOF有阻尼諧和地震輸入跟 OpenSeesPy(rayleigh()+UniformExcitation)一致到機器精度。")


if __name__ == "__main__":
    main()
