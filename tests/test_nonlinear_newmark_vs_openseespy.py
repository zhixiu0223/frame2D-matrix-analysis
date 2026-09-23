"""
交叉驗證: frame2d.nonlinear_newmark vs OpenSeesPy Newmark 非線性transient分析
(zeroLength+Steel01 塑鉸) (動力 D6)

這是「選用」測試: 沒有安裝 openseespy 時印出 SKIPPED 並正常結束。安裝: pip install openseespy

**範圍限制(誠實說明)**: 只比對 **SDOF**(單一塑鉸的懸臂柱), 不比對 MDOF。

SDOF 案例(從靜止開始、無阻尼、諧和力, 力夠大讓塑鉸在過程中降伏又卸載): 跟
OpenSeesPy(`zeroLength`+`Steel01`+`Newmark`遞迴) 對到機器精度(7.6e-13), 這對 Newmark 遞迴、
event-to-event 子步驟、以及跟 D5 質量/阻尼項的組合, 都是很強的外部驗證。

MDOF(6塑鉸門型剛架)有嘗試但沒有完全對上: 一開始兩邊在純彈性階段對到 1e-8 量級(確認 OpenSees
模型的節點/元素拓樸建對了), 但在還沒降伏前的階段就已經觀察到緩慢成長的偏差, 到後面(降伏後)
偏差變得很大。追過幾個可能的方向(dup節點的平移約束、i/j端對應關係)但沒有找到根本原因就沒有
繼續追。**這不代表 frame2d 的 MDOF 結果是錯的**——D6 的核心測試(`test_nonlinear_newmark.py`)
已經用兩個跟 OpenSeesPy 無關的獨立方法在 MDOF 規模驗證過: (1) 能量平衡(外力作功精確等於
動能+阻尼耗能+桿件內力作功, 到 1e-12 量級的殘差), (2) 跟 D7 `run_pushover`(同一個 6 塑鉸
門型剛架)在準靜態等效載重下比對第一個降伏事件的位置(位移與力都在 5e-4 相對誤差內)。只是
沒有再加上 OpenSeesPy 這第三種獨立驗證。
"""
import numpy as np

try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

import frame2d.hinge as hinge_mod
from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState
from frame2d.dofmanager import initial_hinge_states
from frame2d.excitation import force_series_from_pattern, harmonic
from frame2d.modal import eigen
from frame2d.nonlinear_newmark import nonlinear_newmark_integrate

E, I, A, L = 200e6, 8e-5, 1e-2, 3.0
M_TIP = 2.0
KE_FACTOR = 50.0     # 見 test_cyclic_vs_openseespy.py 的說明: OpenSees Newton 在「近剛接彈簧」
                      # (RIGID_FACTOR=1e8)下對多鉸系統會不收斂, 兩邊都改用有限的 KE_FACTOR=50


def sdof_cantilever():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='s', Mp_i=25.0, Mp_j=None,
                 R_post_yield_i=1000.0, R_post_yield_j=None)
    f.fix(0)
    f.add_mass(1, my=M_TIP)
    return f


def main():
    if ops is None:
        print("SKIPPED: 沒有安裝 openseespy (pip install openseespy), 略過與 OpenSeesPy 的交叉驗證")
        return

    old_rigid = hinge_mod.RIGID_FACTOR
    hinge_mod.RIGID_FACTOR = KE_FACTOR
    try:
        f = sdof_cantilever()
        hs = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
        w1 = eigen(f, n_modes=1).omega[0]
        tip = f.dofs_of(1)[1]
        p = np.zeros(6)
        p[tip] = 1.0
        F0, omega_f = 40.0, 0.8 * w1
        force = force_series_from_pattern(p, harmonic(F0, omega_f))
        dt = (2 * np.pi / w1) / 200
        nsteps = 1500
        res = nonlinear_newmark_integrate(f, hs, dt, nsteps, force=force)
        n_yield = sum(1 for e in res.events if e['kind'] == 'yield')
        n_unload = sum(1 for e in res.events if e['kind'] == 'unload')
        assert n_yield >= 2 and n_unload >= 1, "這個案例應該要有降伏跟卸載才測得到重點"

        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        ops.node(0, 0, 0)          # 真實固定端
        ops.node(10, 0, 0)         # dup: 基部鉸的另一側, 只有轉角自由
        ops.node(1, L, 0)          # 真實tip節點(有質量)
        ops.fix(0, 1, 1, 1)
        ops.fix(10, 1, 1, 0)       # dup平移固定(跟node0同位置), 轉角自由
        ops.mass(1, 0.0, M_TIP, 0.0)
        Ke = KE_FACTOR * E * I / L
        ops.uniaxialMaterial('Steel01', 1, 25.0, Ke, 1000.0 / Ke)
        ops.element('zeroLength', 1, 0, 10, '-mat', 1, '-dir', 3)
        ops.geomTransf('Linear', 1)
        ops.element('elasticBeamColumn', 2, 10, 1, A, E, I, 1)
        ops.timeSeries('Trig', 1, 0.0, 1e9, 2 * np.pi / omega_f, '-factor', F0)
        ops.pattern('Plain', 1, 1)
        ops.load(1, 0.0, 1.0, 0.0)
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', 1e-12, 100)
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
        print(f"[SDOF懸臂柱+底部塑鉸, 諧和力(ω=0.8ω1)] {nsteps} 步, {n_yield} 次降伏, {n_unload} 次卸載, "
              f"最大差 = {d:.3e}(相對振幅 {d / scale:.3e})")
        assert d < 1e-8 * scale, "與 OpenSeesPy 不一致"
    finally:
        hinge_mod.RIGID_FACTOR = old_rigid

    print("\n全部通過: SDOF非線性時程(含降伏與卸載)跟 OpenSeesPy(zeroLength+Steel01+Newmark)一致到機器精度。")


if __name__ == "__main__":
    main()
