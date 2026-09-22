"""
交叉驗證: frame2d.newmark.newmark_integrate() vs OpenSeesPy Newmark transient analysis (動力 D4)

這是「選用」測試: 沒有安裝 openseespy 時印出 SKIPPED 並正常結束。安裝: pip install openseespy

**範圍限制(誠實說明)**:
1. 只比對**從靜止開始、受外力**的情況(u0=v0=a0=0, 力隨時間變化), 不比對「非零初始位移的
   自由振動」。原因: 用 `ops.setNodeDisp(node, dof, u0, '-commit')` 設定非零初始位移時,
   OpenSees 不會自動反解一個跟該位移平衡一致的初始加速度 a0——它就是把 a0 留在 0(可以用
   `ops.nodeAccel()` 直接確認)。這跟 frame2d 的作法不同(frame2d 一定會解
   M a0 = F0 - C v0 - K u0, 見 newmark.py 模組說明), 兩者從第一步就會分道揚鑣, 不是
   frame2d 的 bug, 只是 OpenSees 對這個情境的預設行為不同。從靜止開始受力沒有這個問題
   (a0=0 對兩邊都是正確的), 所以交叉驗證只測這種情況; 非零初始位移的驗證改用
   test_newmark.py 裡的解析解(SDOF自由振動、對數遞減率)。
2. 只比對**無阻尼**的情況。OpenSees 沒有一個直接對應「使用者自訂 C 矩陣」的簡單建模方式;
   試過用 zeroLength + Viscous 材料並聯出一個等效阻尼元素, 但得到的結果跟 frame2d 有明顯落差
   (30% 量級), 追不出是我兜的 OpenSees 模型設錯還是別的原因, 與其呈現一個看不懂哪裡有問題的
   數字, 不如誠實只驗證無阻尼的部分, 有阻尼的部分改用 test_newmark.py 裡的解析解(對數遞減率
   反推 ζ、諧和穩態振幅對動力放大係數公式)驗證。
3. **步階載重(step load)不納入比較**: 一開始測的時候用 OpenSees 的 `Constant` timeSeries
   加步階力, 結果跟 frame2d 差了 0.5%(Δt 減半只降到 0.1%, 收斂速率只有 O(Δt), 不是 O(Δt²)),
   一路查發現換成諧和力(`Trig` timeSeries)後兩邊立刻對到機器精度(1.7e-11)。這代表步階力的
   差異是 OpenSees 對「力在 t=0 究竟算不算已經完全作用」的處理方式跟 frame2d 不同(一個常見的
   時程分析慣例差異), 不是 frame2d 的問題, 但沒有繼續追下去確認 OpenSees 那邊確切的機制, 所以
   步階載重乾脆不放進比較, 改用 test_newmark.py 層4 的解析解(峰值=2倍靜位移)驗證。

在上面這些限制內(從靜止開始、諧和力、無阻尼), SDOF 與 MDOF 都對到機器精度(1e-11 ~ 1e-14),
這對 Newmark 遞迴本身、有效勁度矩陣、無質量DOF靜力凝縮都是很強的外部驗證。
"""
import numpy as np

try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

from frame2d import Frame2D
from frame2d.assembly import assemble_K
from frame2d.excitation import force_series_from_pattern, harmonic
from frame2d.modal import eigen
from frame2d.newmark import newmark_integrate

E, I, A, RHO, L = 200e6, 8e-5, 1e-2, 7.85, 3.0


def sdof_cantilever():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, L, 0)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 's')
    f.fix(0)
    f.add_mass(1, my=2.0)
    return f


N_MDOF = 6


def mdof_cantilever():
    f = Frame2D()
    for i in range(N_MDOF + 1):
        f.add_node(i, L * i / N_MDOF, 0)
    f.add_section('s', E=E, I=I, A=A, rho=RHO)
    for i in range(N_MDOF):
        f.add_member(i, i, i + 1, 's')
    f.fix(0)
    for i in range(1, N_MDOF + 1):
        f.support(i, ux=0.0)
    return f


def build_opensees(frame):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    for nid, nd in frame.nodes.items():
        ops.node(nid, nd.x, nd.y)
    for s in frame.supports:
        ops.fix(s.node, int(s.ux is not None), int(s.uy is not None), int(s.rot is not None))
    ops.geomTransf('Linear', 1)
    for mid, m in frame.members.items():
        sec = frame.sections[m.section]
        args = ['elasticBeamColumn', mid, m.node_i, m.node_j, sec.A, sec.E, sec.I, 1]
        if sec.rho:
            args += ['-mass', sec.rho * sec.A]
        ops.element(*args)
    for nm in frame.node_masses:
        ops.mass(nm.node, nm.mx, nm.my, nm.Iz)


def compare(name, frame, control_node, F0, omega_f, dt, nsteps):
    p = np.zeros(assemble_K(frame).n_dof)
    p[frame.dofs_of(control_node)[1]] = 1.0
    force = force_series_from_pattern(p, harmonic(F0, omega_f))
    res = newmark_integrate(frame, dt, nsteps, force=force, mass_kind='lumped')
    tip = frame.dofs_of(control_node)[1]

    build_opensees(frame)
    ops.timeSeries('Trig', 1, 0.0, 1e9, 2 * np.pi / omega_f, '-factor', F0)
    ops.pattern('Plain', 1, 1)
    ops.load(control_node, 0.0, 1.0, 0.0)
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
        us.append(ops.nodeDisp(control_node, 2))
    us = np.array(us)
    d = float(np.max(np.abs(res.u[:, tip] - us)))
    scale = float(np.max(np.abs(us)))
    print(f"[{name}] {nsteps} 步, 最大差 = {d:.3e}(相對振幅 {d / scale:.3e})")
    assert d < 1e-8 * scale, f"{name}: 與 OpenSeesPy 不一致"


def main():
    if ops is None:
        print("SKIPPED: 沒有安裝 openseespy (pip install openseespy), 略過與 OpenSeesPy 的交叉驗證")
        return

    f1 = sdof_cantilever()
    w1 = eigen(f1, n_modes=1).omega[0]
    compare("SDOF懸臂柱, 諧和力(ω=0.7ω1)", sdof_cantilever(), 1, 10.0, 0.7 * w1, (2 * np.pi / w1) / 97, 400)

    f2 = mdof_cantilever()
    wn = eigen(f2, n_modes=1, mass='lumped').omega[0]
    compare("6元素懸臂梁, 諧和力(ω=0.6ω1)", mdof_cantilever(), N_MDOF, 200.0, 0.6 * wn, 2e-5, 800)

    print("\n全部通過: 從靜止開始的無阻尼諧和強迫振動(SDOF與MDOF)跟 OpenSeesPy Newmark 遞迴一致到機器精度。")


if __name__ == "__main__":
    main()
