"""
交叉驗證: frame2d.cyclic (循環塑鉸) vs OpenSeesPy zeroLength 旋轉彈簧 + Steel01 (動力 D7)

這是「選用」測試: 沒有安裝 openseespy 時印出 SKIPPED 並正常結束。安裝: pip install openseespy

OpenSees 模型(跟 frame2d 的物理一模一樣): 樑本身 elasticBeamColumn, 每個桿端經過一個零長度
旋轉彈簧(zeroLength, dir 3)接到節點, 彈簧材料 Steel01 (Fy = Mp, E0 = K_e, b = R_post/K_e,
不開等向硬化 = 雙線性運動硬化)。平動用 equalDOF 綁在一起。位移控制反覆載重
(DisplacementControl), 底剪力取載重係數(單一參考載重 = 1)。

比對: 整條力-位移歷程逐點(每一段載重內, 用 frame2d 記錄的分段線性曲線內插到 OpenSees 的每個
步點)、各塑鉸端在轉折點的彎矩。兩個模型: 帶底部塑鉸的懸臂柱、6 個塑鉸的門型剛架(容量各不相同、
含不對稱與部分卸載的位移歷程)。

為什麼彈簧初始剛度用 K_e = 50·EI/L 而不是 hinge.py 預設的 1e8·EI/L:
  OpenSees 的 Newton 疊代在「近剛接彈簧 + 極小的降伏後剛度」(b = R/K_e ~ 1e-5)時, 多個塑鉸
  同時改變狀態會來回震盪、100 次疊代不收斂(K_e ≥ 1000·EI/L 就會發生; ModifiedNewton、
  NewtonLineSearch、KrylovNewton 都一樣)——這是 OpenSees 求解器的問題, 不是 frame2d 的。
  所以比對時兩邊用**同一個**有限的 K_e(把 hinge.RIGID_FACTOR 暫時設成 50), 兩個程式的模型
  仍然完全相同, 比對照樣有意義。
"""
import numpy as np

try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

import frame2d.hinge as hinge_mod
from frame2d import Frame2D
from frame2d.cyclic import CyclicHingeState, make_protocol, run_cyclic
from frame2d.dofmanager import initial_hinge_states

E, I, A = 200e6, 8e-5, 1e-2
KE_FACTOR = 50.0


def column():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, 4)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='s', Mp_i=100.0, Mp_j=None,
                 R_post_yield_i=2000.0, R_post_yield_j=None)
    f.fix(0)
    return f


def portal():
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4).add_node(3, 6, 0)
    f.add_section('s', E=E, I=I, A=A)
    caps = [(60.0, 68.0), (75.0, 85.0), (63.0, 71.0)]                     # 各鉸容量不同
    hard = [(1500.0, 1800.0), (1200.0, 1500.0), (1600.0, 1400.0)]
    for mid, (i, j) in enumerate([(0, 1), (1, 2), (3, 2)]):
        f.add_member(mid, node_i=i, node_j=j, section='s', Mp_i=caps[mid][0], Mp_j=caps[mid][1],
                     R_post_yield_i=hard[mid][0], R_post_yield_j=hard[mid][1])
    f.fix(0).fix(3)
    return f


def run_opensees(f, protocol, step, control_node):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    for nid, n in f.nodes.items():
        ops.node(nid, n.x, n.y)
    for s in f.supports:
        ops.fix(s.node, int(s.ux is not None), int(s.uy is not None), int(s.rot is not None))
    ops.geomTransf('Linear', 1)
    spring_tags, tag, mat = {}, 1000, 1
    for mid, m in f.members.items():
        sec = f.sections[m.section]
        ni, nj = f.nodes[m.node_i], f.nodes[m.node_j]
        L = float(np.hypot(nj.x - ni.x, nj.y - ni.y))
        Ke = KE_FACTOR * sec.E * sec.I / L
        dup = (100 + 2 * mid, 101 + 2 * mid)
        ops.node(dup[0], ni.x, ni.y)
        ops.node(dup[1], nj.x, nj.y)
        ops.equalDOF(m.node_i, dup[0], 1, 2)
        ops.equalDOF(m.node_j, dup[1], 1, 2)
        for end, real in enumerate((m.node_i, m.node_j)):
            Mp = (m.Mp_i, m.Mp_j)[end]
            R = (m.R_post_yield_i, m.R_post_yield_j)[end]
            if Mp is None:
                Mp, R = 1e9, 0.0                              # 沒有塑鉸容量的端點: 永遠彈性
            ops.uniaxialMaterial('Steel01', mat, Mp, Ke, max(R, 0.0) / Ke)
            ops.element('zeroLength', tag, real, dup[end], '-mat', mat, '-dir', 3)
            spring_tags[(mid, end)] = tag
            tag += 1
            mat += 1
        ops.element('elasticBeamColumn', mid, dup[0], dup[1], sec.A, sec.E, sec.I, 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(control_node, 1.0, 0.0, 0.0)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-10, 100)
    ops.algorithm('Newton')
    ops.analysis('Static')
    u, F, cur = [0.0], [0.0], 0.0
    turn_moments = []                                          # 每個轉折點的各鉸彎矩
    for tgt in protocol:
        n = max(1, int(round(abs(tgt - cur) / step)))
        ops.integrator('DisplacementControl', control_node, 1, (tgt - cur) / n)
        for _ in range(n):
            if ops.analyze(1) != 0:
                raise RuntimeError("OpenSees 這一步不收斂(請調小步長或改用較小的 K_e 因子)")
            u.append(ops.nodeDisp(control_node, 1))
            F.append(ops.getLoadFactor(1))
        turn_moments.append({k: ops.eleResponse(t, 'force')[2] for k, t in spring_tags.items()})   # zeroLength 'force' = [Fx,Fy,M]_i + [..]_j, 取節點 i 的彎矩
        cur = tgt
    return np.array(u), np.array(F), turn_moments


def compare(name, make_frame, protocol, d_nominal, step, control_node, base_nodes):
    f = make_frame()
    hs = CyclicHingeState.from_hinge_states(initial_hinge_states(f))
    dof = f.dofs_of(control_node)[0]
    res = run_cyclic(f, hs, [dof], [1.0], protocol, d_nominal, [f.dofs_of(b)[0] for b in base_nodes])
    uo, Fo, moments = run_opensees(make_frame(), protocol, step, control_node)

    # 每一段載重(相鄰目標之間): 以 frame2d 的分段線性曲線內插到 OpenSees 的每個步點
    li = [0]
    for t in protocol:
        i = li[-1]
        while i < len(res.u) - 1 and abs(res.u[i] - t) > 1e-12:
            i += 1
        li.append(i)
    lo, cur = [0], 0.0
    for t in protocol:
        lo.append(lo[-1] + max(1, int(round(abs(t - cur) / step))))
        cur = t
    worst = 0.0
    for a in range(len(protocol)):
        ua, Fa = res.u[li[a]:li[a + 1] + 1], res.F[li[a]:li[a + 1] + 1]
        ub, Fb = uo[lo[a]:lo[a + 1] + 1], Fo[lo[a]:lo[a + 1] + 1]
        if ua[0] > ua[-1]:
            ua, Fa = ua[::-1], Fa[::-1]
        order = np.argsort(ub)
        Fi = np.empty_like(Fb)
        Fi[order] = np.interp(ub[order], ua, Fa)
        worst = max(worst, float(np.max(np.abs(Fi - Fb))))
    Fmax = float(np.max(np.abs(Fo)))
    # 轉折點的各鉸彎矩(絕對值; 兩個程式對「彎矩正負號」的慣例不同)
    mworst = 0.0
    for a, tm in enumerate(moments):
        for (mid, e), Mo in tm.items():
            Mf = res.hinge_M[(mid, e)][li[a + 1]]
            if (mid, e) in res.hinge_M:
                mworst = max(mworst, abs(abs(Mf) - abs(Mo)))
    n_y = sum(1 for e in res.events if e['kind'] == 'yield')
    n_u = sum(1 for e in res.events if e['kind'] == 'unload')
    print(f"[{name}] frame2d {len(res.u)} 點 / OpenSees {len(uo)} 點; 降伏 {n_y} 次、卸載 {n_u} 次; "
          f"最大力 {Fmax:.4f}; 力-位移逐點最大差 {worst:.2e} (相對 {worst / Fmax:.2e}); 轉折點鉸彎矩最大差 {mworst:.2e}")
    assert n_y >= 2 and n_u >= 2, "測試案例要有降伏與卸載才有意義"
    assert worst < 1e-8 * Fmax, f"{name}: 力-位移與 OpenSeesPy 不一致"
    assert mworst < 1e-6 * max(abs(m) for tm in moments for m in tm.values()), f"{name}: 轉折點的鉸彎矩與 OpenSeesPy 不一致"


def main():
    if ops is None:
        print("SKIPPED: 沒有安裝 openseespy (pip install openseespy), 略過與 OpenSeesPy 的交叉驗證")
        return
    old = hinge_mod.RIGID_FACTOR
    hinge_mod.RIGID_FACTOR = KE_FACTOR
    try:
        L, ke = 4.0, 3 * E * I / 4.0**3
        uy = (100.0 / L) / ke
        compare("懸臂柱+底部塑鉸 等幅3圈", column, make_protocol([6 * uy], 3), uy / 3, uy / 20, 1, [0])
        compare("懸臂柱+底部塑鉸 遞增幅值", column, make_protocol([2 * uy, 4 * uy, 8 * uy], 2), uy / 3, uy / 20, 1, [0])
        compare("6塑鉸門型剛架 遞增幅值", portal, make_protocol([0.02, 0.05], 2), 0.004, 0.001, 1, [0, 3])
        compare("6塑鉸門型剛架 不對稱/部分卸載", portal, [0.03, -0.01, 0.05, 0.02, 0.04, -0.04, 0.0], 0.004, 0.001, 1, [0, 3])
    finally:
        hinge_mod.RIGID_FACTOR = old
    print("\n全部通過: 循環塑鉸的力-位移歷程與轉折點鉸彎矩跟 OpenSeesPy (zeroLength + Steel01) 一致到機器精度。")


if __name__ == "__main__":
    main()
