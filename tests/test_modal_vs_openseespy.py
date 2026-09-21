"""
交叉驗證: frame2d.modal.eigen() vs OpenSeesPy (動力分析 D2)

這是「選用」測試: 沒有安裝 openseespy 時印出 SKIPPED 並正常結束(不算失敗)。
安裝: pip install openseespy

比對內容(集中質量 lumped 與一致質量 consistent 各一遍, 兩個模型):
  1. 特徵值 ω²: 相對誤差
  2. 模態形狀: MAC (模態保證準則 |vᵀu|²/(vᵀv·uᵀu)), 不比數值(正規化與正負號慣例不同)
  3. 質量與參與係數(OpenSees modalProperties):
     - 總質量 totalMass、可動質量 totalFreeMass、|Γ|、有效質量
     - 集中質量: 逐項一致
     - 一致質量: 特徵值與模態形狀仍一致, 總質量一致, 但**參與係數不會逐項一致**:
       OpenSees modalProperties 對一致質量的可動質量/參與係數採用另一套算法(可動質量等於
       集中質量的值, 我沒有完整重現它的定義); frame2d 用完整地面運動向量 r, Γ = φᵀ M r,
       含一致質量與支承DOF的耦合項, 這是均勻地面運動下的精確有效地震力 -M r a_g, 並且
       收斂到連續體解(test_modal.py 層2: 懸臂梁第1模態有效質量比 0.61308, 相對差 6e-10)。
       這裡只確認兩者差異很小(第1模態 |Γ| 在 2% 內), 不宣稱相等。

重要的實務細節(踩過的坑): OpenSeesPy 預設的 eigen 求解器(ARPACK)在模型有「無質量DOF」時
(集中質量下的轉角自由度)會失敗:
    ArpackSolver::Error ... Could not build an Arnoldi factorization
必須改用  ops.eigen('-fullGenLapack', n)  (很慢, 只適合小模型, 但這裡夠用)。
"""
import numpy as np

# 注意: 不能在模組層級 sys.exit()。pytest 收集測試時會 import 每個 test_*.py,
# 模組層級的 SystemExit 會讓整個 pytest 以 INTERNALERROR 中止。所以沒裝就只設 ops=None,
# 所有實際工作放在 main() 裡, 只有直接執行(或被 test_zz 的 subprocess 執行)才會跑。
try:
    import openseespy.opensees as ops
except ImportError:
    ops = None

from frame2d import Frame2D
from frame2d.modal import eigen

E, I, A, RHO, L = 200e6, 8e-5, 1e-2, 7.85, 3.0


def cantilever(n=8):
    f = Frame2D()
    for i in range(n + 1):
        f.add_node(i, L * i / n, 0)
    f.add_section('s', E=E, I=I, A=A, rho=RHO)
    for i in range(n):
        f.add_member(i, i, i + 1, 's')
    f.fix(0)
    for i in range(1, n + 1):
        f.support(i, ux=0.0)
    return f


def mixed_frame():
    """斜桿(frame與truss)、不同斷面、含轉動慣量的節點質量, 結構不對稱(模態不簡併)。"""
    f = Frame2D()
    f.add_node(0, 0, 0).add_node(1, 0, 4).add_node(2, 6, 4.5).add_node(3, 6.5, 0).add_node(4, 3, 7.5)
    f.add_section('c', E=E, I=I, A=A, rho=RHO)
    f.add_section('b', E=E, I=2 * I, A=2 * A, rho=RHO)
    f.add_section('t', E=E, I=I, A=1e-3, rho=RHO)
    f.add_member(0, 0, 1, 'c').add_member(1, 1, 2, 'b').add_member(2, 3, 2, 'c').add_member(3, 1, 4, 'b')
    f.add_truss(4, 4, 2, 't')
    f.fix(0).fix(3)
    f.add_mass(1, mx=3.0, my=1.0, Iz=0.5).add_mass(2, mx=2.0, my=2.5).add_mass(4, mx=1.5, my=1.0, Iz=0.1)
    return f


def build_opensees(f, kind):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    for nid, n in f.nodes.items():
        ops.node(nid, n.x, n.y)
    for s in f.supports:
        ops.fix(s.node, int(s.ux is not None), int(s.uy is not None), int(s.rot is not None))
    ops.geomTransf('Linear', 1)
    ops.uniaxialMaterial('Elastic', 1, E)
    for mid, m in f.members.items():
        sec = f.sections[m.section]
        if m.member_type == 'truss':
            args = ['Truss', mid, m.node_i, m.node_j, sec.A, 1, '-rho', sec.rho * sec.A]
            if kind == 'consistent':
                args += ['-cMass', 1]
        else:
            args = ['elasticBeamColumn', mid, m.node_i, m.node_j, sec.A, sec.E, sec.I, 1, '-mass', sec.rho * sec.A]
            if kind == 'consistent':
                args.append('-cMass')
        ops.element(*args)
    for nm in f.node_masses:
        ops.mass(nm.node, nm.mx, nm.my, nm.Iz)


def mac(v, u):
    return float((v @ u) ** 2 / ((v @ v) * (u @ u)))


def main():
    if ops is None:
        print("SKIPPED: 沒有安裝 openseespy (pip install openseespy), 略過與 OpenSeesPy 的交叉驗證")
        return
    N_MODES = 6
    for name, make in (("懸臂梁(8元素, 只留彎曲)", cantilever), ("斜桿+桁架混合剛架", mixed_frame)):
        for kind in ('lumped', 'consistent'):
            f = make()
            build_opensees(f, kind)
            lam_os = np.array(ops.eigen('-fullGenLapack', N_MODES))
            md = eigen(f, n_modes=N_MODES, mass=kind)
            e = float(np.max(np.abs(md.omega**2 / lam_os - 1.0)))
            ids = sorted(f.nodes)
            macs = []
            for j in range(N_MODES):
                v = np.array([[ops.nodeEigenvector(n, j + 1, d) for d in (1, 2, 3)] for n in ids]).ravel()
                u = np.array([md.node_shape(j, n) for n in ids]).ravel()
                macs.append(mac(v, u))
            print(f"[{name} / {kind}] ω² 最大相對誤差 {e:.2e}, 最小MAC {min(macs):.12f}")
            assert e < 1e-9, "特徵值與 OpenSeesPy 不一致"
            assert min(macs) > 1 - 1e-9, "模態形狀與 OpenSeesPy 不一致(MAC)"

            mp = ops.modalProperties('-return')
            for d, key in (('x', 'MX'), ('y', 'MY')):
                tm = mp['totalMass'][0 if d == 'x' else 1]
                assert abs(md.mass_total[d] - tm) < 1e-10 * tm, f"總質量與OpenSees不一致({d})"
                if md.mass_free[d] <= 0:            # 該方向沒有可動質量(例如全部被拘束), 無參與係數可比
                    continue
                gam_os = np.abs(np.array(mp['partiFactor' + key])[:N_MODES])
                if kind == 'lumped':
                    free_os = mp['totalFreeMass'][0 if d == 'x' else 1]
                    assert abs(md.mass_free[d] - free_os) < 1e-10 * free_os, f"可動質量與OpenSees不一致({d})"
                    assert np.allclose(np.abs(md.gamma[d]), gam_os, rtol=1e-8, atol=1e-10), f"|Γ{d}|與OpenSees不一致"
                    ms_os = np.array(mp['partiMass' + key])[:N_MODES]
                    assert np.allclose(md.eff_mass[d], ms_os, rtol=1e-8, atol=1e-10), f"有效質量({d})與OpenSees不一致"
                else:
                    # 一致質量: 定義不同, 只檢查第1模態在 2% 內, 並印出差異供參考
                    d1 = abs(abs(md.gamma[d][0]) - gam_os[0]) / gam_os[0]
                    print(f"    一致質量 方向{d}: 第1模態 |Γ| frame2d {abs(md.gamma[d][0]):.5f} vs OpenSees {gam_os[0]:.5f} (差 {d1:.2%}, 定義不同)")
                    assert d1 < 0.02, "一致質量的第1模態參與係數與OpenSees相差超過2%, 不只是定義差異"

    print("\n全部通過: 特徵值(1e-9)、模態形狀MAC(1-1e-9)、總質量與 OpenSeesPy 一致; 集中質量的可動質量/參與係數/有效質量逐項一致; 一致質量的參與係數定義不同(差異<2%)。")


if __name__ == "__main__":
    main()
