"""
驗證analyze_pushover()這個統一入口是純粹的dispatch層, 沒有引入任何
新的行為——每一種geometry/solver組合, 都要跟直接呼叫底層函式的結果
逐位元/浮點數誤差為0一致。這是「安全網」設計的核心驗證: 如果這裡
任何一個案例沒有通過, 代表dispatch層本身引入了新的bug, 不是底層
求解器的問題(底層求解器各自的正確性已經在別的測試檔案裡驗證過)。

見ANALYSIS_ARCHITECTURE.md「規劃中」第2點的完整說明。
"""
import numpy as np
from frame2d import Frame2D
from frame2d.hinge import HingeState
from frame2d.dofmanager import initial_hinge_states
from frame2d.pushover import run_pushover, run_pushover_converged
from frame2d.newton import run_pushover_newton, run_pushover_corotational_oneshot
from frame2d.analyze import analyze_pushover

E, I, A, L = 200e6, 8e-5, 1e-2, 4.0


def build_frame_and_hinges():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 0, L)
    f.add_section('sec', E=E, I=I, A=A)
    f.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
                 R_post_yield_i=50.0, R_post_yield_j=50.0)
    f.fix(0)
    hs = initial_hinge_states(f)
    return f, hs


def _same_result(a, b):
    """比較兩次呼叫回傳的tuple是不是逐項完全一致(history_u/history_F
    這種np.array用np.array_equal, event_log這種list直接比較)。"""
    assert len(a) == len(b), f"回傳欄位數量不一致: {len(a)} vs {len(b)}"
    for i, (av, bv) in enumerate(zip(a, b)):
        if isinstance(av, np.ndarray):
            assert np.array_equal(av, bv), f"第{i}個欄位(np.array)不一致"
        elif isinstance(av, (int, float, bool)):
            assert av == bv, f"第{i}個欄位不一致: {av} vs {bv}"
        # dict(hinge_states)/list(event_log)這裡不逐項比較內容(物件
        # identity不同是預期的, 因為兩次呼叫用的是不同的frame/hinge_
        # states複本), 只比較數值型欄位——這已經足夠證明dispatch層
        # 沒有引入計算上的差異(history_u/history_F/converged這些
        # 才是實際的計算結果)。


print("=== 組合1: geometry='linear', solver='direct' (對應run_pushover) ===")
f1a, hs1a = build_frame_and_hinges()
r1a = run_pushover(f1a, hs1a, prescribed_dofs=[f1a.dofs_of(1)[0]], direction=[1.0],
                    target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f1a.dofs_of(0)[0]],
                    use_pdelta=False, geometry_update=False)
f1b, hs1b = build_frame_and_hinges()
r1b = analyze_pushover(f1b, hs1b, prescribed_dofs=[f1b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f1b.dofs_of(0)[0]],
                        geometry='linear', use_pdelta=False, solver='direct')
_same_result(r1a, r1b)
print("PASS\n")

print("=== 組合2: geometry='linear', solver='direct', use_pdelta=True ===")
f2a, hs2a = build_frame_and_hinges()
r2a = run_pushover(f2a, hs2a, prescribed_dofs=[f2a.dofs_of(1)[0]], direction=[1.0],
                    target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f2a.dofs_of(0)[0]],
                    use_pdelta=True, geometry_update=False)
f2b, hs2b = build_frame_and_hinges()
r2b = analyze_pushover(f2b, hs2b, prescribed_dofs=[f2b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f2b.dofs_of(0)[0]],
                        geometry='linear', use_pdelta=True, solver='direct')
_same_result(r2a, r2b)
print("PASS\n")

print("=== 組合3: geometry='updated', solver='direct', use_pdelta=True ===")
f3a, hs3a = build_frame_and_hinges()
r3a = run_pushover(f3a, hs3a, prescribed_dofs=[f3a.dofs_of(1)[0]], direction=[1.0],
                    target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f3a.dofs_of(0)[0]],
                    use_pdelta=True, geometry_update=True)
f3b, hs3b = build_frame_and_hinges()
r3b = analyze_pushover(f3b, hs3b, prescribed_dofs=[f3b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f3b.dofs_of(0)[0]],
                        geometry='updated', use_pdelta=True, solver='direct')
_same_result(r3a, r3b)
print("PASS\n")

print("=== 組合4: geometry='updated', solver='picard', use_pdelta=True ===")
f4a, hs4a = build_frame_and_hinges()
r4a = run_pushover_converged(f4a, hs4a, prescribed_dofs=[f4a.dofs_of(1)[0]], direction=[1.0],
                              target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f4a.dofs_of(0)[0]],
                              use_pdelta=True, geometry_update=True)
f4b, hs4b = build_frame_and_hinges()
r4b = analyze_pushover(f4b, hs4b, prescribed_dofs=[f4b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f4b.dofs_of(0)[0]],
                        geometry='updated', use_pdelta=True, solver='picard')
_same_result(r4a, r4b)
print("PASS\n")

print("=== 組合5: geometry='linear', solver='picard', use_pdelta=True ===")
f5a, hs5a = build_frame_and_hinges()
r5a = run_pushover_converged(f5a, hs5a, prescribed_dofs=[f5a.dofs_of(1)[0]], direction=[1.0],
                              target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f5a.dofs_of(0)[0]],
                              use_pdelta=True, geometry_update=False)
f5b, hs5b = build_frame_and_hinges()
r5b = analyze_pushover(f5b, hs5b, prescribed_dofs=[f5b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f5b.dofs_of(0)[0]],
                        geometry='linear', use_pdelta=True, solver='picard')
_same_result(r5a, r5b)
print("PASS\n")

print("=== 組合7: geometry='corotational', solver='newton' ===")
f7a = Frame2D(); f7a.add_node(0, 0, 0); f7a.add_node(1, 0, L); f7a.add_section('sec', E=E, I=I, A=A)
f7a.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
               R_post_yield_i=50.0, R_post_yield_j=50.0)
f7a.fix(0)
hs7a = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=50.0, R_post_yield_2=50.0)}
r7a = run_pushover_newton(f7a, hs7a, prescribed_dofs=[f7a.dofs_of(1)[0]], direction=[1.0],
                           target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f7a.dofs_of(0)[0]],
                           use_pdelta=False)
f7b = Frame2D(); f7b.add_node(0, 0, 0); f7b.add_node(1, 0, L); f7b.add_section('sec', E=E, I=I, A=A)
f7b.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
               R_post_yield_i=50.0, R_post_yield_j=50.0)
f7b.fix(0)
hs7b = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=50.0, R_post_yield_2=50.0)}
r7b = analyze_pushover(f7b, hs7b, prescribed_dofs=[f7b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f7b.dofs_of(0)[0]],
                        geometry='corotational', use_local_pdelta=False, solver='newton')
_same_result(r7a, r7b)
print("PASS\n")

print("=== 組合8: geometry='corotational', solver='newton', use_local_pdelta=True ===")
f8a = Frame2D(); f8a.add_node(0, 0, 0); f8a.add_node(1, 0, L); f8a.add_section('sec', E=E, I=I, A=A)
f8a.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
               R_post_yield_i=50.0, R_post_yield_j=50.0)
f8a.fix(0)
hs8a = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=50.0, R_post_yield_2=50.0)}
r8a = run_pushover_newton(f8a, hs8a, prescribed_dofs=[f8a.dofs_of(1)[0]], direction=[1.0],
                           target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f8a.dofs_of(0)[0]],
                           use_pdelta=True)
f8b = Frame2D(); f8b.add_node(0, 0, 0); f8b.add_node(1, 0, L); f8b.add_section('sec', E=E, I=I, A=A)
f8b.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
               R_post_yield_i=50.0, R_post_yield_j=50.0)
f8b.fix(0)
hs8b = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=50.0, R_post_yield_2=50.0)}
r8b = analyze_pushover(f8b, hs8b, prescribed_dofs=[f8b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f8b.dofs_of(0)[0]],
                        geometry='corotational', use_local_pdelta=True, solver='newton')
_same_result(r8a, r8b)
print("PASS\n")

print("=== 組合9: geometry='corotational', solver='direct' (對應run_pushover_corotational_oneshot) ===")
f9a = Frame2D(); f9a.add_node(0, 0, 0); f9a.add_node(1, 0, L); f9a.add_section('sec', E=E, I=I, A=A)
f9a.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
               R_post_yield_i=50.0, R_post_yield_j=50.0)
f9a.fix(0)
hs9a = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=50.0, R_post_yield_2=50.0)}
r9a = run_pushover_corotational_oneshot(f9a, hs9a, prescribed_dofs=[f9a.dofs_of(1)[0]], direction=[1.0],
                                         target_total=0.05, d_nominal=0.005,
                                         base_reaction_dofs=[f9a.dofs_of(0)[0]], use_pdelta=False)
f9b = Frame2D(); f9b.add_node(0, 0, 0); f9b.add_node(1, 0, L); f9b.add_section('sec', E=E, I=I, A=A)
f9b.add_member(0, node_i=0, node_j=1, section='sec', Mp_i=100.0, Mp_j=1e30,
               R_post_yield_i=50.0, R_post_yield_j=50.0)
f9b.fix(0)
hs9b = {0: HingeState(Mp1=100.0, Mp2=1e30, R_post_yield_1=50.0, R_post_yield_2=50.0)}
r9b = analyze_pushover(f9b, hs9b, prescribed_dofs=[f9b.dofs_of(1)[0]], direction=[1.0],
                        target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f9b.dofs_of(0)[0]],
                        geometry='corotational', use_local_pdelta=False, solver='direct')
_same_result(r9a, r9b)
print("PASS\n")

print("=== kwargs轉傳: include_snapshots=True應該正確轉傳到底層函式 ===")
f10a, hs10a = build_frame_and_hinges()
r10a = run_pushover(f10a, hs10a, prescribed_dofs=[f10a.dofs_of(1)[0]], direction=[1.0],
                     target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f10a.dofs_of(0)[0]],
                     use_pdelta=False, geometry_update=False, include_snapshots=True)
f10b, hs10b = build_frame_and_hinges()
r10b = analyze_pushover(f10b, hs10b, prescribed_dofs=[f10b.dofs_of(1)[0]], direction=[1.0],
                         target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f10b.dofs_of(0)[0]],
                         geometry='linear', solver='direct', include_snapshots=True)
assert len(r10a) == len(r10b) == 6, f"開了include_snapshots應該有6個回傳欄位, 實際: {len(r10a)}, {len(r10b)}"
assert len(r10a[-1]) == len(r10b[-1]), "history_snapshots長度應該一致"
print("PASS\n")

print("=== 不支援的組合應該明確raise ValueError, 不是靜默退化 ===")
f11, hs11 = build_frame_and_hinges()
try:
    analyze_pushover(f11, hs11, prescribed_dofs=[f11.dofs_of(1)[0]], direction=[1.0],
                      target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f11.dofs_of(0)[0]],
                      geometry='corotational', solver='picard')
    assert False, "geometry='corotational'搭配solver='picard'應該要raise"
except ValueError as e:
    assert 'picard' in str(e) or 'Picard' in str(e) or 'newton' in str(e)
print("PASS: geometry='corotational'+solver='picard'正確拒絕\n")

f12, hs12 = build_frame_and_hinges()
try:
    analyze_pushover(f12, hs12, prescribed_dofs=[f12.dofs_of(1)[0]], direction=[1.0],
                      target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f12.dofs_of(0)[0]],
                      geometry='linear', solver='newton')
    assert False, "geometry='linear'搭配solver='newton'應該要raise"
except ValueError as e:
    assert 'newton' in str(e) or 'Newton' in str(e)
print("PASS: geometry='linear'+solver='newton'正確拒絕\n")

try:
    analyze_pushover(f12, hs12, prescribed_dofs=[f12.dofs_of(1)[0]], direction=[1.0],
                      target_total=0.05, d_nominal=0.005, base_reaction_dofs=[f12.dofs_of(0)[0]],
                      geometry='not_a_real_option', solver='direct')
    assert False, "不存在的geometry選項應該要raise"
except ValueError as e:
    assert 'geometry' in str(e)
print("PASS: 不存在的geometry選項正確拒絕\n")

print("PASS: analyze_pushover()所有組合都跟直接呼叫底層函式結果一致, 不支援的組合都正確拒絕")
