"""
驗證案例30: distributed_load(direction='global_y') -- 全域垂直方向均佈載重
(屋頂重力/雪載重的標準表示方式), 對照SW FEA slop-roof案例

背景: 使用者的斜屋頂模型, 均佈載重(10kN/m)的方向是"沿全域垂直方向,
大小以沿桿件長度量測", 不是我們原本distributed_load()預設的"沿桿件
局部+y方向"(對斜桿件而言, 局部y是垂直於桿件本身, 不是真正的鉛直向下)。
這是README裡"分佈載重方向慣例"章節記錄的"尚未支援: 均佈載重的任意角度"
——這次有了具體案例驅動, 正式實作。

**設計**: distributed_load(member, w, direction='global_y')。全域載重
向量固定是(0,-w)(垂直向下, 大小w), 依桿件角度用transformation_matrix
分解成局部x(軸向)+局部y(橫向)兩個分量, 局部y分量套用既有的
fixed_end_forces_udl(), 局部x分量套用新增的fixed_end_forces_axial_udl()
(這條公式是新推導的, 用對稱結構的基本靜力學驗證: 對稱桿件受均佈軸向
載重, 兩端固定端反力必然對半分, 不需要額外的sympy推導)。兩者的固定端
反力向量直接相加。

目前只支援均佈(w_start=w_end)+整根桿件, 不支援局部段/線性變化
(model.py的DistributedLoad.__post_init__會擋掉其他組合, 需要時再擴充)。
只有主要求解器solve()(=DOFManager)支援, 靜力凝縮版本solve_condensation()
遇到這個方向會直接報錯(還沒實作)。

**開發過程中抓到的兩個bug(已修正)**:
1. postprocess.py的軸向累加項一開始正負號反了(N += 應該是 N -=), 用
   桿件實際的N(x)斜率(對照w_local_x的理論值)才抓出來。
2. 第一次逐點比對時用陣列位置直接對應, 因為桿件內部有均佈載重時
   member_internal_forces()會自動插入額外取樣點導致陣列長度變成13不是
   11, 造成離譜的誤差假象——這是這個專案裡第三次踩到同一種"取樣點位置
   vs陣列索引"的坑(前兩次分別在partial_udl跟combined_release_loads的
   測試裡出現過), 改用內插對齊(np.interp)才是正確做法。

**軸力正負號的note**: SW FEA報告的AxF欄位, 在這個案例裡是"壓力為正"
(跟我們"拉力為正"的慣例相反)——用最簡單的獨立案例(垂直懸臂柱, 頂端加
向下點載重, 物理上必然受壓)驗證過我們自己的N=-Fx1公式是正確的(給出
負值=受壓, 符合物理直覺), 所以這裡的比對統一用-N(我們的)去對SW的AxF。
"""
import numpy as np
from frame2d import Frame2D, solve, solve_condensation
from frame2d.postprocess import member_internal_forces

E, I, A = 200e6, 8e-5, 0.01


def build():
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 0, 4)
    f.add_node(3, 6, 4)
    f.add_node(4, 3, 7)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 2, 's')   # F0: 左柱
    f.add_member(1, 2, 4, 's')   # F7: 左屋頂斜梁
    f.add_member(2, 4, 3, 's')   # F8: 右屋頂斜梁
    f.add_member(3, 3, 1, 's')   # F2: 右柱
    f.fix(0)
    f.fix(1)
    f.distributed_load(1, w=10.0, direction='global_y')
    f.distributed_load(2, w=10.0, direction='global_y')
    return f


print("=== 反力比對 ===")
f = build()
r = solve(f)
sw_reactions = {0: (11.047, 42.426, -20.157), 1: (-11.047, 42.426, 20.157)}
for n in [0, 1]:
    ux, uy, rot = f.dofs_of(n)
    Rx, Ry, M = r.reactions[ux], r.reactions[uy], r.reactions[rot]
    swv = sw_reactions[n]
    print(f"  node{n}: f2d=({Rx:.3f},{Ry:.3f},{M:.3f})  sw={swv}")
    assert abs(Rx - swv[0]) < 0.01 and abs(Ry - swv[1]) < 0.01 and abs(M - swv[2]) < 0.01
print("  PASS: 反力精確吻合SW FEA\n")


print("=== 逐點N/V/M比對 (4桿件x11點, 用內插對齊) ===")
sw_data = {
    0: (4.0, [(42.426, -11.047, 20.157), (42.426, -11.047, 15.738), (42.426, -11.047, 11.319),
              (42.426, -11.047, 6.901), (42.426, -11.047, 2.482), (42.426, -11.047, -1.937),
              (42.426, -11.047, -6.355), (42.426, -11.047, -10.774), (42.426, -11.047, -15.193),
              (42.426, -11.047, -19.612), (42.426, -11.047, -24.030)]),
    1: (4.242640687119285, [(37.811, 22.189, -24.030), (34.811, 19.189, -15.253), (31.811, 16.189, -7.748),
                             (28.811, 13.189, -1.516), (25.811, 10.189, 3.443), (22.811, 7.189, 7.129),
                             (19.811, 4.189, 9.543), (16.811, 1.189, 10.684), (13.811, -1.811, 10.551),
                             (10.811, -4.811, 9.147), (7.811, -7.811, 6.469)]),
    2: (4.242640687119285, [(7.811, 7.811, 6.469), (10.811, 4.811, 9.147), (13.811, 1.811, 10.551),
                             (16.811, -1.189, 10.684), (19.811, -4.189, 9.543), (22.811, -7.189, 7.129),
                             (25.811, -10.189, 3.443), (28.811, -13.189, -1.516), (31.811, -16.189, -7.748),
                             (34.811, -19.189, -15.253), (37.811, -22.189, -24.030)]),
    3: (4.0, [(42.426, 11.047, -24.030), (42.426, 11.047, -19.612), (42.426, 11.047, -15.193),
              (42.426, 11.047, -10.774), (42.426, 11.047, -6.355), (42.426, 11.047, -1.937),
              (42.426, 11.047, 2.482), (42.426, 11.047, 6.901), (42.426, 11.047, 11.319),
              (42.426, 11.047, 15.738), (42.426, 11.047, 20.157)]),
}
max_err = 0
for mid, (L, rows) in sw_data.items():
    x_full, N, V, M = member_internal_forces(f, r, mid, n=201)
    x_query = np.linspace(0, L, 11)
    N_i = np.interp(x_query, x_full, N)
    V_i = np.interp(x_query, x_full, V)
    M_i = np.interp(x_query, x_full, M)
    for i, (n_sw, v_sw, m_sw) in enumerate(rows):
        max_err = max(max_err, abs(-N_i[i] - n_sw), abs(V_i[i] - v_sw), abs(M_i[i] - m_sw))
print(f"  整體最大誤差: {max_err:.4f}")
assert max_err < 0.01, f"誤差過大: {max_err}"
print("  PASS: 4桿件x11點N/V/M全數吻合SW FEA(誤差在報告小數位精度範圍內)\n")


print("=== 靜力凝縮版本目前還不支援這個方向, 應該直接報錯 ===")
try:
    solve_condensation(f)
    raise AssertionError("預期應該報錯, 但沒有")
except ValueError:
    print("  確認: 報錯符合預期(已知限制, solve()才支援)\n")

print("PASS: distributed_load(direction='global_y') 驗證完成")
