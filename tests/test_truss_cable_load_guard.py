"""
驗證: truss/cable桿件加桿件內部載重(distributed_load/member_point_load)
必須明確報錯, 不能被靜默接受算出錯誤答案

背景: 使用者問「truss/cable上加桿件內部載重, 是真實的bug還是本來就不
給加」。查證後發現: solve.py(靜力凝縮版本)本來就有明確擋掉這個組合的
ValueError, 但dofmanager.py(主要求解器solve())沒有對應的檢查——只在
truss/cable桿件的"兩端都是外部支承、旋轉自由度整個懸空"這種特殊情況
才會被更通用的「自由度沒有勁度卻被施加外力」安全網意外攔下來; 一旦
truss/cable桿件的兩端都接到其他frame桿件(旋轉自由度被"借用", 不再
懸空), 均佈載重/桿件內部集中力的固定端彎矩項就會被靜默塞進求解,
算出一個看起來合理、但truss/cable桿端彎矩M1/M2不等於0(明確違反truss/
cable定義, 兩端絕對不可能有彎矩)的錯誤答案。

**這是真實的bug, 已經修正**(在dofmanager.py補上跟solve.py同樣的
ValueError檢查, 直接在載重組裝迴圈開頭擋掉)。這裡的測試專門針對
「兩端都接到其他frame桿件」這個容易被忽略、原本會繞過安全網的情況,
而不是重測"兩端都是外部支承"這種本來就會被攔下來的簡單情況。
"""
from frame2d import Frame2D, solve

E, A, I = 200e6, 0.01, 8e-5


def build_truss_between_frames(load_kind, member_type):
    """node0-node1: frame; node1-node2: truss或cable(在這裡加載重);
    node2-node3: frame。兩端都接frame, 旋轉自由度不懸空。"""
    f = Frame2D()
    f.add_node(0, 0, 0)
    f.add_node(1, 6, 0)
    f.add_node(2, 12, 0)
    f.add_node(3, 12, 4)
    f.add_section('s', E=E, I=I, A=A)
    f.add_member(0, 0, 1, 's')
    if member_type == 'truss':
        f.add_truss(1, 1, 2, 's')
    else:
        f.add_cable(1, 1, 2, 's')
    f.add_member(2, 2, 3, 's')
    f.fix(0)
    f.fix(3)
    if load_kind == 'distributed':
        f.distributed_load(1, w=-5.0)
    else:
        f.member_point_load(1, a=3.0, fy=-10.0)
    return f


print("=== 雙端接frame的truss/cable, 加桿件內部載重必須明確報錯 ===")
for member_type in ['truss', 'cable']:
    for load_kind in ['distributed', 'point']:
        f = build_truss_between_frames(load_kind, member_type)
        try:
            solve(f)
            raise AssertionError(
                f"{member_type}+{load_kind}: 應該要報錯, 但沒有(bug還在, 會靜默"
                "算出truss/cable桿端有非零彎矩的錯誤答案)")
        except ValueError as e:
            assert member_type in str(e), f"錯誤訊息應該提到{member_type}"
            print(f"  {member_type} + {load_kind}載重: 正確報錯 ✓")

print("\nPASS: dofmanager.py正確擋掉truss/cable+桿件內部載重的無效組合,")
print("不會再靜默算出truss/cable桿端彎矩不等於0的錯誤答案")
