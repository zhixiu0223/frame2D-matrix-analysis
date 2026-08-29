from .model import Frame2D, Node, Section, Member, Support, PointLoad, DistributedLoad
from .result import SolveResult, MemberResult
from .solve import solve_condensation
from .dofmanager import solve_dofmanager as solve   # 順序很重要: 見下方note

__all__ = [
    "Frame2D", "Node", "Section", "Member", "Support", "PointLoad", "DistributedLoad",
    "solve", "solve_condensation", "SolveResult", "MemberResult",
]

# note: `from .solve import ...` 這一行本身有個Python import機制的副作用
# ——會把套件本身的 `solve` 屬性設成「solve.py這個子模組物件」, 蓋掉手動
# 指定的名字綁定。所以`from .dofmanager import solve_dofmanager as solve`
# 這行一定要排在`from .solve import ...`後面, 才能讓solve_dofmanager
# 真正贏過那個子模組副作用, 變成frame2d.solve指向的東西。