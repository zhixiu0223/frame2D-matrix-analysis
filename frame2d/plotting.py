"""
繪圖模組: 結構圖 / 受力圖 / 軸力圖 / 剪力圖 / 彎矩圖 / 變形圖

只負責畫圖, 不重算任何力學 —— 所有數字都來自 SolveResult 跟 postprocess.py。
在 Pydroid3 上呼叫 plt.show() 會跳出內建繪圖檢視畫面; 在無顯示器環境
(Termux/伺服器)請改用 plt.savefig(...)。
"""
import numpy as np
import matplotlib.pyplot as plt

from .postprocess import member_internal_forces, member_deformed_shape
from .elements import member_geometry


def _member_endpoints(frame, member_id):
    m = frame.members[member_id]
    ni, nj = frame.nodes[m.node_i], frame.nodes[m.node_j]
    return ni, nj


def _draw_moment_arc(ax, x, y, m_value, radius, color):
    """在(x,y)畫一個帶箭頭的弧形, 標示集中力矩的旋轉方向:
    m>0(逆時針, CCW)畫↺, m<0(順時針, CW)畫↻。用matplotlib的FancyArrowPatch
    畫270度弧線+箭頭, 弧線方向依m的正負決定。"""
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.path import Path
    theta1, theta2 = (-60, 210) if m_value > 0 else (210, -60)
    n_pts = 30
    thetas = np.linspace(np.radians(theta1), np.radians(theta2), n_pts)
    arc_x = x + radius * np.cos(thetas)
    arc_y = y + radius * np.sin(thetas)
    ax.plot(arc_x, arc_y, color=color, lw=1.5, zorder=3)
    # 箭頭畫在弧線末端, 方向沿著弧線切線
    tip_x, tip_y = arc_x[-1], arc_y[-1]
    prev_x, prev_y = arc_x[-2], arc_y[-2]
    ax.annotate('', xy=(tip_x, tip_y), xytext=(prev_x, prev_y),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=1.5, mutation_scale=10), zorder=3)


def plot_structure(frame, ax=None, show_node_ids=True, show_member_ids=True, show_dimensions=True):
    """① 結構尺寸圖 (可標註每根桿件長度 + 結構總寬/總高)"""
    if ax is None:
        fig, ax = plt.subplots()
    for mid, m in frame.members.items():
        ni, nj = _member_endpoints(frame, mid)
        if m.member_type == 'truss':
            ax.plot([ni.x, nj.x], [ni.y, nj.y], 'k--', lw=1.3, zorder=1)      # 桁架: 黑虛線
        elif m.member_type == 'cable':
            ax.plot([ni.x, nj.x], [ni.y, nj.y], color='teal', ls=':', lw=1.3, zorder=1)  # 纜線: 藍綠點線
        else:
            ax.plot([ni.x, nj.x], [ni.y, nj.y], 'k-', lw=2, zorder=1)        # 樑柱: 實線
        if show_member_ids:
            L, _ = member_geometry(ni, nj)
            xm, ym = (ni.x + nj.x) / 2, (ni.y + nj.y) / 2
            tag = {'frame': f'M{mid}', 'truss': f'T{mid}', 'cable': f'C{mid}'}[m.member_type]
            label = f'{tag} (L={L:.2f})' if show_dimensions else tag
            ax.annotate(label, (xm, ym), color='blue', fontsize=8,
                        ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
        # 內部鉸接(release): 在釋放的那一端畫個空心紅圈標記(跟SW FEA app的
        # 慣例一致), offset一點點避免蓋住節點本身的黑點
        if getattr(m, 'release_i', False) or getattr(m, 'release_j', False):
            L, angle = member_geometry(ni, nj)
            c, s = np.cos(angle), np.sin(angle)
            off = min(L * 0.06, 0.15)
            if m.release_i:
                ax.plot(ni.x + off * c, ni.y + off * s, 'o', mfc='none', mec='red', mew=1.3, ms=7, zorder=3)
            if m.release_j:
                ax.plot(nj.x - off * c, nj.y - off * s, 'o', mfc='none', mec='red', mew=1.3, ms=7, zorder=3)
    for nid, n in frame.nodes.items():
        ax.plot(n.x, n.y, 'ko', ms=5, zorder=2)
        if show_node_ids:
            ax.annotate(f'N{nid}', (n.x, n.y), textcoords="offset points",
                        xytext=(6, 6), fontsize=8, color='dimgray')

    for s in frame.supports:
        n = frame.nodes[s.node]
        _draw_support_symbol(ax, n.x, n.y, s, _auto_scale(frame) * 0.05)

    if show_dimensions:
        _draw_overall_dimensions(frame, ax)

    ax.set_aspect('equal')
    ax.set_title('Geometry')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(True, alpha=0.3)
    return ax


def _draw_support_symbol(ax, x, y, support, size):
    """依工程圖標準慣例區分三種支承符號 (假設支承畫在結構下方,排線朝下):
    - Fixed (ux,uy,rot全拘束): 橫線 + 排線(牆面) — 不是三角形, 避免跟鉸接搞混
    - Hinge/Pin (ux,uy拘束,rot自由): 三角形 + 排線
    - Roller (只拘束1個平移方向): 圓圈 + 排線
    """
    n_hatch = 5
    hatch_y0 = y - size
    rot_c = support.rot is not None
    ux_c = support.ux is not None
    uy_c = support.uy is not None

    def _draw_hatch(base_y):
        ax.plot([x - size, x + size], [base_y, base_y], 'k-', lw=1.3, zorder=1)
        for i in range(n_hatch + 1):
            hx = x - size + 2 * size * i / n_hatch
            ax.plot([hx, hx - size * 0.4], [base_y, base_y - size * 0.5], 'k-', lw=0.8, zorder=1)

    if rot_c and ux_c and uy_c:
        # Fixed: 直接在節點位置畫橫線+排線(牆面), 沒有三角形
        _draw_hatch(y)
    elif ux_c and uy_c and not rot_c:
        # Hinge/Pin: 三角形(頂點在節點) + 排線
        ax.plot([x - size, x, x + size, x - size], [hatch_y0, y, hatch_y0, hatch_y0],
                'k-', lw=1.3, zorder=1)
        _draw_hatch(hatch_y0)
    else:
        # Roller: 三角形 + 圓圈 + 排線 (只拘束1個方向)
        ax.plot([x - size, x, x + size, x - size], [hatch_y0, y, hatch_y0, hatch_y0],
                'k-', lw=1.3, zorder=1)
        circ = plt.Circle((x, hatch_y0 - size * 0.3), size * 0.3, fill=False, color='k', lw=1.0, zorder=1)
        ax.add_patch(circ)
        _draw_hatch(hatch_y0 - size * 0.6)


def _draw_overall_dimensions(frame, ax):
    """在結構左側/下方畫總高/總寬的尺寸線(雙箭頭+數字)"""
    xs = [n.x for n in frame.nodes.values()]
    ys = [n.y for n in frame.nodes.values()]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    scale = _auto_scale(frame)
    off_x = x_min - scale * 0.12
    off_y = y_min - scale * 0.12

    if y_max - y_min > 1e-9:
        ax.annotate('', xy=(off_x, y_max), xytext=(off_x, y_min),
                    arrowprops=dict(arrowstyle='<->', color='gray', lw=1))
        ax.plot([off_x, x_min], [y_min, y_min], color='gray', lw=0.5, ls=':')
        ax.plot([off_x, x_min], [y_max, y_max], color='gray', lw=0.5, ls=':')
        ax.annotate(f'H={y_max - y_min:.2f}', (off_x, (y_min + y_max) / 2),
                    ha='right', va='center', fontsize=8, color='gray', rotation=90)

    if x_max - x_min > 1e-9:
        ax.annotate('', xy=(x_max, off_y), xytext=(x_min, off_y),
                    arrowprops=dict(arrowstyle='<->', color='gray', lw=1))
        ax.plot([x_min, x_min], [off_y, y_min], color='gray', lw=0.5, ls=':')
        ax.plot([x_max, x_max], [off_y, y_min], color='gray', lw=0.5, ls=':')
        ax.annotate(f'L={x_max - x_min:.2f}', ((x_min + x_max) / 2, off_y),
                    ha='center', va='top', fontsize=8, color='gray')


def _auto_scale(frame):
    xs = [n.x for n in frame.nodes.values()]
    ys = [n.y for n in frame.nodes.values()]
    return max(max(xs) - min(xs), max(ys) - min(ys), 1.0)


def plot_loads(frame, ax=None):
    """② 受力圖"""
    if ax is None:
        fig, ax = plt.subplots()
    plot_structure(frame, ax=ax, show_node_ids=False, show_member_ids=False, show_dimensions=False)
    scale = _auto_scale(frame)
    extra_pts_x, extra_pts_y = [], []

    for pl in frame.point_loads:
        n = frame.nodes[pl.node]
        mag = np.hypot(pl.fx, pl.fy)
        if mag > 1e-9:
            L = scale * 0.15
            dx, dy = pl.fx / mag * L, pl.fy / mag * L
            tail_x, tail_y = n.x - dx, n.y - dy
            ax.annotate('', xy=(n.x, n.y), xytext=(tail_x, tail_y),
                        arrowprops=dict(arrowstyle='->', color='red', lw=2))
            parts = []
            if abs(pl.fx) > 1e-9:
                parts.append(f'Fx={pl.fx:.3g} N')
            if abs(pl.fy) > 1e-9:
                parts.append(f'Fy={pl.fy:.3g} N')
            ax.annotate(', '.join(parts), (tail_x, tail_y), color='red', fontsize=8,
                        ha='center', xytext=(0, -10 if dy >= 0 else 10), textcoords='offset points')
            extra_pts_x += [tail_x]
            extra_pts_y += [tail_y]
        if abs(pl.m) > 1e-9:
            arc_r = scale * 0.06
            _draw_moment_arc(ax, n.x, n.y, pl.m, arc_r, 'purple')
            ax.annotate(f'M={pl.m:.3g} N*m', (n.x, n.y), color='purple', fontsize=8,
                        xytext=(8 + arc_r * 60, -12), textcoords='offset points')

    for dl in frame.distributed_loads:
        ni, nj = _member_endpoints(frame, dl.member)
        L, angle = member_geometry(ni, nj)
        c, s = np.cos(angle), np.sin(angle)
        x_start = 0.0 if dl.x_start is None else dl.x_start
        x_end = L if dl.x_end is None else dl.x_end
        n_arrows = 6
        max_w = max(abs(dl.w_start), abs(dl.w_end), 1e-9)
        arrow_len = scale * 0.08
        # 箭頭方向: direction='local'(預設)畫在垂直於桿件的方向(局部+y);
        # direction='global_y'或'global'畫成真正的全域方向, 不管桿件本身
        # 斜不斜, 才不會誤導成"垂直於斜屋頂表面"
        if dl.direction == 'global_y':
            perp_x, perp_y = 0.0, -1.0   # 全域垂直向下, 固定方向不隨桿件角度變
        elif dl.direction == 'global':
            ang = np.radians(dl.angle_deg)
            perp_x, perp_y = np.cos(ang), np.sin(ang)   # 全域任意角度, 固定方向不隨桿件角度變
        else:
            perp_x, perp_y = -s, c       # 局部+y方向(垂直於桿件)
        for k in range(n_arrows + 1):
            t = k / n_arrows
            wx = dl.w_start + (dl.w_end - dl.w_start) * t
            xt = x_start + (x_end - x_start) * t   # 只在[x_start,x_end]範圍內畫箭頭
            px, py = ni.x + xt * c, ni.y + xt * s
            length = arrow_len * (wx / max_w) if max_w > 0 else 0
            tail_x, tail_y = px - perp_x * length, py - perp_y * length
            ax.annotate('', xy=(px, py), xytext=(tail_x, tail_y),
                        arrowprops=dict(arrowstyle='->', color='orange', lw=1))
            extra_pts_x += [tail_x]
            extra_pts_y += [tail_y]
        if abs(dl.w_start - dl.w_end) < 1e-9:
            w_label = f'w={dl.w_start:.3g} N/m'
        else:
            w_label = f'w={dl.w_start:.3g}~{dl.w_end:.3g} N/m'
        mx = ni.x + (x_start + x_end) / 2 * c
        my = ni.y + (x_start + x_end) / 2 * s
        ax.annotate(w_label, (mx, my), color='orange', fontsize=8,
                    xytext=(0, 10), textcoords='offset points', ha='center')

    for pl in frame.member_point_loads:
        ni, nj = _member_endpoints(frame, pl.member)
        L, angle = member_geometry(ni, nj)
        c, s = np.cos(angle), np.sin(angle)
        px, py = ni.x + pl.a * c, ni.y + pl.a * s   # 桿件內部位置(局部座標a換算成全域座標)
        if pl.direction == 'global':
            # direction='global': F+angle_deg直接就是全域方向, 不用像
            # local模式那樣把局部fx/fy轉回全域(這裡沒有局部分量可轉,
            # pl.fx/pl.fy在這個模式下固定是0)
            mag = pl.F
            ang = np.radians(pl.angle_deg)
            gdx, gdy = np.cos(ang) * mag, np.sin(ang) * mag
            label = f'F={pl.F:.3g} N @{pl.angle_deg:.0f}°'
        else:
            mag = np.hypot(pl.fx, pl.fy)
            # fx,fy是局部座標分量, 換算成全域方向畫箭頭
            gdx = pl.fx * c - pl.fy * s
            gdy = pl.fx * s + pl.fy * c
            parts = []
            if abs(pl.fx) > 1e-9:
                parts.append(f'fx={pl.fx:.3g} N')
            if abs(pl.fy) > 1e-9:
                parts.append(f'fy={pl.fy:.3g} N')
            label = ', '.join(parts)
        if mag > 1e-9:
            Larrow = scale * 0.15
            dx, dy = gdx / mag * Larrow, gdy / mag * Larrow
            tail_x, tail_y = px - dx, py - dy
            ax.annotate('', xy=(px, py), xytext=(tail_x, tail_y),
                        arrowprops=dict(arrowstyle='->', color='crimson', lw=2))
            ax.annotate(label, (tail_x, tail_y), color='crimson', fontsize=8,
                        ha='center', xytext=(0, -10 if dy >= 0 else 10), textcoords='offset points')
            extra_pts_x += [tail_x]
            extra_pts_y += [tail_y]
        if abs(pl.m) > 1e-9:
            arc_r = scale * 0.06
            ax.plot(px, py, 'D', color='darkmagenta', ms=5, zorder=2)
            _draw_moment_arc(ax, px, py, pl.m, arc_r, 'darkmagenta')
            ax.annotate(f'm={pl.m:.3g} N*m', (px, py), color='darkmagenta', fontsize=8,
                        xytext=(8 + arc_r * 60, -12), textcoords='offset points')

    # 直接手動計算範圍再 set_xlim/ylim, 比依賴 relim/autoscale 對 annotate 更可靠
    if extra_pts_x:
        xs = [n.x for n in frame.nodes.values()] + extra_pts_x
        ys = [n.y for n in frame.nodes.values()] + extra_pts_y
        pad_x = (max(xs) - min(xs)) * 0.1 + 1e-6
        pad_y = (max(ys) - min(ys)) * 0.1 + 1e-6
        ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
        ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)

    ax.set_title('Loads')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    return ax


def plot_diagram(frame, result, kind, ax=None, scale=None):
    """③④⑤ 軸力圖/剪力圖/彎矩圖. kind: 'N', 'V', 或 'M' """
    if ax is None:
        fig, ax = plt.subplots()
    plot_structure(frame, ax=ax, show_node_ids=False, show_member_ids=False, show_dimensions=False)

    idx = {'N': 1, 'V': 2, 'M': 3}[kind]
    label = {'N': 'Axial Force N (units: N)', 'V': 'Shear Force V (units: N)', 'M': 'Bending Moment M (units: N*m)'}[kind]

    n_samples = 41  # 加密取樣點, 讓極值偵測(轉折點)更準確
    # 先算全部member的最大值, 用來自動定比例尺
    all_vals = []
    per_member = {}
    for mid in frame.members:
        x, N, V, M = member_internal_forces(frame, result, mid, n=n_samples)
        vals = (N, V, M)[idx - 1]
        per_member[mid] = (x, vals)
        all_vals.append(np.max(np.abs(vals)))
    max_val = max(all_vals) if all_vals else 1.0
    if max_val < 1e-9:
        max_val = 1.0

    if scale is None:
        scale = _auto_scale(frame) * 0.15 / max_val

    for mid in frame.members:
        ni, nj = _member_endpoints(frame, mid)
        L, angle = member_geometry(ni, nj)
        c, s = np.cos(angle), np.sin(angle)
        x, vals = per_member[mid]
        perp_x, perp_y = -s, c
        px = ni.x + x * c - perp_x * vals * scale
        py = ni.y + x * s - perp_y * vals * scale
        base_x = ni.x + x * c
        base_y = ni.y + x * s
        ax.plot(px, py, color='crimson', lw=1)
        ax.fill(np.concatenate([base_x, px[::-1]]), np.concatenate([base_y, py[::-1]]),
                color='crimson', alpha=0.2)

        # 標示: 兩端點 + 桿件內部所有局部極值(轉折點, 即幾何不連續/曲線變號的地方)
        # 標籤往外(遠離桿件軸線的方向)偏移一點, 避免壓在曲線上被誤讀
        for i in _label_indices(vals):
            sign = 1 if vals[i] >= 0 else -1
            off_x = -perp_x * sign * 10
            off_y = -perp_y * sign * 10
            ax.annotate(f'{vals[i]:.3g}', (px[i], py[i]), fontsize=7, color='crimson',
                        ha='center', va='center', xytext=(off_x, off_y), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='crimson', lw=0.5, alpha=0.9))

    ax.set_title(label)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    return ax


def _label_indices(vals, atol=1e-6):
    """找出陣列裡該標數值的index: 兩端點 + 內部所有局部極值(轉折點)"""
    idxs = {0, len(vals) - 1}
    d = np.diff(vals)
    for i in range(1, len(d)):
        if d[i - 1] * d[i] < -atol:  # 相鄰斜率變號 -> 局部極大/極小
            idxs.add(i)
    return sorted(idxs)


def plot_deformed(frame, result, ax=None, scale=None, show_values=True):
    """⑥ 變形圖 (show_values=True: 在每個節點旁標出實際位移量, 方便對照SW FEA
    這類會直接印出位移數字的工具, 快速確認斷面/材料設定是否正確)"""
    if ax is None:
        fig, ax = plt.subplots()
    plot_structure(frame, ax=ax, show_node_ids=False, show_member_ids=False, show_dimensions=False)

    # 用實際內插出來的變形曲線(不是只看節點平移量)算最大偏移量, 才不會漏掉
    # "兩端節點位移都是0、但轉角很大"的情況(例如簡支梁純轉角撐起中間撓度)
    if scale is None:
        max_offset = 1e-9
        for mid in frame.members:
            ni, nj = _member_endpoints(frame, mid)
            L, angle = member_geometry(ni, nj)
            base_x = ni.x + np.linspace(0, L, 21) * np.cos(angle)
            base_y = ni.y + np.linspace(0, L, 21) * np.sin(angle)
            X, Y = member_deformed_shape(frame, result, mid, scale=1.0, n=21)
            offset = np.hypot(X - base_x, Y - base_y)
            max_offset = max(max_offset, np.max(offset))
        scale = _auto_scale(frame) * 0.1 / max_offset

    for mid in frame.members:
        X, Y = member_deformed_shape(frame, result, mid, scale=scale, n=41)
        ax.plot(X, Y, color='blue', lw=2)

    if show_values:
        for nid, n in frame.nodes.items():
            ux, uy, _ = frame.dofs_of(nid)
            dx, dy = result.displacements[ux], result.displacements[uy]
            mag = np.hypot(dx, dy)
            if mag < 1e-12:
                continue   # 支承等位移=0的節點不標, 避免畫面雜亂
            ax.annotate(f'Δ={mag:.4g} m', (n.x, n.y), fontsize=7, color='darkblue',
                        xytext=(5, -10), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='darkblue', lw=0.5, alpha=0.85))

        # 每根桿件內部(不只節點)的最大偏移量也標出來 -- 簡支梁這類兩端節點
        # 位移剛好=0(支承)、真正最大撓度發生在跨間的情況, 不然會漏標
        for mid in frame.members:
            ni, nj = _member_endpoints(frame, mid)
            L, angle = member_geometry(ni, nj)
            n_sample = 41
            t = np.linspace(0, L, n_sample)
            base_x = ni.x + t * np.cos(angle)
            base_y = ni.y + t * np.sin(angle)
            X1, Y1 = member_deformed_shape(frame, result, mid, scale=scale, n=n_sample)      # 放大過的(畫圖用)
            X1r, Y1r = member_deformed_shape(frame, result, mid, scale=1.0, n=n_sample)       # 真實值(標籤用)
            offset_real = np.hypot(X1r - base_x, Y1r - base_y)
            i_max = np.argmax(offset_real)
            if offset_real[i_max] < 1e-9 or not (0 < i_max < n_sample - 1):
                continue   # 端點已經在上面節點迴圈標過, 這裡只標"桿件內部"的最大值避免重複
            ax.annotate(f'Δmax={offset_real[i_max]:.4g} m', (X1[i_max], Y1[i_max]), fontsize=7, color='darkblue',
                        xytext=(5, -12), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='darkblue', lw=0.5, alpha=0.85))

    ax.set_title(f'Deformation (scale x{scale:.1f})')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    return ax


def _simpson(y, x):
    """辛普森法數值積分(需要偶數個區間, 即奇數個點)。之前這裡用手動
    梯形法, 但梯形法對「力矩=∫w(x)*x dx」這種積分不是精確的——w(x)
    本身是線性的梯形載重, x*w(x)是二次式, 梯形法對二次式有O(h²)的
    誤差, 實際測出來400個細分還是有0.06這種級別的殘留誤差(手算過
    梯形法誤差公式-  (b-a)³/(12n²)*f'', 算出來剛好等於實測的殘留量,
    確認就是這個原因, 不是邏輯寫錯)。辛普森法對三次以下多項式是
    精確積分, 能完全消除這個誤差來源, 不只是縮小它。"""
    y, x = np.asarray(y, dtype=float), np.asarray(x, dtype=float)
    n = len(x) - 1
    if n % 2 == 1:
        # 奇數個區間(偶數個點)時退化成加一個梯形法分段, 這個情況在
        # 這個模組裡不會發生(呼叫端都固定用偶數個區間), 保留只是
        # 避免未來有人改了取樣點數忘記這個限制時直接算錯而不自知。
        return _trapz(y, x)
    h = (x[-1] - x[0]) / n
    s = y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2])
    return float(s * h / 3)


def _trapz(y, x):
    """手動實作梯形法數值積分, 不用numpy.trapz/trapezoid——numpy 2.0
    把trapz改名成trapezoid, 為了不管部署環境的numpy版本新舊都能跑,
    這裡自己算(輸入都已經是等距取樣的array, 邏輯很單純)。只有
    _simpson()在區間數是奇數這個理論上不會發生的邊界情況才會退化
    呼叫這個。"""
    y, x = np.asarray(y), np.asarray(x)
    return float(np.sum((y[1:] + y[:-1]) * 0.5 * np.diff(x)))


def _member_inspan_load_resultant(frame, member_id):
    """算這根桿件身上「跨間載重」(distributed_load + member_point_load,
    不包含端點的節點集中力, 那些已經包含在end_forces_local裡了)在
    局部座標系下的合力/合力矩(對i端取矩, 逆時針為正)——用細分數值
    積分, 跟member_internal_forces()同一套局部x/y分解邏輯(重用同一種
    direction='local'/'global_y'/'global'的分解方式), 只是這裡積分
    出「總和」而不是逐點的N/V/M曲線。用途: 畫自由體圖時驗證平衡
    (Fx1+Fx2+這裡的fx_total應該=0, 以此類推), 純粹拿來驗證/繪圖用,
    不是拿去解方程式, 數值積分精度已經足夠。"""
    m = frame.members[member_id]
    ni, nj = frame.nodes[m.node_i], frame.nodes[m.node_j]
    L, angle = member_geometry(ni, nj)
    c_ang, s_ang = np.cos(angle), np.sin(angle)

    fx_total, fy_total, m_about_i = 0.0, 0.0, 0.0

    for dl in frame.distributed_loads:
        if dl.member != member_id:
            continue
        x0 = 0.0 if dl.x_start is None else dl.x_start
        x1 = L if dl.x_end is None else dl.x_end
        w0 = dl.w_start
        w1 = dl.w_start if dl.w_end is None else dl.w_end
        if dl.direction == 'local':
            wx0 = wx1 = 0.0
            wy0, wy1 = w0, w1
        elif dl.direction == 'global_y':
            wx0 = wx1 = s_ang * (-w0)
            wy0 = wy1 = c_ang * (-w0)
        else:  # 'global'
            ang = np.radians(dl.angle_deg)
            ux, uy = np.cos(ang), np.sin(ang)
            R = np.array([[c_ang, s_ang], [-s_ang, c_ang]])
            local0 = R @ (np.array([ux, uy]) * w0)
            local1 = R @ (np.array([ux, uy]) * w1)
            wx0, wy0 = local0
            wx1, wy1 = local1
        # 細分梯形法數值積分(w在x0~x1之間線性變化)
        n = 400
        xs = np.linspace(x0, x1, n + 1)
        wxs = np.linspace(wx0, wx1, n + 1)
        wys = np.linspace(wy0, wy1, n + 1)
        fx_total += _simpson(wxs, xs)
        fy_total += _simpson(wys, xs)
        m_about_i += _simpson(wys * xs, xs)  # 橫向分量對i端的力矩貢獻

    for pl in frame.member_point_loads:
        if pl.member != member_id:
            continue
        if pl.direction == 'global' and pl.F is not None:
            ang = np.radians(pl.angle_deg)
            ux, uy = np.cos(ang), np.sin(ang)
            R = np.array([[c_ang, s_ang], [-s_ang, c_ang]])
            fx, fy = R @ (np.array([ux, uy]) * pl.F)
        else:
            fx, fy = pl.fx, pl.fy
        fx_total += fx
        fy_total += fy
        m_about_i += fy * pl.a + pl.m

    return fx_total, fy_total, m_about_i


def plot_member_fbd(frame, result, member_id, ax=None):
    """畫這根桿件的自由體圖(free body diagram): 桿件本身(不含其他
    桿件/支承), 兩端各畫「其餘結構對這一端的作用力+力矩」(直接來自
    end_forces_local, 局部->全域座標轉換), 桿件身上的跨間載重(均佈
    載重/桿件集中力)也一起畫出來, 並且在圖上標出ΣFx/ΣFy/ΣM(對i端
    取矩)驗證平衡(理論上應該非常接近0, 只有浮點數值誤差)。

    正負號慣例已經用兩個已知答案的案例實際驗證過(懸臂梁尖端載重、
    簡支梁跨間均佈載重), 確認 Fx1+Fx2+跨間load的fx合力=0(以此類推
    Fy, M)這個等式成立。"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    m = frame.members[member_id]
    ni, nj = _member_endpoints(frame, member_id)
    L, angle = member_geometry(ni, nj)
    c_ang, s_ang = np.cos(angle), np.sin(angle)

    ax.plot([ni.x, nj.x], [ni.y, nj.y], color='black', lw=3, zorder=2,
            solid_capstyle='round')
    ax.plot([ni.x, nj.x], [ni.y, nj.y], 'o', color='black', ms=6, zorder=3)
    ax.annotate(f'N{m.node_i}', (ni.x, ni.y),
                xytext=(-10, -14), textcoords='offset points', fontsize=8, ha='right')
    ax.annotate(f'N{m.node_j}', (nj.x, nj.y), xytext=(10, -14), textcoords='offset points', fontsize=8)
    ax.set_title(f'Free Body Diagram: Member {member_id} ({m.member_type}), L = {L:.4g} m', fontsize=11)

    mr = result.member_results[member_id]
    Fx1, Fy1, M1, Fx2, Fy2, M2 = mr.end_forces_local
    # 局部->全域: R(angle) = [[c,-s],[s,c]] 作用在局部座標的力向量上
    Fx1_g = c_ang * Fx1 - s_ang * Fy1
    Fy1_g = s_ang * Fx1 + c_ang * Fy1
    Fx2_g = c_ang * Fx2 - s_ang * Fy2
    Fy2_g = s_ang * Fx2 + c_ang * Fy2

    arrow_len = L * 0.28

    def _draw_end_force(x, y, fx, fy, mm, color):
        # 分量畫法(不合成單一斜向箭頭): Fx/Fy各自畫成獨立的水平/垂直
        # 箭頭, 這樣才是「大小+方向」的完整分量表示, 直接就能跟手算
        # 對照, 不用另外標角度; 跟webapi前端GUI裡反力/外力的「分量」
        # 顯示模式是同一套設計理念。
        if abs(fx) > 1e-6:
            sign = 1 if fx > 0 else -1
            ax.annotate('', xy=(x, y), xytext=(x - sign * arrow_len, y),
                        arrowprops=dict(arrowstyle='-|>', color=color, lw=2, mutation_scale=14), zorder=4)
            ax.annotate(f'Fx={fx:.4g} N', (x - sign * arrow_len, y), fontsize=7.5, color=color,
                        xytext=(4, 6), textcoords='offset points')
        if abs(fy) > 1e-6:
            sign = 1 if fy > 0 else -1
            ax.annotate('', xy=(x, y), xytext=(x, y - sign * arrow_len),
                        arrowprops=dict(arrowstyle='-|>', color=color, lw=2, mutation_scale=14), zorder=4)
            ax.annotate(f'Fy={fy:.4g} N', (x, y - sign * arrow_len), fontsize=7.5, color=color,
                        xytext=(4, -10), textcoords='offset points')
        if abs(mm) > 1e-6:
            _draw_moment_arc(ax, x, y, mm, L * 0.12, color)
            ax.annotate(f'M={mm:.4g} N*m', (x, y), fontsize=8, color=color,
                        xytext=(6, -16), textcoords='offset points')

    _draw_end_force(ni.x, ni.y, Fx1_g, Fy1_g, M1, 'crimson')
    _draw_end_force(nj.x, nj.y, Fx2_g, Fy2_g, M2, 'darkorange')

    # 跨間載重: 重用plot_loads()同一套均佈/集中力畫法會牽涉到整個
    # frame的迴圈跟座標系, 這裡直接針對這一根桿件簡化畫(只畫這根
    # 桿件身上的, 不畫別根桿件或節點的)
    for dl in frame.distributed_loads:
        if dl.member != member_id:
            continue
        x0 = 0.0 if dl.x_start is None else dl.x_start
        x1 = L if dl.x_end is None else dl.x_end
        # 沿桿軸取幾個點畫示意箭頭(不用畫到跟主圖一樣精緻, 這裡只是
        # 標示「這裡有均佈載重」讓自由體圖看得懂受力來源)
        n_arrows = 5
        for t in np.linspace(x0, x1, n_arrows):
            px = ni.x + t * c_ang
            py = ni.y + t * s_ang
            ax.annotate('', xy=(px, py), xytext=(px, py + L * 0.08),
                        arrowprops=dict(arrowstyle='-|>', color='orange', lw=1, mutation_scale=8), zorder=1)
        ax.annotate(f'w={dl.w_start:.4g}' + (f'~{dl.w_end:.4g}' if dl.w_end is not None else '') + ' N/m',
                    (ni.x + (x0 + x1) / 2 * c_ang, ni.y + (x0 + x1) / 2 * s_ang + L * 0.1),
                    fontsize=7, color='darkorange', ha='center')

    for pl in frame.member_point_loads:
        if pl.member != member_id:
            continue
        px = ni.x + pl.a * c_ang
        py = ni.y + pl.a * s_ang
        ax.plot([px], [py], 's', color='purple', ms=5, zorder=3)
        ax.annotate(f'a={pl.a:.3g} m', (px, py), fontsize=7, color='purple', xytext=(4, 8), textcoords='offset points')

    # 平衡驗證: Fx1+Fx2+跨間載重fx合力 應該=0(以此類推Fy, M對i端取矩)
    fx_load, fy_load, m_load = _member_inspan_load_resultant(frame, member_id)
    sum_fx = Fx1 + Fx2 + fx_load
    sum_fy = Fy1 + Fy2 + fy_load
    sum_m = M1 + M2 + Fy2 * L + m_load
    ax.text(0.02, 0.02,
            f'L = {L:.4g} m (member length, node {m.node_i} to node {m.node_j})\n'
            f'Equilibrium check (local coords, units: N, N*m):\n'
            f'ΣFx = {sum_fx:.3g}   ΣFy = {sum_fy:.3g}   ΣM(about node i) = {sum_m:.3g}\n'
            f'(should be ~0, small residual is floating-point noise)',
            transform=ax.transAxes, fontsize=7.5, va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', fc='#f0fdf4', ec='#16a34a', lw=0.8))

    margin = max(L * 0.5, 0.5)
    ax.set_xlim(min(ni.x, nj.x) - margin, max(ni.x, nj.x) + margin)
    ax.set_ylim(min(ni.y, nj.y) - margin, max(ni.y, nj.y) + margin)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(alpha=0.2)
    return ax


def plot_member_own_diagrams(frame, result, member_id, figsize=(10, 8)):
    """單一桿件自己的N/V/M/變形四合一圖(不是整個結構, 只有這一根
    桿件的內力沿桿長分佈+變形形狀), 給查詢/自由體圖功能搭配使用,
    方便針對特定桿件抓出來單獨檢查或設計用。"""
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    x, N, V, M = member_internal_forces(frame, result, member_id, n=101)

    # 變形圖的縮放倍率: 跟plot_deformed()同一套auto-scale邏輯, 只是
    # 這裡只需要考慮這一根桿件自己的偏移量(不用掃整個結構), 讓變形
    # 形狀相對於這根桿件自己的長度L有合理的視覺放大幅度, 不會小到
    # 看起來像一條直線。
    ni, nj = _member_endpoints(frame, member_id)
    L, angle = member_geometry(ni, nj)
    fig.suptitle(f'Member {member_id}: Internal Forces & Deformation (L = {L:.4g} m)', fontsize=11)
    base_x = ni.x + np.linspace(0, L, 41) * np.cos(angle)
    base_y = ni.y + np.linspace(0, L, 41) * np.sin(angle)
    Xr, Yr = member_deformed_shape(frame, result, member_id, scale=1.0, n=41)
    max_offset = max(np.max(np.hypot(Xr - base_x, Yr - base_y)), 1e-9)
    scale = L * 0.15 / max_offset
    X, Y = member_deformed_shape(frame, result, member_id, scale=scale, n=101)

    for ax, vals, label, color, unit in [
        (axes[0, 0], N, 'N (Axial Force)', '#2563eb', 'N'),
        (axes[0, 1], V, 'V (Shear Force)', '#2563eb', 'N'),
        (axes[1, 0], M, 'M (Bending Moment)', '#dc2626', 'N*m'),
    ]:
        ax.plot(x, vals, color=color, lw=1.5)
        ax.fill_between(x, vals, 0, color=color, alpha=0.15)
        ax.axhline(0, color='black', lw=0.5)
        i_max = int(np.argmax(np.abs(vals)))
        ax.annotate(f'{vals[i_max]:.4g} {unit}', (x[i_max], vals[i_max]), fontsize=8, color=color)
        ax.set_title(f'Member {member_id}: {label}', fontsize=10)
        ax.set_xlabel('x along member (m)')
        ax.set_ylabel(f'{label.split()[0]} ({unit})')
        ax.grid(alpha=0.2)

    ax4 = axes[1, 1]
    ax4.plot(base_x, base_y, '--', color='gray', lw=1, alpha=0.6)
    ax4.plot(X, Y, color='#16a34a', lw=1.5)
    ax4.set_title(f'Member {member_id}: Deformed Shape (scale x{scale:.2g})', fontsize=10)
    ax4.set_xlabel('x (m)')
    ax4.set_ylabel('y (m)')
    ax4.set_aspect('equal')
    ax4.grid(alpha=0.2)

    fig.tight_layout()
    return fig


def plot_all(frame, result, figsize=(14, 9)):
    """六合一總覽圖"""
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    plot_structure(frame, ax=axes[0, 0])
    plot_loads(frame, ax=axes[0, 1])
    plot_deformed(frame, result, ax=axes[0, 2])
    plot_diagram(frame, result, 'N', ax=axes[1, 0])
    plot_diagram(frame, result, 'V', ax=axes[1, 1])
    plot_diagram(frame, result, 'M', ax=axes[1, 2])
    fig.tight_layout()
    return fig
