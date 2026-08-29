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
                parts.append(f'Fx={pl.fx:.3g}')
            if abs(pl.fy) > 1e-9:
                parts.append(f'Fy={pl.fy:.3g}')
            ax.annotate(', '.join(parts), (tail_x, tail_y), color='red', fontsize=8,
                        ha='center', xytext=(0, -10 if dy >= 0 else 10), textcoords='offset points')
            extra_pts_x += [tail_x]
            extra_pts_y += [tail_y]
        if abs(pl.m) > 1e-9:
            ax.annotate(f'M={pl.m:.3g}', (n.x, n.y), color='purple', fontsize=8,
                        xytext=(8, -12), textcoords='offset points')

    for dl in frame.distributed_loads:
        ni, nj = _member_endpoints(frame, dl.member)
        L, angle = member_geometry(ni, nj)
        c, s = np.cos(angle), np.sin(angle)
        x_start = 0.0 if dl.x_start is None else dl.x_start
        x_end = L if dl.x_end is None else dl.x_end
        n_arrows = 6
        max_w = max(abs(dl.w_start), abs(dl.w_end), 1e-9)
        arrow_len = scale * 0.08
        for k in range(n_arrows + 1):
            t = k / n_arrows
            wx = dl.w_start + (dl.w_end - dl.w_start) * t
            xt = x_start + (x_end - x_start) * t   # 只在[x_start,x_end]範圍內畫箭頭
            px, py = ni.x + xt * c, ni.y + xt * s
            # 箭頭方向沿局部+y (perp to 桿件), 長度依載重大小比例
            perp_x, perp_y = -s, c
            length = arrow_len * (wx / max_w) if max_w > 0 else 0
            tail_x, tail_y = px - perp_x * length, py - perp_y * length
            ax.annotate('', xy=(px, py), xytext=(tail_x, tail_y),
                        arrowprops=dict(arrowstyle='->', color='orange', lw=1))
            extra_pts_x += [tail_x]
            extra_pts_y += [tail_y]
        if abs(dl.w_start - dl.w_end) < 1e-9:
            w_label = f'w={dl.w_start:.3g}'
        else:
            w_label = f'w={dl.w_start:.3g}~{dl.w_end:.3g}'
        mx = ni.x + (x_start + x_end) / 2 * c
        my = ni.y + (x_start + x_end) / 2 * s
        ax.annotate(w_label, (mx, my), color='orange', fontsize=8,
                    xytext=(0, 10), textcoords='offset points', ha='center')

    for pl in frame.member_point_loads:
        ni, nj = _member_endpoints(frame, pl.member)
        L, angle = member_geometry(ni, nj)
        c, s = np.cos(angle), np.sin(angle)
        px, py = ni.x + pl.a * c, ni.y + pl.a * s   # 桿件內部位置(局部座標a換算成全域座標)
        mag = np.hypot(pl.fx, pl.fy)
        if mag > 1e-9:
            # fx,fy是局部座標分量, 換算成全域方向畫箭頭
            gdx = pl.fx * c - pl.fy * s
            gdy = pl.fx * s + pl.fy * c
            Larrow = scale * 0.15
            dx, dy = gdx / mag * Larrow, gdy / mag * Larrow
            tail_x, tail_y = px - dx, py - dy
            ax.annotate('', xy=(px, py), xytext=(tail_x, tail_y),
                        arrowprops=dict(arrowstyle='->', color='crimson', lw=2))
            parts = []
            if abs(pl.fx) > 1e-9:
                parts.append(f'fx={pl.fx:.3g}')
            if abs(pl.fy) > 1e-9:
                parts.append(f'fy={pl.fy:.3g}')
            ax.annotate(', '.join(parts), (tail_x, tail_y), color='crimson', fontsize=8,
                        ha='center', xytext=(0, -10 if dy >= 0 else 10), textcoords='offset points')
            extra_pts_x += [tail_x]
            extra_pts_y += [tail_y]
        if abs(pl.m) > 1e-9:
            ax.annotate(f'm={pl.m:.3g}', (px, py), color='darkmagenta', fontsize=8,
                        xytext=(8, -12), textcoords='offset points')
            ax.plot(px, py, 'D', color='darkmagenta', ms=5, zorder=2)

    # 直接手動計算範圍再 set_xlim/ylim, 比依賴 relim/autoscale 對 annotate 更可靠
    if extra_pts_x:
        xs = [n.x for n in frame.nodes.values()] + extra_pts_x
        ys = [n.y for n in frame.nodes.values()] + extra_pts_y
        pad_x = (max(xs) - min(xs)) * 0.1 + 1e-6
        pad_y = (max(ys) - min(ys)) * 0.1 + 1e-6
        ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
        ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)

    ax.set_title('Loads')
    return ax


def plot_diagram(frame, result, kind, ax=None, scale=None):
    """③④⑤ 軸力圖/剪力圖/彎矩圖. kind: 'N', 'V', 或 'M' """
    if ax is None:
        fig, ax = plt.subplots()
    plot_structure(frame, ax=ax, show_node_ids=False, show_member_ids=False, show_dimensions=False)

    idx = {'N': 1, 'V': 2, 'M': 3}[kind]
    label = {'N': 'Axial Force (N)', 'V': 'Shear Force (V)', 'M': 'Bending Moment (M)'}[kind]

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
            ax.annotate(f'Δ={mag:.4g}', (n.x, n.y), fontsize=7, color='darkblue',
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
            ax.annotate(f'Δmax={offset_real[i_max]:.4g}', (X1[i_max], Y1[i_max]), fontsize=7, color='darkblue',
                        xytext=(5, -12), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='darkblue', lw=0.5, alpha=0.85))

    ax.set_title(f'Deformation (scale x{scale:.1f})')
    return ax


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
