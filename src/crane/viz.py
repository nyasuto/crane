"""phase portrait と stick-figure アニメーション（slope frame 描画）。"""

from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
from matplotlib.colors import ListedColormap  # noqa: E402

from crane.stride import StrideResult


def link_points_abs(theta_st: float, psi_sw: float, foot_x: float) -> tuple[np.ndarray, np.ndarray]:
    """slope frame での hip / swing 足の座標（絶対角版、脚長 l=1）。

    stance 足を (foot_x, 0) に置く。脚角 ψ の足先は hip + (sin ψ, −cos ψ)。
    theta_st: stance 脚の絶対角、psi_sw: swing 脚の絶対角。
    """
    hip = np.array([foot_x - np.sin(theta_st), np.cos(theta_st)])
    swing = hip + np.array([np.sin(psi_sw), -np.cos(psi_sw)])
    return hip, swing


def link_points(theta: float, phi: float, foot_x: float) -> tuple[np.ndarray, np.ndarray]:
    """slope frame での hip / swing 足の座標 (脚長 l=1)。後方互換ラッパー。

    stance 足を (foot_x, 0) に置く。脚角 ψ の足先は hip + (sin ψ, −cos ψ)。
    stance: ψ=θ、swing: ψ=θ−φ。
    """
    return link_points_abs(theta, theta - phi, foot_x)


def plot_phase_portrait(strides: list[StrideResult], out: Path) -> None:
    """θ–θ̇ 平面の軌道。リミットサイクルなら閉軌道が重なる。"""
    fig, ax = plt.subplots(figsize=(6, 5))
    for s in strides:
        ax.plot(s.x[0], s.x[2], lw=0.8)
    ax.set_xlabel("theta [rad]")
    ax.set_ylabel("theta_dot [rad/s]")
    ax.set_title(f"Phase portrait ({len(strides)} strides)")
    fig.savefig(out, dpi=150)
    plt.close(fig)


def animate_walk(
    strides: list[StrideResult],
    gamma: float,
    out: Path,
    fps: int = 30,
    angles_of: Callable[[float, float], tuple[float, float]] | None = None,
) -> None:
    """stick-figure アニメ。slope frame で描き、全体を −γ 回転して坂を見せる。

    angles_of: (q0, q1) -> (stance 絶対角, swing 絶対角)。
    None なら simplest 規約 lambda q0, q1: (q0, q0 - q1)。
    """
    if angles_of is None:
        angles_of = lambda q0, q1: (q0, q0 - q1)  # noqa: E731

    rot = np.array([[np.cos(-gamma), -np.sin(-gamma)], [np.sin(-gamma), np.cos(-gamma)]])

    # 各 stride の stance 足アンカーを heel-strike 位置で更新しながらフレーム列を作る
    frames: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    foot_x = 0.0
    dt_frame = 1.0 / fps
    for s in strides:
        # strike 瞬間を最終フレームとして含める（arange は末端を除くため append）
        t_resampled = np.append(np.arange(0.0, s.t[-1], dt_frame), s.t[-1])
        q0_i = np.interp(t_resampled, s.t, s.x[0])
        q1_i = np.interp(t_resampled, s.t, s.x[1])
        for q0, q1 in zip(q0_i, q1_i):
            theta_st, psi_sw = angles_of(float(q0), float(q1))
            hip, swing = link_points_abs(theta_st, psi_sw, foot_x)
            foot = np.array([foot_x, 0.0])
            frames.append((rot @ foot, rot @ hip, rot @ swing))
        # heel-strike: x_strike（衝突前状態）の swing 足位置を次 stride の stance アンカーに使う。
        # 衝突後の leg-swap より前の着地ジオメトリが正しい foot_x を与えるため x_strike を参照する。
        theta_st_s, psi_sw_s = angles_of(float(s.x_strike[0]), float(s.x_strike[1]))
        _, swing_end = link_points_abs(theta_st_s, psi_sw_s, foot_x)
        foot_x = float(swing_end[0])

    fig, ax = plt.subplots(figsize=(8, 4))
    (stance_line,) = ax.plot([], [], "o-", lw=2, color="tab:blue")
    (swing_line,) = ax.plot([], [], "o-", lw=2, color="tab:orange")
    span = max(foot_x + 2.0, 4.0)
    ground = rot @ np.array([[-1.0, span], [0.0, 0.0]])
    ax.plot(ground[0], ground[1], "k-", lw=1)
    ax.set_xlim(-1.0, span)
    ax.set_ylim(-span * np.sin(gamma) - 0.5, 1.5)
    ax.set_aspect("equal")

    def update(i: int):
        foot, hip, swing = frames[i]
        stance_line.set_data([foot[0], hip[0]], [foot[1], hip[1]])
        swing_line.set_data([hip[0], swing[0]], [hip[1], swing[1]])
        return stance_line, swing_line

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)
    try:
        anim.save(out, writer=animation.FFMpegWriter(fps=fps))
    except (FileNotFoundError, RuntimeError):
        # ffmpeg が無い環境では GIF にフォールバック
        out = out.with_suffix(".gif")
        anim.save(out, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def rocker_joints(x: np.ndarray, L: float, R: float, contact_x: float):
    """rocker compass の関節座標（slope frame）。stance 接触点を (contact_x, 0) に。

    返り値: (stance_contact, stance_C, hip, swing_C, swing_lowest)
    """
    th_st, th_sw = x[0], x[1]

    def down(theta):
        return np.array([np.sin(theta), -np.cos(theta)])

    C_st = np.array([contact_x, R])
    hip = C_st - (L - R) * down(th_st)
    C_sw = hip + (L - R) * down(th_sw)
    swing_lowest = C_sw + np.array([0.0, -R])
    return np.array([contact_x, 0.0]), C_st, hip, C_sw, swing_lowest


def _arc_foot(center: np.ndarray, theta: float, R: float, rot: np.ndarray, n: int = 24):
    """曲率中心 center・半径 R の円弧足を polyline でサンプル（脚軸 down(θ) を中心に半円弧）。

    脚軸方向（θ=0 で接地点が真下）を中心に ±90° の半円弧を描いて足裏の丸みを見せる。
    rot で −γ 回転して返す。
    """
    base = np.arctan2(-np.cos(theta), np.sin(theta))  # down(θ) 方向の角度（+x から）
    angs = base + np.linspace(-np.pi / 2, np.pi / 2, n)
    pts = center[:, None] + R * np.array([np.cos(angs), np.sin(angs)])
    return rot @ pts


def animate_rocker(
    strides: list[StrideResult],
    L: float,
    R: float,
    gamma: float,
    out: Path,
    fps: int = 30,
) -> None:
    """rocker-foot compass の stick-figure アニメ（slope frame、−γ 回転）。

    stance 脚（hip→曲率中心 C_st、青）と swing 脚（hip→C_sw、橙）を描き、各脚先に
    半径 R の円弧足を polyline で描いて丸足の転がりを見せる。

    点足 compass/kneed と違い、円弧足は step 内で接地点が転がって移動する
    （モデルでは P_st_x = −R·θ_st）。各フレームで rolling-without-slip
    Δcontact = −R·Δθ_st を累積して stance 接触 x を更新し、丸足が脚の下を
    前転する様子を描く（R→0 で転がり項が消え点足と一致）。heel-strike 時には
    着地ジオメトリ（x_strike, leg-swap 前）の swing 円弧足の接地 x（曲率中心
    C_sw の真下＝C_sw_x）を次 step の開始 contact_x にして滑らかに連結する。
    """
    rot = np.array([[np.cos(-gamma), -np.sin(-gamma)], [np.sin(-gamma), np.cos(-gamma)]])

    frames: list[tuple[np.ndarray, ...]] = []
    contact_x = 0.0
    dt_frame = 1.0 / fps
    for s in strides:
        t_resampled = np.append(np.arange(0.0, s.t[-1], dt_frame), s.t[-1])
        q0_i = np.interp(t_resampled, s.t, s.x[0])
        q1_i = np.interp(t_resampled, s.t, s.x[1])
        prev_q0 = None
        for q0, q1 in zip(q0_i, q1_i):
            # step 内の転がり: Δcontact = −R·Δθ_st（rolling-without-slip）を累積
            if prev_q0 is not None:
                contact_x += -R * (float(q0) - prev_q0)
            prev_q0 = float(q0)
            _, C_st, hip, C_sw, _ = rocker_joints(np.array([q0, q1, 0, 0]), L, R, contact_x)
            st_foot = _arc_foot(C_st, float(q0), R, rot)
            sw_foot = _arc_foot(C_sw, float(q1), R, rot)
            frames.append((rot @ C_st, rot @ hip, rot @ C_sw, st_foot, sw_foot))
        # heel-strike 着地ジオメトリ（leg-swap 前）の swing 円弧足の接地 x を次 step の
        # 開始 contact_x に。最後のフレームまで転がした contact_x で C_sw_x を評価して連結。
        _, _, _, C_sw_s, _ = rocker_joints(s.x_strike, L, R, contact_x)
        contact_x = float(C_sw_s[0])

    fig, ax = plt.subplots(figsize=(8, 4))
    (stance_line,) = ax.plot([], [], "o-", lw=2, color="tab:blue")
    (swing_line,) = ax.plot([], [], "o-", lw=2, color="tab:orange")
    (stance_foot,) = ax.plot([], [], "-", lw=2.5, color="tab:blue")
    (swing_foot,) = ax.plot([], [], "-", lw=2.5, color="tab:orange")
    span = max(contact_x + 2.0, 4.0)
    ground = rot @ np.array([[-1.0, span], [0.0, 0.0]])
    ax.plot(ground[0], ground[1], "k-", lw=1)
    ax.set_xlim(-1.0, span)
    ax.set_ylim(-span * np.sin(gamma) - 0.5, 1.5)
    ax.set_aspect("equal")

    def update(i: int):
        C_st, hip, C_sw, st_foot, sw_foot = frames[i]
        stance_line.set_data([hip[0], C_st[0]], [hip[1], C_st[1]])
        swing_line.set_data([hip[0], C_sw[0]], [hip[1], C_sw[1]])
        stance_foot.set_data(st_foot[0], st_foot[1])
        swing_foot.set_data(sw_foot[0], sw_foot[1])
        return stance_line, swing_line, stance_foot, swing_foot

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)
    try:
        anim.save(out, writer=animation.FFMpegWriter(fps=fps))
    except (FileNotFoundError, RuntimeError):
        # ffmpeg が無い環境では GIF にフォールバック
        out = out.with_suffix(".gif")
        anim.save(out, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def kneed_joints(x: np.ndarray, l_t: float, l_s: float, foot_x: float):
    """kneed の関節座標列（slope frame）。stance 足を (foot_x, 0) に置く。

    状態 x = [θ_st, θ_th, θ_sh, ...]（絶対角）。stance 脚は直線（膝ロック）。
    返り値: (stance_foot, stance_knee, hip, swing_knee, swing_foot)
    """
    th_st, th_th, th_sh = x[0], x[1], x[2]
    length = l_t + l_s

    def down(theta):
        return np.array([np.sin(theta), -np.cos(theta)])

    foot = np.array([foot_x, 0.0])
    hip = foot - length * down(th_st)
    knee_st = hip + l_t * down(th_st)
    knee_sw = hip + l_t * down(th_th)
    foot_sw = knee_sw + l_s * down(th_sh)
    return foot, knee_st, hip, knee_sw, foot_sw


def animate_kneed(
    strides: list[StrideResult],
    l_t: float,
    l_s: float,
    gamma: float,
    out: Path,
    fps: int = 30,
) -> None:
    """kneed walker の 4 セグメント stick-figure アニメ（slope frame、−γ 回転）。

    stance 脚（foot–knee–hip、青）と swing 脚（hip–knee–foot、橙）の 2 本の
    polyline。膝が屈曲する unlocked 相を含め、3 角度をそれぞれ補間して描く。
    stance 足アンカーは heel-strike 時（x_strike）の swing 足位置で前進させる。
    """
    rot = np.array([[np.cos(-gamma), -np.sin(-gamma)], [np.sin(-gamma), np.cos(-gamma)]])

    frames: list[tuple[np.ndarray, ...]] = []
    foot_x = 0.0
    dt_frame = 1.0 / fps
    for s in strides:
        t_resampled = np.append(np.arange(0.0, s.t[-1], dt_frame), s.t[-1])
        q0_i = np.interp(t_resampled, s.t, s.x[0])
        q1_i = np.interp(t_resampled, s.t, s.x[1])
        q2_i = np.interp(t_resampled, s.t, s.x[2])
        for q0, q1, q2 in zip(q0_i, q1_i, q2_i):
            pts = kneed_joints(np.array([q0, q1, q2]), l_t, l_s, foot_x)
            frames.append(tuple(rot @ p for p in pts))
        # heel-strike 着地ジオメトリ（leg-swap 前）の swing 足 x を次アンカーに
        pts_s = kneed_joints(s.x_strike, l_t, l_s, foot_x)
        foot_x = float(pts_s[4][0])

    fig, ax = plt.subplots(figsize=(8, 4))
    (stance_line,) = ax.plot([], [], "o-", lw=2, color="tab:blue")
    (swing_line,) = ax.plot([], [], "o-", lw=2, color="tab:orange")
    span = max(foot_x + 2.0, 4.0)
    ground = rot @ np.array([[-1.0, span], [0.0, 0.0]])
    ax.plot(ground[0], ground[1], "k-", lw=1)
    ax.set_xlim(-1.0, span)
    ax.set_ylim(-span * np.sin(gamma) - 0.5, 1.5)
    ax.set_aspect("equal")

    def update(i: int):
        foot, knee_st, hip, knee_sw, foot_sw = frames[i]
        stance_line.set_data([foot[0], knee_st[0], hip[0]], [foot[1], knee_st[1], hip[1]])
        swing_line.set_data([hip[0], knee_sw[0], foot_sw[0]], [hip[1], knee_sw[1], foot_sw[1]])
        return stance_line, swing_line

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)
    try:
        anim.save(out, writer=animation.FFMpegWriter(fps=fps))
    except (FileNotFoundError, RuntimeError):
        out = out.with_suffix(".gif")
        anim.save(out, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def rocker_kneed_joints(x: np.ndarray, l_t: float, l_s: float, R: float, contact_x: float):
    """rocker kneed の関節座標（slope frame）。stance 接触点を (contact_x, 0) に。

    返り値: (stance_contact, stance_C, knee_st, hip, knee_sw, swing_C, swing_lowest)
    """
    th_st, th_th, th_sh = x[0], x[1], x[2]
    length = l_t + l_s

    def down(theta):
        return np.array([np.sin(theta), -np.cos(theta)])

    C_st = np.array([contact_x, R])
    hip = C_st - (length - R) * down(th_st)
    knee_st = hip + l_t * down(th_st)
    knee_sw = hip + l_t * down(th_th)
    C_sw = knee_sw + (l_s - R) * down(th_sh)
    swing_lowest = C_sw + np.array([0.0, -R])
    return np.array([contact_x, 0.0]), C_st, knee_st, hip, knee_sw, C_sw, swing_lowest


def animate_rocker_kneed(
    strides: list[StrideResult],
    l_t: float,
    l_s: float,
    R: float,
    gamma: float,
    out: Path,
    fps: int = 30,
) -> None:
    """rocker-foot kneed walker の 4 セグメント + 円弧足 stick-figure アニメ。

    animate_kneed（4 セグメント描画）と animate_rocker（円弧足の step 内転がり）を
    合成する。stance 脚は polyline contact→C_st→knee_st→hip（青）、swing 脚は
    hip→knee_sw→C_sw（橙）として描き、C_st / C_sw それぞれに半径 R の円弧足
    （`_arc_foot`）を添える。膝は swing 大腿角 θ_th と脛角 θ_sh の差で屈曲する。

    stance 接触点 contact_x は step 内で rolling-without-slip
    Δcontact = −R·Δθ_st を累積して前転し、heel-strike 時には着地ジオメトリ
    （x_strike, leg-swap 前）の swing 円弧足の接地 x（曲率中心 C_sw の真下）を
    次 step の開始 contact_x にして連結する。全体を −γ 回転して坂を見せる。
    """
    rot = np.array([[np.cos(-gamma), -np.sin(-gamma)], [np.sin(-gamma), np.cos(-gamma)]])

    frames: list[tuple[np.ndarray, ...]] = []
    contact_x = 0.0
    dt_frame = 1.0 / fps
    for s in strides:
        t_resampled = np.append(np.arange(0.0, s.t[-1], dt_frame), s.t[-1])
        q0_i = np.interp(t_resampled, s.t, s.x[0])
        q1_i = np.interp(t_resampled, s.t, s.x[1])
        q2_i = np.interp(t_resampled, s.t, s.x[2])
        prev_q0 = None
        for q0, q1, q2 in zip(q0_i, q1_i, q2_i):
            # step 内の転がり: Δcontact = −R·Δθ_st（rolling-without-slip）を累積
            if prev_q0 is not None:
                contact_x += -R * (float(q0) - prev_q0)
            prev_q0 = float(q0)
            _, C_st, knee_st, hip, knee_sw, C_sw, _ = rocker_kneed_joints(
                np.array([q0, q1, q2]), l_t, l_s, R, contact_x
            )
            st_foot = _arc_foot(C_st, float(q0), R, rot)
            sw_foot = _arc_foot(C_sw, float(q2), R, rot)
            frames.append(
                (rot @ C_st, rot @ knee_st, rot @ hip, rot @ knee_sw, rot @ C_sw, st_foot, sw_foot)
            )
        # heel-strike 着地ジオメトリ（leg-swap 前）の swing 円弧足の接地 x を次 step の
        # 開始 contact_x に。最後のフレームまで転がした contact_x で C_sw_x を評価して連結。
        _, _, _, _, _, C_sw_s, _ = rocker_kneed_joints(s.x_strike, l_t, l_s, R, contact_x)
        contact_x = float(C_sw_s[0])

    fig, ax = plt.subplots(figsize=(8, 4))
    # 関節マーカーは小さめに（大きな点が円弧足と膝屈曲を隠さないように）。膝の屈曲が
    # 見えるよう knee を強調し、円弧足は太い線で丸みを見せる。
    (stance_line,) = ax.plot([], [], "o-", lw=2, ms=5, color="tab:blue")
    (swing_line,) = ax.plot([], [], "o-", lw=2, ms=5, color="tab:orange")
    (stance_foot,) = ax.plot([], [], "-", lw=2.5, color="tab:blue")
    (swing_foot,) = ax.plot([], [], "-", lw=2.5, color="tab:orange")
    span = max(contact_x + 2.0, 4.0)
    ground = rot @ np.array([[-1.0, span], [0.0, 0.0]])
    ax.plot(ground[0], ground[1], "k-", lw=1)
    # 歩行帯に合わせて縦横を確保しつつ脚長スケールで walker が見えるようにする
    ax.set_xlim(-1.0, span)
    ax.set_ylim(-span * np.sin(gamma) - 0.5, 1.5)
    ax.set_aspect("equal")

    def update(i: int):
        C_st, knee_st, hip, knee_sw, C_sw, st_foot, sw_foot = frames[i]
        stance_line.set_data([C_st[0], knee_st[0], hip[0]], [C_st[1], knee_st[1], hip[1]])
        swing_line.set_data([hip[0], knee_sw[0], C_sw[0]], [hip[1], knee_sw[1], C_sw[1]])
        stance_foot.set_data(st_foot[0], st_foot[1])
        swing_foot.set_data(sw_foot[0], sw_foot[1])
        return stance_line, swing_line, stance_foot, swing_foot

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)
    try:
        anim.save(out, writer=animation.FFMpegWriter(fps=fps))
    except (FileNotFoundError, RuntimeError):
        # ffmpeg が無い環境では GIF にフォールバック
        out = out.with_suffix(".gif")
        anim.save(out, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def plot_basin(result, out: Path) -> None:
    """basin 分類グリッドを描画。CONVERGED=緑, FELL=赤, UNDECIDED=灰。"""
    cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])  # 0,1,2
    fig, ax = plt.subplots(figsize=(6, 5))
    extent = [
        result.ax0_vals[0],
        result.ax0_vals[-1],
        result.ax1_vals[0],
        result.ax1_vals[-1],
    ]
    ax.imshow(
        result.grid,
        origin="lower",
        extent=extent,
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=2,
        interpolation="nearest",
    )
    a0, a1 = result.axes
    ax.plot(result.fixed_point[a0], result.fixed_point[a1], "k*", markersize=12)
    ax.set_xlabel(f"section[{a0}]")
    ax.set_ylabel(f"section[{a1}]")
    ax.set_title(f"{result.model_name}  basin_fraction={result.basin_fraction:.3f}")
    fig.savefig(out, dpi=150)
    plt.close(fig)
