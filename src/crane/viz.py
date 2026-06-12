"""phase portrait と stick-figure アニメーション（slope frame 描画）。"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np

from crane.stride import StrideResult


def link_points(theta: float, phi: float, foot_x: float) -> tuple[np.ndarray, np.ndarray]:
    """slope frame での hip / swing 足の座標 (脚長 l=1)。

    stance 足を (foot_x, 0) に置く。脚角 ψ の足先は hip + (sin ψ, −cos ψ)。
    stance: ψ=θ、swing: ψ=θ−φ。
    """
    hip = np.array([foot_x - np.sin(theta), np.cos(theta)])
    psi = theta - phi
    swing = hip + np.array([np.sin(psi), -np.cos(psi)])
    return hip, swing


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


def animate_walk(strides: list[StrideResult], gamma: float, out: Path, fps: int = 30) -> None:
    """stick-figure アニメ。slope frame で描き、全体を −γ 回転して坂を見せる。"""
    rot = np.array([[np.cos(-gamma), -np.sin(-gamma)], [np.sin(-gamma), np.cos(-gamma)]])

    # 各 stride の stance 足アンカーを heel-strike 位置で更新しながらフレーム列を作る
    frames: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    foot_x = 0.0
    dt_frame = 1.0 / fps
    for s in strides:
        t_resampled = np.arange(0.0, s.t[-1], dt_frame)
        for ti in t_resampled:
            idx = int(np.searchsorted(s.t, ti))
            idx = min(idx, s.t.size - 1)
            theta, phi = s.x[0, idx], s.x[1, idx]
            hip, swing = link_points(theta, phi, foot_x)
            foot = np.array([foot_x, 0.0])
            frames.append((rot @ foot, rot @ hip, rot @ swing))
        # heel-strike: swing 足の着地点が次の stance アンカー
        _, swing_end = link_points(s.x_strike[0], s.x_strike[1], foot_x)
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
