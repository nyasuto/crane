# scripts/wheel_limit.py
"""Phase 4a.1: rocker_compass の R を掃引し「車輪への漸近」を検証。

固定 γ=0.030 で R∈(0,L) を continuation 追跡。各 R で δ(R)・max|λ|(R) を測り、
粗い R で basin_fraction(R) を測る。曲線・JSON・basin montage を data/runs/ に出力。
Usage: uv run python scripts/wheel_limit.py [--dr 0.05] [--basin-res 50] [--workers N]
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from crane import references_mcgeer as ref  # noqa: E402
from crane.basin import basin_slice  # noqa: E402
from crane.efficiency import mechanical_cot, relative_loss, step_collision_loss  # noqa: E402
from crane.models.rocker_compass import (  # noqa: E402
    RockerCompassParams,
    kinetic_energy,
    make_rocker_compass,
)
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402
from crane.stride import StrideError, stride  # noqa: E402


def _params(R):
    return RockerCompassParams(
        m=ref.M_LEG,
        m_h=ref.M_HIP,
        c=ref.C_HIP_TO_COM,
        rho=ref.RHO_GYR,
        L=ref.L_LEG,
        R=R,
        gamma=ref.GAMMA_GAIT,
        g=ref.G,
    )


def _hip_x(theta_st, R, L):
    """stance 接触フレームでの hip 水平位置: -R*theta - (L-R)*sin(theta)。"""
    return -R * theta_st - (L - R) * np.sin(theta_st)


def _measure(R, guess):
    """R で limit cycle を求め (fp, delta, cot, lam_max) を返す。失敗時 None。"""
    p = _params(R)
    model = make_rocker_compass(p)
    fp = find_limit_cycle(model, np.array(guess))
    if not fp.converged:
        return None
    lam_max = float(np.max(np.abs(fp.eigenvalues)))
    try:
        result = stride(model, model.lift(fp.y))
    except StrideError:
        return None
    loss, ke_pre, ke_post = step_collision_loss(
        result.x_strike, result.x_end, lambda x: kinetic_energy(x, p)
    )
    delta = relative_loss(ke_pre, ke_post)
    theta_start = float(result.x[0, 0])
    theta_strike = float(result.x_strike[0])
    step_len = abs(_hip_x(theta_strike, R, p.L) - _hip_x(theta_start, R, p.L))
    m_total = 2.0 * p.m + p.m_h  # 総質量（2脚 + hip）。重力入力に対応
    cot = mechanical_cot(loss, m=m_total, g=p.g, step_length=step_len) if step_len > 1e-9 else None
    return fp, delta, cot, lam_max


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dr", type=float, default=0.05)
    parser.add_argument("--basin-res", type=int, default=50)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    records = {}
    base = _measure(0.3, ref.SECTION_GUESS)
    if base is None:
        raise SystemExit("ERROR: R=0.3 limit cycle not found")
    records[0.3] = base

    def sweep(direction):
        guess = records[0.3][0].y
        R = 0.3
        while True:
            R = round(R + direction * args.dr, 4)
            if R <= 0.0 or R >= ref.L_LEG:
                break
            m = _measure(R, guess)
            if m is None or m[3] >= 1.0:
                records[R] = None
                break
            records[R] = m
            guess = m[0].y

    sweep(+1)
    sweep(-1)

    Rs = sorted(r for r, v in records.items() if v is not None)
    delta_curve = [records[r][1] for r in Rs]
    cot_curve = [records[r][2] for r in Rs]
    lam_curve = [records[r][3] for r in Rs]
    print("R sweep:", [(r, round(records[r][1], 4), round(records[r][3], 3)) for r in Rs])

    basin_Rs = [r for r in (0.05, 0.2, 0.4, 0.6, 0.8) if r in records and records[r] is not None]
    if Rs and Rs[-1] not in basin_Rs:
        basin_Rs.append(Rs[-1])
    basin_fracs = {}
    basin_grids = {}
    for r in basin_Rs:
        fp = records[r][0]
        res = basin_slice(
            make_rocker_compass,
            _params(r),
            fp.y,
            axes=(0, 1),
            half_widths=(0.16, 0.65),
            resolution=args.basin_res,
            model_name=f"rocker_compass_R{r}",
            n_workers=args.workers,
        )
        basin_fracs[r] = res.basin_fraction
        basin_grids[r] = res
        print(f"basin R={r}: fraction={res.basin_fraction:.3f}")

    run_dir = new_run_dir("wheel_limit")
    fig, axs = plt.subplots(3, 1, figsize=(7, 10), sharex=True)
    axs[0].plot(Rs, delta_curve, "o-")
    axs[0].set_ylabel("delta (relative loss)")
    axs[0].axhline(0.0, color="gray", lw=0.5)
    axs[1].plot(Rs, lam_curve, "o-")
    axs[1].axhline(1.0, color="red", lw=0.5)
    axs[1].set_ylabel("max|lambda|")
    if basin_fracs:
        axs[2].plot(list(basin_fracs), list(basin_fracs.values()), "s-")
    axs[2].set_ylabel("basin_fraction")
    axs[2].set_xlabel("R (foot radius); R->L=1 is the wheel limit")
    fig.suptitle("Approaching the wheel: rocker_compass R-continuation (gamma=0.030)")
    fig.savefig(run_dir / "R_sweep.png", dpi=150)
    plt.close(fig)

    if basin_grids:
        cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])
        order = sorted(basin_grids)
        fig2, axs2 = plt.subplots(1, len(order), figsize=(4 * len(order), 4))
        if len(order) == 1:
            axs2 = [axs2]
        for ax, r in zip(axs2, order):
            g = basin_grids[r]
            ax.imshow(
                g.grid,
                origin="lower",
                extent=[g.ax0_vals[0], g.ax0_vals[-1], g.ax1_vals[0], g.ax1_vals[-1]],
                aspect="auto",
                cmap=cmap,
                vmin=0,
                vmax=2,
                interpolation="nearest",
            )
            ax.plot(g.fixed_point[0], g.fixed_point[1], "k*", markersize=10)
            ax.set_title(f"R={r}\nfrac={g.basin_fraction:.3f}")
        fig2.tight_layout()
        fig2.savefig(run_dir / "basin_R_montage.png", dpi=150)
        plt.close(fig2)

    data = {
        "gamma": ref.GAMMA_GAIT,
        "sin_gamma": float(np.sin(ref.GAMMA_GAIT)),
        "R": Rs,
        "delta": delta_curve,
        "cot": cot_curve,
        "max_lambda": lam_curve,
        "fixed_points": [records[r][0].y.tolist() for r in Rs],
        "basin_fraction": basin_fracs,
        "walking_vanishes_at": (min((r for r, v in records.items() if v is None), default=None)),
    }
    (run_dir / "R_sweep.json").write_text(json.dumps(data, indent=2, default=str))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
