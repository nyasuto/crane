# scripts/active_basin_sweep.py
"""Phase 5c: γ=0.030 で R×push_off 格子の能動 basin を測り受動と比較。

受動不動点は R=0.3 から R-continuation で各 R に供給し、能動不動点はそれを seed に直接収束。
basin_fraction(R,P) の曲線・montage・JSON（順位・変動係数つき）を data/runs/ に出力。
Usage: uv run python scripts/active_basin_sweep.py [--basin-res 50] [--workers N]
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from crane.basin import basin_slice  # noqa: E402
from crane.models.powered_rocker_compass import (  # noqa: E402
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass  # noqa: E402
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402

NOM = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, gamma=0.030, g=9.81)
R_VALUES = [0.05, 0.20, 0.40, 0.60]
P_VALUES = [0.0, 0.04, 0.08]


def passive_fp(R, n_steps=12):
    """R=0.3 から R-continuation で R の受動不動点を返す。失敗時 None。"""
    guess = np.array([0.30844, -1.26256, -0.87914])
    for r in np.linspace(0.3, R, n_steps):
        fp = find_limit_cycle(make_rocker_compass(RockerCompassParams(**NOM, R=float(r))), guess)
        if not fp.converged:
            return None
        guess = fp.y
    return guess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basin-res", type=int, default=50)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    fracs = {}
    grids = {}
    for R in R_VALUES:
        seed = passive_fp(R)
        if seed is None:
            print(f"R={R}: passive fp not found, skip")
            continue
        for P in P_VALUES:
            params = PoweredRockerCompassParams(**NOM, R=R, push_off=P)
            fp = find_limit_cycle(make_powered_rocker_compass(params), seed)
            if not fp.converged:
                print(f"R={R} P={P}: cycle not found, skip")
                continue
            res = basin_slice(
                make_powered_rocker_compass,
                params,
                fp.y,
                axes=(0, 1),
                half_widths=(0.16, 0.65),
                resolution=args.basin_res,
                model_name=f"R{R}_P{P}",
                n_workers=args.workers,
            )
            fracs[(R, P)] = res.basin_fraction
            grids[(R, P)] = res
            seed = fp.y
            print(f"R={R} P={P}: basin_fraction={res.basin_fraction:.3f}")

    run_dir = new_run_dir("active_basin")

    fig, ax = plt.subplots(figsize=(7, 5))
    for P in P_VALUES:
        xs = [R for R in R_VALUES if (R, P) in fracs]
        ys = [fracs[(R, P)] for R in xs]
        ax.plot(xs, ys, "o-", label=f"push_off={P}")
    ax.set_xlabel("R (foot radius)")
    ax.set_ylabel("basin_fraction")
    ax.set_title("Active vs passive basin: basin(R) for each push-off (gamma=0.030)")
    ax.legend()
    fig.savefig(run_dir / "active_basin_curves.png", dpi=150)
    plt.close(fig)

    cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])
    rows = [R for R in R_VALUES if any((R, P) in grids for P in P_VALUES)]
    if rows:
        nP = len(P_VALUES)
        nR = len(rows)
        fig2, axs = plt.subplots(nR, nP, figsize=(3 * nP, 3 * nR))
        axs = np.atleast_2d(axs)
        for i, R in enumerate(rows):
            for j, P in enumerate(P_VALUES):
                a = axs[i, j]
                if (R, P) not in grids:
                    a.axis("off")
                    continue
                g = grids[(R, P)]
                a.imshow(
                    g.grid,
                    origin="lower",
                    extent=[g.ax0_vals[0], g.ax0_vals[-1], g.ax1_vals[0], g.ax1_vals[-1]],
                    aspect="auto",
                    cmap=cmap,
                    vmin=0,
                    vmax=2,
                    interpolation="nearest",
                )
                a.set_title(f"R={R} P={P}\nfrac={g.basin_fraction:.3f}", fontsize=8)
        fig2.tight_layout()
        fig2.savefig(run_dir / "active_basin_montage.png", dpi=150)
        plt.close(fig2)

    summary = {}
    for P in P_VALUES:
        ys = {R: fracs[(R, P)] for R in R_VALUES if (R, P) in fracs}
        if not ys:
            continue
        vals = np.array(list(ys.values()))
        argmax_R = max(ys, key=ys.get)
        cv = float(np.std(vals) / np.mean(vals)) if np.mean(vals) > 0 else None
        summary[P] = {"argmax_R": argmax_R, "cv_across_R": cv, "fractions": ys}

    data = {
        "gamma": 0.030,
        "R_values": R_VALUES,
        "P_values": P_VALUES,
        "basin_fraction": {f"R{R}_P{P}": fracs[(R, P)] for (R, P) in fracs},
        "summary_by_pushoff": {str(P): summary[P] for P in summary},
    }
    (run_dir / "active_basin.json").write_text(json.dumps(data, indent=2, default=str))
    print(
        "summary:",
        {
            P: (summary[P]["argmax_R"], round(summary[P]["cv_across_R"], 3))
            for P in summary
            if summary[P]["cv_across_R"] is not None
        },
    )
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
