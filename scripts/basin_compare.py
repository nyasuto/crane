# scripts/basin_compare.py
"""Phase 4a: 5モデルの basin を比較。montage + basin_fraction 表を出力。

統制ペア（point-foot vs rocker-foot）で「円弧足が basin を広げるか」を検証する。
仮説が偽でも結果を正直に報告する（強制しない）。
Usage: uv run python scripts/basin_compare.py [--resolution 120] [--workers N]
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from crane.basin import basin_slice  # noqa: E402
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle  # noqa: E402

from basin_map import _registry  # noqa: E402  同 scripts/ 内


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=120)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    reg = _registry()
    order = ["simplest", "compass", "rocker_compass", "kneed", "rocker_kneed"]
    results = {}
    for name in order:
        make_fn, params, guess, axes, half_widths = reg[name]
        fp = find_limit_cycle(make_fn(params), guess)
        if not fp.converged:
            raise SystemExit(f"ERROR: limit cycle not found for {name}")
        res = basin_slice(
            make_fn,
            params,
            fp.y,
            axes=axes,
            half_widths=half_widths,
            resolution=args.resolution,
            model_name=name,
            n_workers=args.workers,
        )
        results[name] = res
        print(f"{name}: basin_fraction={res.basin_fraction:.3f}")

    # montage
    cmap = ListedColormap(["#2ca02c", "#d62728", "#999999"])
    fig, axs = plt.subplots(1, len(order), figsize=(4 * len(order), 4))
    for ax, name in zip(axs, order):
        r = results[name]
        ax.imshow(
            r.grid,
            origin="lower",
            extent=[r.ax0_vals[0], r.ax0_vals[-1], r.ax1_vals[0], r.ax1_vals[-1]],
            aspect="auto",
            cmap=cmap,
            vmin=0,
            vmax=2,
            interpolation="nearest",
        )
        a0, a1 = r.axes
        ax.plot(r.fixed_point[a0], r.fixed_point[a1], "k*", markersize=10)
        ax.set_title(f"{name}\nfrac={r.basin_fraction:.3f}")
    fig.tight_layout()
    run_dir = new_run_dir("basin_compare")
    fig.savefig(run_dir / "basin_compare.png", dpi=150)
    plt.close(fig)

    # 方向性主張: rocker > point
    table = {name: results[name].basin_fraction for name in order}
    claims = {
        "compass_rocker_gt_point": table["rocker_compass"] > table["compass"],
        "kneed_rocker_gt_point": table["rocker_kneed"] > table["kneed"],
    }
    (run_dir / "basin_fractions.json").write_text(
        json.dumps({"fractions": table, "claims": claims}, indent=2)
    )
    print("claims:", claims)
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
