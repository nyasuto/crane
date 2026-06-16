"""Phase 5b デモ: 動力付き rocker_compass（γ=0.030 + push-off 増強）の動力サイクル → 多歩 → アニメ。

Usage: uv run python scripts/walk_powered_rocker.py [--push-off 0.08] [--strides 30] [--perturb 0.005]
"""

import argparse
import json
import sys

import numpy as np

from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_rocker, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--push-off", type=float, default=0.08)
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.005)
    args = parser.parse_args()

    p = PoweredRockerCompassParams(
        m=1.0,
        m_h=0.0,
        c=0.37,
        rho=0.32,
        L=1.0,
        R=0.3,
        gamma=0.030,
        push_off=args.push_off,
        g=9.81,
    )
    model = make_powered_rocker_compass(p)
    fp = find_limit_cycle(model, np.array([0.30844, -1.26256, -0.87914]))
    if not fp.converged:
        sys.exit("ERROR: 動力サイクルが見つからない（--push-off を調整）")
    lam = np.abs(fp.eigenvalues)
    print(f"converged (gamma=0.030, push_off={args.push_off}): y*={fp.y} |lambda|={lam}")

    run_dir = new_run_dir(f"powered_rocker_P{args.push_off:g}")
    y0 = fp.y * (1.0 + args.perturb)
    x = model.lift(y0)
    strides = []
    for i in range(args.strides):
        try:
            result = stride(model, x)
        except StrideError as e:
            print(f"stride {i}: FELL ({e})")
            break
        deviation = float(np.linalg.norm(model.project(result.x_end) - fp.y))
        print(f"stride {i}: t={result.t_step:.4f} deviation={deviation:.3e}")
        strides.append(result)
        x = result.x_end

    plot_phase_portrait(strides, run_dir / "phase_portrait.png")
    animate_rocker(strides, p.L, p.R, p.gamma, run_dir / "walk.mp4")

    meta = {
        "gamma": 0.030,
        "push_off": args.push_off,
        "R": p.R,
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0),
        "n_strides_completed": len(strides),
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
