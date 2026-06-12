"""Phase 2 デモ: compass gait のリミットサイクル発見 → 多歩シミュ → アニメ生成。

Usage: uv run python scripts/walk_compass.py [--gamma-deg <published>] [--strides 30] [--perturb 0.005]
"""

import argparse
import json
import sys

import numpy as np

from crane import references_goswami as ref
from crane.models.compass import CompassParams, make_compass
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_walk, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma-deg", type=float, default=np.rad2deg(ref.GAMMA_GAIT))
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.005)
    args = parser.parse_args()

    gamma = np.deg2rad(args.gamma_deg)
    p = CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=gamma, g=ref.G)
    model = make_compass(p)
    run_dir = new_run_dir(f"compass_g{args.gamma_deg:g}deg")

    fp = find_limit_cycle(model, np.array(ref.SECTION_GUESS))
    if not fp.converged:
        for i, (y, r) in enumerate(fp.history):
            print(f"  newton[{i}] y={y} residual={r:.3e}")
        sys.exit("ERROR: limit cycle not found")
    lam = np.abs(fp.eigenvalues)
    print(f"converged: y*={fp.y}  |lambda|={lam}")
    if lam.max() > 1.0:
        print(f"WARNING: cycle is UNSTABLE (max|lambda|={lam.max():.3f}) — walk will fall")

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
    animate_walk(strides, gamma, run_dir / "walk.mp4", angles_of=lambda q0, q1: (q0, q1))

    meta = {
        "params": {"m": p.m, "m_h": p.m_h, "a": p.a, "b": p.b, "gamma": p.gamma, "g": p.g},
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0),
        "n_strides_completed": len(strides),
        "perturb": args.perturb,
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
