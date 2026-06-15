"""Phase 3 デモ: 点足 kneed walker のリミットサイクル発見 → 多歩シミュ → アニメ生成。

Usage: uv run python scripts/walk_kneed.py [--gamma-deg <published>] [--strides 30] [--perturb 0.005]
"""

import argparse
import json
import sys

import numpy as np

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, make_kneed
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_kneed, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma-deg", type=float, default=np.rad2deg(ref.GAMMA_GAIT))
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.005)
    args = parser.parse_args()

    gamma = np.deg2rad(args.gamma_deg)
    p = KneedParams(
        m_h=ref.M_HIP,
        m_t=ref.M_THIGH,
        m_s=ref.M_SHANK,
        l_t=ref.L_THIGH,
        l_s=ref.L_SHANK,
        b_t=ref.B_THIGH,
        b_s=ref.B_SHANK,
        gamma=gamma,
        g=ref.G,
    )
    model = make_kneed(p)
    run_dir = new_run_dir(f"kneed_g{args.gamma_deg:g}deg")

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
    flexion_maxes = []
    for i in range(args.strides):
        try:
            result = stride(model, x)
        except StrideError as e:
            print(f"stride {i}: FELL ({e})")
            break
        deviation = float(np.linalg.norm(model.project(result.x_end) - fp.y))
        flexion_max = float(np.max(np.abs(result.x[2] - result.x[1])))
        flexion_maxes.append(flexion_max)
        print(
            f"stride {i}: t={result.t_step:.4f} deviation={deviation:.3e} "
            f"knee_flexion_max={np.degrees(flexion_max):.1f}deg"
        )
        strides.append(result)
        x = result.x_end

    plot_phase_portrait(strides, run_dir / "phase_portrait.png")
    animate_kneed(strides, p.l_t, p.l_s, gamma, run_dir / "walk.mp4")

    meta = {
        "params": {
            "m_h": p.m_h,
            "m_t": p.m_t,
            "m_s": p.m_s,
            "l_t": p.l_t,
            "l_s": p.l_s,
            "b_t": p.b_t,
            "b_s": p.b_s,
            "gamma": p.gamma,
            "g": p.g,
        },
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0),
        "n_strides_completed": len(strides),
        "knee_flexion_max_deg": float(np.degrees(max(flexion_maxes))) if flexion_maxes else None,
        "perturb": args.perturb,
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
