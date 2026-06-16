"""Phase 4a: 単体モデルの basin of attraction 地図を生成。

Usage: uv run python scripts/basin_map.py <model> [--resolution 120] [--workers N]
  model: simplest | compass | kneed | rocker_compass | rocker_kneed
"""

import argparse
import json
import sys

import numpy as np

from crane.basin import CONVERGED, FELL, UNDECIDED, basin_slice
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.viz import plot_basin


def _registry():
    """name -> (make_fn, params, guess, axes, half_widths) を返す。

    half_widths は実装時に basin を窓に収めるよう較正（収束領域が窓端に
    接していたら広げる）。下記は初期値。
    """
    from crane import references as ref_s
    from crane import references_goswami as ref_c
    from crane import references_kneed as ref_k
    from crane import references_mcgeer as ref_rc
    from crane import references_mcgeer_knees as ref_rk
    from crane.models.compass import CompassParams, make_compass
    from crane.models.kneed import KneedParams, make_kneed
    from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
    from crane.models.rocker_kneed import RockerKneedParams, make_rocker_kneed
    from crane.models.simplest import SimplestParams, make_simplest

    return {
        "simplest": (
            make_simplest,
            SimplestParams(gamma=ref_s.GAMMA_REF),
            np.array([ref_s.LONG_PERIOD_THETA, ref_s.LONG_PERIOD_THETA_DOT]),
            (0, 1),
            (0.08, 0.08),
        ),
        "compass": (
            make_compass,
            CompassParams(
                m=ref_c.M_LEG,
                m_h=ref_c.M_HIP,
                a=ref_c.A,
                b=ref_c.B,
                gamma=ref_c.GAMMA_GAIT,
                g=ref_c.G,
            ),
            np.array(ref_c.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
        "kneed": (
            make_kneed,
            KneedParams(
                m_h=ref_k.M_HIP,
                m_t=ref_k.M_THIGH,
                m_s=ref_k.M_SHANK,
                l_t=ref_k.L_THIGH,
                l_s=ref_k.L_SHANK,
                b_t=ref_k.B_THIGH,
                b_s=ref_k.B_SHANK,
                gamma=ref_k.GAMMA_GAIT,
                g=ref_k.G,
            ),
            np.array(ref_k.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
        "rocker_compass": (
            make_rocker_compass,
            RockerCompassParams(
                m=ref_rc.M_LEG,
                m_h=ref_rc.M_HIP,
                c=ref_rc.C_HIP_TO_COM,
                rho=ref_rc.RHO_GYR,
                L=ref_rc.L_LEG,
                R=ref_rc.R_FOOT,
                gamma=ref_rc.GAMMA_GAIT,
                g=ref_rc.G,
            ),
            np.array(ref_rc.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
        "rocker_kneed": (
            make_rocker_kneed,
            RockerKneedParams(
                m_h=ref_rk.M_HIP,
                m_t=ref_rk.M_THIGH,
                m_s=ref_rk.M_SHANK,
                l_t=ref_rk.L_THIGH,
                l_s=ref_rk.L_SHANK,
                b_t=ref_rk.B_THIGH,
                b_s=ref_rk.B_SHANK,
                R=ref_rk.R_FOOT,
                gamma=ref_rk.GAMMA_GAIT,
                g=ref_rk.G,
            ),
            np.array(ref_rk.SECTION_GUESS),
            (0, 1),
            (0.12, 0.5),
        ),
    }


def main() -> None:
    reg = _registry()
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=sorted(reg))
    parser.add_argument("--resolution", type=int, default=120)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    make_fn, params, guess, axes, half_widths = reg[args.model]
    model = make_fn(params)
    fp = find_limit_cycle(model, guess)
    if not fp.converged:
        sys.exit(f"ERROR: limit cycle not found for {args.model}")
    print(f"{args.model}: y*={fp.y}  |lambda|={np.abs(fp.eigenvalues)}")

    res = basin_slice(
        make_fn,
        params,
        fp.y,
        axes=axes,
        half_widths=half_widths,
        resolution=args.resolution,
        model_name=args.model,
        n_workers=args.workers,
    )
    run_dir = new_run_dir(f"basin_{args.model}")
    plot_basin(res, run_dir / "basin.png")
    meta = {
        "model": args.model,
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": np.abs(fp.eigenvalues).tolist(),
        "axes": list(axes),
        "half_widths": list(half_widths),
        "resolution": args.resolution,
        "basin_fraction": res.basin_fraction,
        "counts": {
            "converged": int(np.sum(res.grid == CONVERGED)),
            "fell": int(np.sum(res.grid == FELL)),
            "undecided": int(np.sum(res.grid == UNDECIDED)),
        },
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"basin_fraction={res.basin_fraction:.3f}  outputs -> {run_dir}")


if __name__ == "__main__":
    main()
