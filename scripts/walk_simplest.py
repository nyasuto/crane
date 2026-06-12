"""Phase 1 デモ: リミットサイクル発見 → 摂動から多歩シミュ → アニメ生成。

Usage: uv run python scripts/walk_simplest.py [--gamma 0.009] [--strides 30] [--perturb 0.01]
"""

import argparse
import json
import sys

import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams, lift
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_walk, plot_phase_portrait


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma", type=float, default=ref.GAMMA_REF)
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.01, help="不動点への相対摂動")
    args = parser.parse_args()

    p = SimplestParams(gamma=args.gamma)
    run_dir = new_run_dir(f"simplest_g{args.gamma:g}")

    fp = find_limit_cycle(p, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    print(
        f"converged={fp.converged}  y*={fp.y}  |lambda|={np.abs(fp.eigenvalues) if fp.eigenvalues is not None else None}"
    )
    for i, (y, r) in enumerate(fp.history):
        print(f"  newton[{i}] y={y} residual={r:.3e}")

    if not fp.converged:
        print(
            "ERROR: Newton 法が収束しませんでした。"
            f" --gamma {args.gamma} が基準推測のニュートン域外の可能性があります。"
            " 収束履歴を確認してください。",
            file=sys.stderr,
        )
        sys.exit(1)

    # 不動点からわずかに摂動した初期条件で多歩シミュ（basin 内なら収束するはず）
    y0 = fp.y * (1.0 + args.perturb)
    x = lift(y0[0], y0[1])
    strides = []
    for i in range(args.strides):
        try:
            result = stride(p, x)
        except StrideError as e:
            print(f"stride {i}: FELL ({e})")
            break
        deviation = float(np.linalg.norm(np.array([result.x_end[0], result.x_end[2]]) - fp.y))
        print(f"stride {i}: t={result.t_step:.4f} deviation={deviation:.3e}")
        strides.append(result)
        x = result.x_end

    plot_phase_portrait(strides, run_dir / "phase_portrait.png")
    animate_walk(strides, p.gamma, run_dir / "walk.mp4")

    meta = {
        "params": {"gamma": p.gamma},
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": np.abs(fp.eigenvalues).tolist(),
        "n_strides_completed": len(strides),
        "perturb": args.perturb,
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
