"""Phase 5a デモ: 動力付き simplest walker の平地リミットサイクル → 多歩シミュ → アニメ。

平地 γ=0 では直接 shoot が静止解に落ちるため、受動 γ=0.009 サイクルから
γ→0・push_off→target へ continuation して平地サイクルを求める。
Usage: uv run python scripts/walk_powered.py [--push-off 0.115] [--strides 30] [--perturb 0.01]
"""

import argparse
import json
import sys

import numpy as np

from crane.models.powered_simplest import PoweredSimplestParams, make_powered_simplest
from crane.runs import new_run_dir
from crane.search import find_limit_cycle
from crane.stride import StrideError, stride
from crane.viz import animate_walk, plot_phase_portrait


def find_level_cycle(push_off_target: float, n_steps: int = 6):
    """受動 γ=0.009 から γ→0・push_off→target へ continuation。

    各ステップで前解を seed に Newton。静止解（|θ|<0.05）は棄却。
    最終 (γ=0, push_off=target) の model と FixedPoint を返す。失敗時 (None, None)。
    """
    gammas = np.linspace(0.009, 0.0, n_steps)
    pushes = np.linspace(0.0, push_off_target, n_steps)
    guess = np.array([0.2003109, -0.1998325])
    model = fp = None
    for g, P in zip(gammas, pushes):
        model = make_powered_simplest(PoweredSimplestParams(gamma=float(g), push_off=float(P)))
        fp = find_limit_cycle(model, guess)
        if not fp.converged or abs(fp.y[0]) < 0.05:
            return None, None
        guess = fp.y
    return model, fp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--push-off", type=float, default=0.115)
    parser.add_argument("--strides", type=int, default=30)
    parser.add_argument("--perturb", type=float, default=0.01)
    args = parser.parse_args()

    model, fp = find_level_cycle(args.push_off)
    if model is None or fp is None or not fp.converged:
        sys.exit("ERROR: 平地動力サイクルが見つからない（--push-off を調整）")
    lam = np.abs(fp.eigenvalues)
    print(f"converged (gamma=0, push_off={args.push_off}): y*={fp.y}  |lambda|={lam}")
    if lam.max() > 1.0:
        print(f"WARNING: UNSTABLE (max|lambda|={lam.max():.3f})")

    run_dir = new_run_dir(f"powered_P{args.push_off:g}")
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
    animate_walk(strides, 0.0, run_dir / "walk.mp4")

    meta = {
        "gamma": 0.0,
        "push_off": args.push_off,
        "pushoff_work_per_step": 0.5 * args.push_off**2,
        "fixed_point": fp.y.tolist(),
        "eigenvalues_abs": lam.tolist(),
        "stable": bool(lam.max() < 1.0),
        "n_strides_completed": len(strides),
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
