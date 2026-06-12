"""slope continuation で compass gait family を追跡し分岐図を描く。

Usage: uv run python scripts/bifurcation_compass.py [--gamma-max-deg 6.0] [--step-deg 0.05]
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from crane import references_goswami as ref  # noqa: E402
from crane.models.compass import CompassParams, make_compass  # noqa: E402
from crane.runs import new_run_dir  # noqa: E402
from crane.search import find_limit_cycle, poincare_map  # noqa: E402


def make_params(gamma: float) -> CompassParams:
    return CompassParams(m=ref.M_LEG, m_h=ref.M_HIP, a=ref.A, b=ref.B, gamma=gamma, g=ref.G)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma-max-deg", type=float, default=6.0)
    parser.add_argument("--step-deg", type=float, default=0.05)
    args = parser.parse_args()

    run_dir = new_run_dir("bifurcation_compass")
    rows = []

    # period-1 branch: 公表 slope から上向きに continuation
    y = np.array(ref.SECTION_GUESS)
    gamma = ref.GAMMA_GAIT
    step = np.deg2rad(args.step_deg)
    y2 = None  # period-2 branch の seed
    y2_gamma = None  # 不安定化点の gamma
    while gamma <= np.deg2rad(args.gamma_max_deg):
        fp = find_limit_cycle(make_compass(make_params(gamma)), y)
        if not fp.converged:
            print(f"gamma={np.rad2deg(gamma):.3f}deg: period-1 lost")
            break
        lam = float(np.max(np.abs(fp.eigenvalues)))
        rows.append([gamma, 1, fp.y[0], lam])
        print(f"gamma={np.rad2deg(gamma):.3f}deg p1 theta*={fp.y[0]:.5f} max|l|={lam:.4f}")
        if lam > 1.0 and y2 is None:
            y2 = fp.y.copy()  # 不安定化を検知 → ここから period-2 を試す
            y2_gamma = gamma
        y = fp.y
        gamma += step

    # period-2 branch（period-1 不安定化点から上向きに continuation）
    # 1e-3 摂動で収束しない場合は S 写像を数回反復してからシードする（plan の指示）
    if y2 is not None and y2_gamma is not None:
        gamma2 = y2_gamma
        # まず不安定化点で period-2 シードを確立する
        model2 = make_compass(make_params(gamma2))
        # 1e-3 摂動を試し、period-1 と同じ点に落ちる場合は S 写像反復で脱出する
        fp_seed = None
        for eps in [1e-3, 1e-2]:
            y_seed = y2 * (1.0 + eps)
            fp_seed = find_limit_cycle(model2, y_seed, n_strides=2)
            if fp_seed.converged and not np.allclose(fp_seed.y, y2, atol=1e-6):
                break
            # S 写像を最大 20 回反復してからシード
            y_iter = y_seed.copy()
            for _ in range(20):
                try:
                    y_iter = poincare_map(model2, y_iter)
                except Exception:
                    break
            fp_seed = find_limit_cycle(model2, y_iter, n_strides=2)
            if fp_seed.converged and not np.allclose(fp_seed.y, y2, atol=1e-6):
                break
        if fp_seed is not None and fp_seed.converged and not np.allclose(fp_seed.y, y2, atol=1e-6):
            y = fp_seed.y
            # y2 を更新して continuation を開始
            while gamma2 <= np.deg2rad(args.gamma_max_deg):
                fp = find_limit_cycle(make_compass(make_params(gamma2)), y, n_strides=2)
                if not fp.converged:
                    print(f"gamma={np.rad2deg(gamma2):.3f}deg: period-2 lost")
                    break
                lam = float(np.max(np.abs(fp.eigenvalues)))
                rows.append([gamma2, 2, fp.y[0], lam])
                print(f"gamma={np.rad2deg(gamma2):.3f}deg p2 theta*={fp.y[0]:.5f} max|l|={lam:.4f}")
                y = fp.y
                gamma2 += step
        else:
            print("period-2 seed not found near bifurcation point")

    with (run_dir / "bifurcation.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gamma_rad", "period", "theta_star", "max_abs_lambda"])
        writer.writerows(rows)

    data = np.array([[r[0], r[1], r[2], r[3]] for r in rows])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    for period, marker in [(1, "o"), (2, "s")]:
        sel = data[data[:, 1] == period]
        if sel.size:
            stable = sel[sel[:, 3] < 1.0]
            unstable = sel[sel[:, 3] >= 1.0]
            ax1.plot(
                np.rad2deg(stable[:, 0]), stable[:, 2], marker, ms=3, label=f"p{period} stable"
            )
            if unstable.size:
                ax1.plot(
                    np.rad2deg(unstable[:, 0]),
                    unstable[:, 2],
                    marker,
                    ms=3,
                    mfc="none",
                    label=f"p{period} unstable",
                )
    ax1.set_ylabel("theta* [rad]")
    ax1.legend()
    sel1 = data[data[:, 1] == 1]
    ax2.plot(np.rad2deg(sel1[:, 0]), sel1[:, 3], "-")
    ax2.axhline(1.0, color="k", lw=0.5)
    ax2.set_xlabel("slope [deg]")
    ax2.set_ylabel("max|lambda| (period-1)")
    fig.savefig(run_dir / "bifurcation.png", dpi=150)
    print(f"outputs -> {run_dir}")


if __name__ == "__main__":
    main()
