import numpy as np

from crane.basin import BasinResult
from crane.viz import plot_basin


def test_plot_basin_writes_png(tmp_path):
    grid = np.array([[0, 1, 2], [0, 0, 1], [2, 1, 0]], dtype=int)
    res = BasinResult(
        grid=grid,
        ax0_vals=np.linspace(-0.1, 0.1, 3),
        ax1_vals=np.linspace(-0.1, 0.1, 3),
        axes=(0, 1),
        fixed_point=np.array([0.0, 0.0]),
        basin_fraction=float(np.mean(grid == 0)),
        model_name="dummy",
    )
    out = tmp_path / "basin.png"
    plot_basin(res, out)
    assert out.exists() and out.stat().st_size > 0
