import numpy as np

from crane.basin import basin_slice
from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.search import find_limit_cycle

NOM = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)
GUESS = np.array([0.30, -1.0, -0.40])


def test_powered_pushoff_zero_basin_matches_passive():
    """P=0 の能動 basin が受動 rocker_compass の basin に一致（内部ゲート、低解像度）。"""
    passive_params = RockerCompassParams(**NOM)
    powered_params = PoweredRockerCompassParams(**NOM, push_off=0.0)

    fp_passive = find_limit_cycle(make_rocker_compass(passive_params), GUESS)
    fp_powered = find_limit_cycle(make_powered_rocker_compass(powered_params), GUESS)
    assert fp_passive.converged and fp_powered.converged
    np.testing.assert_allclose(fp_powered.y, fp_passive.y, atol=1e-8)

    kw = dict(axes=(0, 1), half_widths=(0.16, 0.65), resolution=15, n_workers=1)
    b_passive = basin_slice(
        make_rocker_compass, passive_params, fp_passive.y, model_name="passive", **kw
    )
    b_powered = basin_slice(
        make_powered_rocker_compass, powered_params, fp_powered.y, model_name="powered", **kw
    )
    assert np.array_equal(b_powered.grid, b_passive.grid)
    assert b_powered.basin_fraction == b_passive.basin_fraction
