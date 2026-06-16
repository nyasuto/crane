import numpy as np

from crane.models.simplest import heelstrike_map as passive_heelstrike
from crane.models.powered_simplest import (
    PoweredSimplestParams,
    make_powered_simplest,
    powered_heelstrike_map,
)


def test_pushoff_zero_equals_passive_map():
    for x in [
        np.array([-0.2, -0.4, -0.2, 0.1]),
        np.array([-0.15, -0.3, -0.25, 0.05]),
    ]:
        np.testing.assert_allclose(
            powered_heelstrike_map(x, 0.0), passive_heelstrike(x), atol=1e-14
        )


def test_pushoff_adds_sin2theta_P_to_stance_velocity():
    x = np.array([-0.2, -0.4, -0.2, 0.1])
    P = 0.1
    out0 = powered_heelstrike_map(x, 0.0)
    outP = powered_heelstrike_map(x, P)
    s = np.sin(2.0 * x[0])
    assert np.isclose(outP[2] - out0[2], s * P, atol=1e-14)
    c = np.cos(2.0 * x[0])
    assert np.isclose(outP[3], (1.0 - c) * outP[2], atol=1e-14)


def test_factory_builds_model_with_passive_dynamics():
    p = PoweredSimplestParams(gamma=0.0, push_off=0.1)
    model = make_powered_simplest(p)
    assert len(model.phases) == 1
    y = np.array([0.2, -0.2])
    x = model.lift(y)
    np.testing.assert_allclose(model.project(x), y, atol=1e-14)
