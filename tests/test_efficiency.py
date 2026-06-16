import numpy as np

from crane.efficiency import mechanical_cot, relative_loss, step_collision_loss
from crane.models.rocker_compass import RockerCompassParams, kinetic_energy, make_rocker_compass
from crane.search import find_limit_cycle
from crane.stride import stride
from crane import references_mcgeer as ref


def test_relative_loss_basic():
    assert relative_loss(1.0, 0.75) == 0.25
    assert relative_loss(2.0, 2.0) == 0.0


def test_mechanical_cot_basic():
    assert np.isclose(mechanical_cot(0.5, m=1.0, g=10.0, step_length=0.5), 0.1)


def test_step_collision_loss_on_rocker_limit_cycle():
    p = RockerCompassParams(
        m=ref.M_LEG,
        m_h=ref.M_HIP,
        c=ref.C_HIP_TO_COM,
        rho=ref.RHO_GYR,
        L=ref.L_LEG,
        R=ref.R_FOOT,
        gamma=ref.GAMMA_GAIT,
        g=ref.G,
    )
    model = make_rocker_compass(p)
    fp = find_limit_cycle(model, np.array(ref.SECTION_GUESS))
    assert fp.converged
    result = stride(model, model.lift(fp.y))
    loss, ke_pre, ke_post = step_collision_loss(
        result.x_strike, result.x_end, lambda x: kinetic_energy(x, p)
    )
    assert ke_pre > ke_post > 0.0
    assert loss > 0.0
    delta = relative_loss(ke_pre, ke_post)
    assert 0.0 < delta < 1.0
