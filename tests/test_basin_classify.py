import numpy as np

from crane.basin import CONVERGED, FELL, classify_ic
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane import references as ref


def _simplest_fp():
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    return model, fp.y


def test_fixed_point_classifies_converged():
    model, y_star = _simplest_fp()
    assert classify_ic(model, y_star, y_star, max_strides=20) == CONVERGED


def test_far_point_classifies_fell():
    model, y_star = _simplest_fp()
    far = y_star + np.array([0.5, 0.5])  # basin 外（転倒）
    assert classify_ic(model, far, y_star, max_strides=20) == FELL
