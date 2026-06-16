import numpy as np

from crane.basin import CONVERGED, FELL, UNDECIDED, classify_ic
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


def test_neither_converged_nor_fell_is_undecided():
    model, y_star = _simplest_fp()
    near = y_star * 1.005  # 固定点から 0.5% 摂動
    # tol が極端に厳しく strides が少ない: 未収束だが転倒もしていない
    assert classify_ic(model, near, y_star, max_strides=2, converge_tol=1e-12) == UNDECIDED
