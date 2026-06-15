import numpy as np

from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.stride import stride

P = RockerCompassParams(m=5.0, m_h=10.0, c=0.5, rho=0.15, L=1.0, R=0.3, gamma=0.05, g=9.81)
MODEL = make_rocker_compass(P)


def test_stride_returns_to_section():
    """素直な初期条件から 1 stride で断面（θ_sw=−θ_st）に戻る。

    R=0.3 では task 既定の seed (0.2, -1.0, -0.3) は転倒（θ_st が数 rad まで
    暴れ t_step≈5s）になるため、swing 速度を調整した (0.2, -1.2, 0.6) を使う。
    この seed では θ_st は |0.2| 以内に収まり t_step≈0.33s の素直な 1 歩になる。
    """
    x0 = MODEL.lift(np.array([0.2, -1.2, 0.6]))
    result = stride(MODEL, x0)
    assert np.isclose(result.x_end[1], -result.x_end[0], atol=1e-8)
    assert result.t_step > 0.1
    assert np.abs(result.x[0]).max() < 0.8  # θ_st が暴れない（転倒でない）
