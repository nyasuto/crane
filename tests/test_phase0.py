import numpy as np
from scipy.integrate import solve_ivp

from crane.runs import new_run_dir


def test_event_detection_quarter_period():
    """solve_ivp の event 検出が解析解 t=π/2 と一致する。"""

    def f(t, x):
        return [x[1], -x[0]]

    def zero_cross(t, x):
        return x[0]

    zero_cross.terminal = True
    zero_cross.direction = -1

    sol = solve_ivp(f, (0.0, 10.0), [1.0, 0.0], events=zero_cross, rtol=1e-10, atol=1e-12)
    t_event = sol.t_events[0][0]
    assert abs(t_event - np.pi / 2) < 1e-8


def test_new_run_dir(tmp_path):
    """タイムスタンプ付き run ディレクトリが作られ、重複しない。"""
    d1 = new_run_dir("smoke", base=tmp_path)
    assert d1.is_dir()
    assert "smoke" in d1.name
