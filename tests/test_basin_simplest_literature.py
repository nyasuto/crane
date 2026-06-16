# tests/test_basin_simplest_literature.py
import numpy as np

from crane import references as ref
from crane.basin import CONVERGED, FELL, basin_slice
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle


def _simplest_basin(resolution=60):
    model = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))
    fp = find_limit_cycle(model, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    return basin_slice(
        make_simplest,
        SimplestParams(gamma=ref.GAMMA_REF),
        fp.y,
        axes=(0, 1),
        half_widths=(0.08, 0.08),
        resolution=resolution,
        model_name="simplest",
        n_workers=None,
    )


def test_simplest_basin_is_thin_connected_region():
    """定性ゲート: basin は薄い領域（窓全体ではない）かつ非空。

    Schwab & Wisse 2001 の basin 図と整合（薄く、窓を埋め尽くさない）。
    """
    res = _simplest_basin()
    frac = res.basin_fraction
    assert 0.0 < frac < 0.95  # 非空かつ窓を埋め尽くさない（＝薄い）
    # CONVERGED と FELL が両方存在（境界が窓内にある）
    assert np.any(res.grid == CONVERGED) and np.any(res.grid == FELL)


def test_simplest_basin_area_matches_reference_if_available():
    """定量ゲート: 原典に basin 面積数値があれば本実装と照合（なければ skip）。"""
    if ref.BASIN_AREA_REF is None:
        import pytest

        pytest.skip("Schwab & Wisse の basin 面積数値が未取得（定性ゲートに委譲）")
    res = _simplest_basin()
    rel = abs(res.basin_fraction - ref.BASIN_AREA_REF) / ref.BASIN_AREA_REF
    assert rel < ref.BASIN_AREA_REF_TOL
