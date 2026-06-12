# tests/test_references.py
from crane import references as ref


def test_references_are_filled():
    """文献値が provenance 付きで記入済みであること（Phase 1 ゲートの前提）。

    isinstance 検査なのは、未記入の Ellipsis (...) が `is not None` を
    すり抜けるため。
    """
    assert ref.PROVENANCE.startswith("Garcia")
    assert "http" in ref.PROVENANCE
    assert ref.GAMMA_REF == 0.009
    assert isinstance(ref.LONG_PERIOD_THETA, float)
    assert isinstance(ref.LONG_PERIOD_THETA_DOT, float)
    assert 0.0 < ref.LONG_PERIOD_THETA < 0.5
    assert isinstance(ref.SHORT_PERIOD_THETA, float)
    assert isinstance(ref.SHORT_PERIOD_THETA_DOT, float)
    assert isinstance(ref.STABLE_GAMMA_MAX, float)
