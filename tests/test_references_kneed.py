# tests/test_references_kneed.py（TDD: 先に書く）
from crane import references_kneed as ref


def test_kneed_references_are_filled():
    """文献値が provenance 付きで記入済み（isinstance は Ellipsis 対策）。"""
    assert "http" in ref.PROVENANCE
    for name in ["M_HIP", "M_THIGH", "M_SHANK", "L_THIGH", "L_SHANK", "B_THIGH", "B_SHANK"]:
        assert isinstance(getattr(ref, name), float), name
    assert isinstance(ref.GAMMA_GAIT, float)
    assert 0.0 < ref.GAMMA_GAIT < 0.3
    assert len(ref.SECTION_GUESS) == 4
