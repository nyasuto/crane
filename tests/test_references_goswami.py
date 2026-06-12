# tests/test_references_goswami.py（先に書く: TDD）
from crane import references_goswami as ref


def test_goswami_references_are_filled():
    """文献値が provenance 付きで記入済み（isinstance 検査は Ellipsis 対策）。"""
    assert ref.PROVENANCE.startswith("Goswami")
    assert "http" in ref.PROVENANCE
    assert isinstance(ref.M_LEG, float)
    assert isinstance(ref.M_HIP, float)
    assert isinstance(ref.A, float)
    assert isinstance(ref.B, float)
    assert isinstance(ref.GAMMA_GAIT, float)  # 公表 gait のある slope [rad]
    assert 0.0 < ref.GAMMA_GAIT < 0.2
