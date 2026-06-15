# tests/test_references_mcgeer.py（TDD: 先に書く）
from crane import references_mcgeer as ref


def test_mcgeer_references_are_filled():
    assert "http" in ref.PROVENANCE
    for name in ["M_LEG", "M_HIP", "C_HIP_TO_COM", "RHO_GYR", "L_LEG", "R_FOOT"]:
        assert isinstance(getattr(ref, name), float), name
    assert isinstance(ref.GAMMA_GAIT, float)
    assert 0.0 < ref.GAMMA_GAIT < 0.3
    assert len(ref.SECTION_GUESS) == 3  # (θ_st, θ̇_st, θ̇_sw)
    assert isinstance(ref.PUBLISHED_STABLE, bool)
    # §5 直脚モデルには印刷固有値が無い（Table 1 は膝あり test machine）— None を固定
    assert ref.POINCARE_EIGENVALUES is None
