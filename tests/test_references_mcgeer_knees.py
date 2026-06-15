from crane import references_mcgeer_knees as ref


def test_mcgeer_knees_references_are_filled():
    assert "http" in ref.PROVENANCE
    for name in [
        "M_HIP",
        "M_THIGH",
        "M_SHANK",
        "L_THIGH",
        "L_SHANK",
        "B_THIGH",
        "B_SHANK",
        "R_FOOT",
    ]:
        assert isinstance(getattr(ref, name), float), name
    assert isinstance(ref.GAMMA_GAIT, float)
    assert 0.0 < ref.GAMMA_GAIT < 0.3
    assert len(ref.SECTION_GUESS) == 3  # (θ_st, θ̇_st, θ̇_sw)
    assert isinstance(ref.PUBLISHED_STABLE, bool)
