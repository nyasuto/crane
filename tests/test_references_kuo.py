from crane import references_kuo as ref


def test_kuo_provenance_present():
    assert isinstance(ref.PROVENANCE, str) and len(ref.PROVENANCE) > 0
    assert (
        ("http" in ref.PROVENANCE) or ("未取得" in ref.PROVENANCE) or ("取得失敗" in ref.PROVENANCE)
    )
