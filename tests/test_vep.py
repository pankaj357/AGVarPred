from agvarpred_core.vep import VEP_IDX, parse_vep


def test_parse_vep_missense():
    parts = [""] * 65
    parts[VEP_IDX["Allele"]] = "T"
    parts[VEP_IDX["Consequence"]] = "missense_variant"
    parts[VEP_IDX["IMPACT"]] = "MODERATE"
    parts[VEP_IDX["SYMBOL"]] = "BRCA1"
    parts[VEP_IDX["SIFT"]] = "deleterious(0.01)"
    parts[VEP_IDX["PolyPhen"]] = "probably_damaging(0.999)"
    parts[VEP_IDX["Protein_position"]] = "61"
    csq = "|".join(parts)

    result = parse_vep(csq)
    assert result["vep_is_missense"] == 1
    assert result["vep_is_synonymous"] == 0
    assert result["vep_IMPACT_score"] == 3
    assert result["vep_SIFT_score"] == 0.01
    assert result["vep_PolyPhen_score"] == 0.999
    assert result["vep_Protein_position"] == 61


def test_parse_vep_empty():
    assert parse_vep("") == {}
    assert parse_vep(None) == {}
