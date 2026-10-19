"""Focused unit tests for the reduced gnomAD TSV query path."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from agvarpred_core.gnomad import _query_gnomad_tsv


def _vep_string(
    consequence: str = "missense_variant",
    impact: str = "MODERATE",
    sift: str = "deleterious(0.01)",
    polyphen: str = "probably_damaging(0.999)",
    protein_pos: str = "100",
    amino_acids: str = "A/V",
    lof: str = "HC",
) -> str:
    """Build a 70-field VEP CSQ string matching gnomAD r2.1.1 field order."""
    parts = [""] * 70
    # Field indices from agvarpred_core.vep.VEP_FIELDS
    parts[0] = "A"  # Allele
    parts[1] = consequence
    parts[2] = impact
    parts[3] = "GENE1"  # SYMBOL
    parts[14] = protein_pos
    parts[15] = amino_acids
    parts[35] = sift
    parts[36] = polyphen
    parts[53] = lof
    return "|".join(parts)


def _build_synthetic_tsv(tmp_path: Path) -> Path:
    import pysam

    tsv_path = tmp_path / "gnomad_reduced.tsv.gz"
    lines = [
        "#CHROM\tSTART\tEND\tREF\tALT\tAF\tVEP\n",
        # Present SNV (multiallelic record split into two rows)
        f"chr1\t99\t100\tA\tG\t0.1\t{_vep_string(consequence='missense_variant', impact='MODERATE')}\n",
        f"chr1\t99\t100\tA\tT\t0.2\t{_vep_string(consequence='synonymous_variant', impact='LOW')}\n",
        # Deletion
        f"chr1\t199\t200\tAGCG\tA\t0.05\t{_vep_string(consequence='inframe_deletion', impact='MODERATE')}\n",
        # Insertion
        "chr1\t299\t300\tA\tAT\t0.0\t.\n",
        # Missing AF and VEP
        "chr2\t49\t50\tC\tT\t0.0\t.\n",
        # chrM variant
        f"chrM\t99\t100\tA\tG\t0.01\t{_vep_string(consequence='missense_variant', impact='MODERATE')}\n",
    ]

    with pysam.BGZFile(str(tsv_path), "wb") as fh:
        for line in lines:
            fh.write(line.encode("utf-8"))

    subprocess.run(
        ["tabix", "-s", "1", "-b", "2", "-e", "3", "-S", "1", str(tsv_path)],
        check=True,
    )
    return tsv_path


@pytest.fixture
def synthetic_tsv(tmp_path: Path) -> Path:
    return _build_synthetic_tsv(tmp_path)


def test_present_snv(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chr1_100_A_G"], synthetic_tsv)
    assert af_map == {"chr1_100_A_G": pytest.approx(0.1)}
    assert vep_map["chr1_100_A_G"]["vep_Consequence"] == "missense_variant"
    assert vep_map["chr1_100_A_G"]["vep_IMPACT"] == "MODERATE"


def test_multiallelic_second_alt(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chr1_100_A_T"], synthetic_tsv)
    assert af_map == {"chr1_100_A_T": pytest.approx(0.2)}
    assert vep_map["chr1_100_A_T"]["vep_Consequence"] == "synonymous_variant"


def test_deletion(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chr1_200_AGCG_A"], synthetic_tsv)
    assert af_map == {"chr1_200_AGCG_A": pytest.approx(0.05)}
    assert vep_map["chr1_200_AGCG_A"]["vep_Consequence"] == "inframe_deletion"


def test_insertion(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chr1_300_A_AT"], synthetic_tsv)
    assert af_map == {"chr1_300_A_AT": pytest.approx(0.0)}
    assert "chr1_300_A_AT" not in vep_map


def test_missing_af(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chr2_50_C_T"], synthetic_tsv)
    assert af_map == {"chr2_50_C_T": pytest.approx(0.0)}
    assert "chr2_50_C_T" not in vep_map


def test_absent_variant(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chr1_100_A_C"], synthetic_tsv)
    assert af_map == {}
    assert vep_map == {}


def test_chrom_m(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chrM_100_A_G"], synthetic_tsv)
    assert af_map == {"chrM_100_A_G": pytest.approx(0.01)}
    assert vep_map["chrM_100_A_G"]["vep_Consequence"] == "missense_variant"


def test_batch_query(synthetic_tsv: Path):
    ids = ["chr1_100_A_G", "chr1_100_A_T", "chr1_200_AGCG_A", "chr1_300_A_AT"]
    af_map, vep_map = _query_gnomad_tsv(ids, synthetic_tsv)
    assert af_map == {
        "chr1_100_A_G": pytest.approx(0.1),
        "chr1_100_A_T": pytest.approx(0.2),
        "chr1_200_AGCG_A": pytest.approx(0.05),
        "chr1_300_A_AT": pytest.approx(0.0),
    }
    assert set(vep_map.keys()) == {"chr1_100_A_G", "chr1_100_A_T", "chr1_200_AGCG_A"}


def test_unknown_chromosome(synthetic_tsv: Path):
    af_map, vep_map = _query_gnomad_tsv(["chrZ_100_A_G"], synthetic_tsv)
    assert af_map == {}
    assert vep_map == {}


def test_tabix_tsv_gnomad_source(synthetic_tsv: Path):
    from agvarpred_core.af_source import TabixTSVGnomADSource

    source = TabixTSVGnomADSource(str(synthetic_tsv))
    assert source.is_available()
    af_map, vep_map = source.query(["chr1_100_A_G", "chr1_100_A_T"])
    assert af_map == {
        "chr1_100_A_G": pytest.approx(0.1),
        "chr1_100_A_T": pytest.approx(0.2),
    }


def test_resolve_af_source_prefers_tsv_over_vcf(tmp_path: Path, monkeypatch):
    from agvarpred_core.af_source import resolve_af_source

    tsv_path = _build_synthetic_tsv(tmp_path)
    monkeypatch.setenv("GNOMAD_TSV", str(tsv_path))
    # Ensure no VCF env var is set.
    monkeypatch.delenv("GNOMAD_VCF", raising=False)
    source = resolve_af_source()
    assert source.name == "tsv_gnomad"


def test_source_parity_compare_maps_counts():
    """Unit test for the source-parity counting logic."""
    import sys

    scripts_dir = str(Path(__file__).parent.parent / "scripts")
    sys.path.insert(0, scripts_dir)
    from validate_gnomad_tsv_source_parity import _compare_maps, _make_absent_ids
    sys.path.pop(0)

    requested = [
        "chr1_100_A_G",   # present in both, AF/VEP match
        "chr1_200_C_T",   # present in both, VEP missing (empty dicts)
        "chr1_300_G_A",   # absent from both
        "chr2_500_T_C",   # present in VCF, missing in TSV
    ]
    absent_ids = _make_absent_ids([], 1, seed=0)
    # Replace the generic absent ID with the one in requested for the test.
    absent_ids = ["chr1_300_G_A"]

    af_vcf = {
        "chr1_100_A_G": 0.1,
        "chr1_200_C_T": 0.0,
        "chr2_500_T_C": 0.01,
    }
    vep_vcf = {
        "chr1_100_A_G": {"vep_IMPACT": "MODERATE"},
    }
    af_tsv = {
        "chr1_100_A_G": 0.1,
        "chr1_200_C_T": 0.0,
    }
    vep_tsv = {
        "chr1_100_A_G": {"vep_IMPACT": "MODERATE"},
    }

    report = _compare_maps(af_vcf, vep_vcf, af_tsv, vep_tsv, requested, absent_ids)

    assert report["requested"] == 4
    assert report["present_in_both"] == 2
    assert report["absent_from_both"] == 1
    assert report["absent_found_in_vcf"] == 0
    assert report["absent_found_in_tsv"] == 0
    assert report["missing_in_vcf"] == 0
    assert report["missing_in_tsv"] == 1
    assert report["af_mismatches"] == 0
    assert report["vep_mismatches"] == 0


def test_source_parity_make_absent_ids_are_non_canonical():
    """Generated absent IDs must use a chromosome not present in gnomAD."""
    import re
    import sys

    scripts_dir = str(Path(__file__).parent.parent / "scripts")
    sys.path.insert(0, scripts_dir)
    from validate_gnomad_tsv_source_parity import _make_absent_ids, CANONICAL_CHROM_RE
    sys.path.pop(0)

    absent = _make_absent_ids(["chr1_100_A_G"], 10, seed=0)
    assert len(absent) == 10
    for vid in absent:
        chrom = vid.split("_")[0]
        assert not CANONICAL_CHROM_RE.match(chrom)
