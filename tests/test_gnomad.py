from pathlib import Path

import pytest

from agvarpred_core.gnomad import query_gnomad_by_variants
from agvarpred_core.utils import get_gnomad_vcf_path, make_variant_id


def _gnomad_available() -> bool:
    try:
        get_gnomad_vcf_path()
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _gnomad_available(), reason="gnomAD VCF not available")
def test_query_gnomad_for_sample_variants():
    vcf_path = get_gnomad_vcf_path()
    # Use the first variant from examples/sample.vcf
    variant_ids = [make_variant_id("15", 67231868, "G", "A")]
    af_map, vep_map = query_gnomad_by_variants(variant_ids, vcf_path)
    # We only assert the function runs and returns dictionaries.
    assert isinstance(af_map, dict)
    assert isinstance(vep_map, dict)


def test_make_variant_id():
    assert make_variant_id("15", 67231868, "G", "A") == "chr15_67231868_G_A"
    assert make_variant_id("chr15", "67231868", "G", "A") == "chr15_67231868_G_A"
