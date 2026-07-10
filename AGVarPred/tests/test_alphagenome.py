import pandas as pd

from agvarpred_core.alphagenome import PrecomputedAlphaGenomeSource


def test_precomputed_alpha_source(tmp_path):
    base = tmp_path / "alpha"
    gene_dir = base / "GENE1_VEP"
    gene_dir.mkdir(parents=True)

    df = pd.DataFrame(
        {
            "ATAC__cell__track": [1.0, 2.0],
            "DNASE__cell__track": [3.0, 4.0],
        },
        index=pd.Index(["chr1_100_A_G", "chr1_200_C_T"], name="variant_id"),
    )
    df.to_parquet(gene_dir / "GENE1_ALL_VEP_RAW_SCORE_MATRIX.parquet")

    variants = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos": [100, 200],
            "ref": ["A", "C"],
            "alt": ["G", "T"],
            "gene": ["GENE1", "GENE1"],
            "variant_id": ["chr1_100_A_G", "chr1_200_C_T"],
        }
    )

    source = PrecomputedAlphaGenomeSource(base)
    out = source.get_features(variants)
    assert len(out) == 2
    # Column names are normalized by clean_name.
    assert "ATAC_cell_track" in out.columns
    assert "DNASE_cell_track" in out.columns
