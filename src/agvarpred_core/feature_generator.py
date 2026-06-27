"""FeatureGenerator: VCF -> full engineered feature matrix."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import pandas as pd

from .af_source import AFSource, NoAFSource, resolve_af_source
from .alphagenome import resolve_alphagenome_source
from .preprocessing import add_af_features, clean_dataframe_columns, encode_vep_features
from .utils import parse_vcf_to_dataframe, valid_dna_alleles


class FeatureGenerator:
    """Generate the complete engineered feature matrix for a set of variants.

    The generator is model-independent: it always returns the full set of
    AlphaGenome functional scores plus VEP and gnomAD-derived features. A
    downstream :class:`agvarpred_core.feature_selector.FeatureSelector` can
    then choose the subset required by a specific model.

    Parameters
    ----------
    af_source:
        An :class:`AFSource` instance for AF/VEP annotations. If None, an
        AF source is resolved automatically (local gnomAD → online → no-AF).
    alpha_mode:
        ``auto`` (default), ``sdk``, or ``precomputed``.
    alpha_api_key:
        AlphaGenome API key. Required for ``sdk`` mode; falls back to
        ``ALPHAGENOME_API_KEY`` environment variable.
    alpha_dir:
        Directory of precomputed AlphaGenome feature matrices. Required for
        ``precomputed`` mode and optional for ``auto``.
    """

    def __init__(
        self,
        af_source: AFSource | None = None,
        alpha_mode: Literal["auto", "sdk", "precomputed"] = "auto",
        alpha_api_key: str | None = None,
        alpha_dir: str | Path | None = None,
    ):
        self.af_source = af_source or resolve_af_source()
        self.alpha_mode = alpha_mode
        self.alpha_api_key = alpha_api_key
        self.alpha_dir = alpha_dir

    def from_vcf(self, vcf_path: str | Path) -> pd.DataFrame:
        """Generate features from a VCF file."""
        variants_df = parse_vcf_to_dataframe(vcf_path)
        return self.from_variants(variants_df)

    def from_variants(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Generate features from a DataFrame of variants.

        Required columns: chrom, pos, ref, alt, gene. A ``variant_id`` column
        will be created if not present.
        """
        if variants_df.empty:
            return pd.DataFrame(index=pd.Index([], name="variant_id"))

        variants_df = variants_df.copy()
        variants_df["ref"] = variants_df["ref"].astype(str).str.upper()
        variants_df["alt"] = variants_df["alt"].astype(str).str.upper()
        variants_df = variants_df[
            variants_df.apply(
                lambda r: valid_dna_alleles(r["ref"], r["alt"]), axis=1
            )
        ]
        if variants_df.empty:
            return pd.DataFrame(index=pd.Index([], name="variant_id"))

        # AlphaGenome functional scores
        alpha_source = resolve_alphagenome_source(
            self.alpha_mode,
            api_key=self.alpha_api_key,
            precomputed_dir=self.alpha_dir,
        )
        alpha_df = alpha_source.get_features(variants_df)

        # AF + VEP annotations (provider may be no-AF)
        variant_ids = variants_df["variant_id"].tolist()
        af_map, vep_map = self.af_source.query(variant_ids)

        # Build annotation DataFrame
        annotation_rows = []
        for vid in variant_ids:
            row = {"variant_id": vid}
            row["gnomAD_AF"] = af_map.get(vid, 0.0)
            vep = vep_map.get(vid, {})
            row.update(vep)
            annotation_rows.append(row)
        annotation_df = pd.DataFrame(annotation_rows)
        if not annotation_df.empty:
            annotation_df = annotation_df.set_index("variant_id")

        # Merge AlphaGenome + annotations
        # If AlphaGenome matrices happen to contain AF/VEP columns, prefer the
        # values from the annotation provider.
        if alpha_df.empty:
            combined = annotation_df
        else:
            overlap = alpha_df.columns.intersection(annotation_df.columns)
            if not overlap.empty:
                alpha_df = alpha_df.drop(columns=overlap)
            combined = alpha_df.join(annotation_df, how="left")

        # Encode VEP and AF features
        combined = encode_vep_features(combined)
        combined = add_af_features(combined)

        # Clean column names (again, to be safe)
        combined = clean_dataframe_columns(combined)

        # Preserve original variant order
        combined = combined.reindex(variants_df["variant_id"].unique())
        combined.index.name = "variant_id"
        return combined

    def __call__(self, vcf_path: str | Path) -> pd.DataFrame:
        """Convenience alias for ``from_vcf``."""
        return self.from_vcf(vcf_path)
