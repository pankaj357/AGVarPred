"""AlphaGenome feature retrieval: SDK or precomputed matrices."""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Literal

import pandas as pd

from .utils import clean_name


class AlphaGenomeSource:
    """Base class for AlphaGenome feature sources."""

    def get_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Return a DataFrame indexed by variant_id with AlphaGenome score columns."""
        raise NotImplementedError


class AlphaGenomeSDKSource(AlphaGenomeSource):
    """Query the AlphaGenome Python SDK for functional genomic scores."""

    def __init__(self, api_key: str | None = None):
        try:
            from alphagenome.data import genome
            from alphagenome.models import dna_client, variant_scorers
        except ImportError as exc:
            raise ImportError(
                "AlphaGenome SDK not installed. Install it or use precomputed features."
            ) from exc

        self._genome = genome
        self._dna_client = dna_client
        self._variant_scorers = variant_scorers
        self.api_key = api_key or os.environ.get("ALPHAGENOME_API_KEY")
        if not self.api_key:
            raise ValueError(
                "AlphaGenome API key required. Set ALPHAGENOME_API_KEY or pass api_key."
            )
        self.dna_model = dna_client.create(self.api_key)
        self.all_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
        self.seq_len = dna_client.SEQUENCE_LENGTH_1MB

    def get_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        rows = []
        all_columns = set()
        for _, row in variants_df.iterrows():
            chrom = str(row["chrom"])
            if chrom == "chrMT":
                chrom = "chrM"
            pos = int(row["pos"])
            ref = str(row["ref"])
            alt = str(row["alt"])
            variant_id = str(row["variant_id"])

            try:
                variant_obj = self._genome.Variant(
                    chromosome=chrom,
                    position=pos,
                    reference_bases=ref,
                    alternate_bases=alt,
                    name=variant_id,
                )
                interval = variant_obj.reference_interval.resize(self.seq_len)
                scores = self.dna_model.score_variant(
                    interval=interval,
                    variant=variant_obj,
                    variant_scorers=self.all_scorers,
                    organism=self._dna_client.Organism.HOMO_SAPIENS,
                )
                df_scores = self._variant_scorers.tidy_scores(scores)
                vector = {"variant_id": variant_id}
                for _, s in df_scores.iterrows():
                    col = f"{s.get('output_type', 'unk')}__{s.get('biosample_name', 'unk')}__{s.get('track_name', 'unk')}"
                    vector[col] = s.get("raw_score", None)
                    all_columns.add(col)
                rows.append(vector)
            except Exception as exc:
                warnings.warn(f"AlphaGenome failed for {variant_id}: {exc}")
                continue

        if not rows:
            return pd.DataFrame(index=pd.Index([], name="variant_id"))

        df = pd.DataFrame(rows).set_index("variant_id")
        df = df.reindex(columns=sorted(all_columns))
        df.columns = [clean_name(c) for c in df.columns]
        return df


class PrecomputedAlphaGenomeSource(AlphaGenomeSource):
    """Load AlphaGenome scores from a directory of precomputed matrices.

    Expected layout:
        <base_dir>/
            <gene>_VEP/
                <gene>_ALL_VEP_RAW_SCORE_MATRIX.parquet
    """

    def __init__(self, base_dir: str | Path):
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"Precomputed AlphaGenome directory not found: {self.base_dir}")

    def get_features(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        if variants_df.empty:
            return pd.DataFrame(index=pd.Index([], name="variant_id"))

        required_genes = variants_df["gene"].dropna().unique()
        dfs = []
        for gene in required_genes:
            gene_dir = self.base_dir / f"{gene}_VEP"
            if not gene_dir.exists():
                continue
            parquet_files = list(gene_dir.glob("*.parquet"))
            if not parquet_files:
                continue
            for f in parquet_files:
                try:
                    df = pd.read_parquet(f)
                    df.columns = [clean_name(c) for c in df.columns]
                    if "variant_id" in df.columns:
                        df = df.set_index("variant_id")
                    elif df.index.name is None:
                        df.index.name = "variant_id"
                    dfs.append(df)
                except Exception as exc:
                    warnings.warn(f"Failed to read {f}: {exc}")

        if not dfs:
            return pd.DataFrame(index=pd.Index([], name="variant_id"))

        combined = pd.concat(dfs, ignore_index=False)
        combined = combined[~combined.index.duplicated(keep="first")]
        # Keep only requested variants
        wanted = set(variants_df["variant_id"])
        combined = combined[combined.index.isin(wanted)]
        return combined


def resolve_alphagenome_source(
    mode: Literal["auto", "sdk", "precomputed"],
    api_key: str | None = None,
    precomputed_dir: str | Path | None = None,
) -> AlphaGenomeSource:
    """Return an AlphaGenome feature source according to the chosen mode.

    Priority in ``auto`` mode:
        1. AlphaGenome SDK + API key if available.
        2. Precomputed matrices if ``precomputed_dir`` is provided.
        3. Clear error if neither is available.
    """
    if mode == "sdk":
        return AlphaGenomeSDKSource(api_key=api_key)

    if mode == "precomputed":
        if not precomputed_dir:
            raise ValueError("precomputed_dir is required when alpha-mode=precomputed")
        return PrecomputedAlphaGenomeSource(precomputed_dir)

    if mode == "auto":
        sdk_available = _sdk_available()
        if sdk_available and (api_key or os.environ.get("ALPHAGENOME_API_KEY")):
            return AlphaGenomeSDKSource(api_key=api_key)
        if precomputed_dir:
            return PrecomputedAlphaGenomeSource(precomputed_dir)
        raise RuntimeError(
            "AlphaGenome features cannot be generated automatically. "
            "Either install the AlphaGenome SDK and set ALPHAGENOME_API_KEY, "
            "or provide a precomputed AlphaGenome feature directory with "
            "--alpha-dir / alpha_mode=precomputed."
        )

    raise ValueError(f"Unknown alpha mode: {mode}")


def _sdk_available() -> bool:
    try:
        import alphagenome  # noqa: F401
        return True
    except ImportError:
        return False
