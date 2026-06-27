"""Shared utilities for variant naming, file I/O, and manifest handling."""

from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


def clean_name(s) -> str:
    """Normalize column names to match the training pipeline.

    Replaces colons, hyphens, and spaces with underscores, removes all
    non-alphanumeric/underscore characters, and collapses repeated
    underscores. This must stay identical to the cleaning used during
    model training.
    """
    s = str(s)
    s = s.replace(":", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9_]", "", s)
    s = re.sub(r"_+", "_", s)
    return s


def make_variant_id(chrom: str, pos, ref: str, alt: str) -> str:
    """Create the canonical variant_id used throughout the pipeline."""
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    pos = int(pos)
    ref = str(ref).upper()
    alt = str(alt).upper()
    return f"{chrom}_{pos}_{ref}_{alt}"


def parse_variant_id(variant_id: str) -> tuple[str, int, str, str]:
    """Parse a canonical variant_id into (chrom, pos, ref, alt)."""
    parts = variant_id.split("_")
    if len(parts) != 4:
        raise ValueError(f"Invalid variant_id: {variant_id}")
    return parts[0], int(parts[1]), parts[2], parts[3]


def normalize_allele(allele: str) -> str:
    """Uppercase and strip whitespace from an allele string."""
    return str(allele).strip().upper()


def valid_dna_alleles(*alleles: str) -> bool:
    """Return True if all alleles contain only ACGTN characters."""
    return all(set(a) <= set("ACGTN") for a in alleles if a)


def sha256_file(path: str | Path) -> str:
    """Compute the SHA256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def load_selected_features(path: str | Path) -> list[str]:
    """Load a newline-delimited feature list and clean each name."""
    with open(path, "r") as fh:
        return [clean_name(line.strip()) for line in fh if line.strip()]


def parse_vcf_to_dataframe(vcf_path: str | Path) -> pd.DataFrame:
    """Parse a VCF into a minimal DataFrame with variant/gene columns.

    Required columns in the output:
        - chrom, pos, ref, alt, gene

    The VCF is parsed with pysam. ``gene`` is extracted from (in order):
    ``INFO/GENE``, ``INFO/SYMBOL``, or the first VEP consequence's
    ``SYMBOL`` field. If no gene can be determined, the value is left
    empty; callers should handle missing genes appropriately.
    """
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "pysam is required for VCF parsing. Install it with: pip install pysam"
        ) from exc

    rows = []
    vcf_path = Path(vcf_path)
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            chrom = str(rec.chrom)
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            pos = rec.pos
            ref = str(rec.ref).upper()
            if not rec.alts:
                continue
            for alt in rec.alts:
                alt = str(alt).upper()
                gene = _extract_gene_from_record(rec)
                rows.append(
                    {
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": alt,
                        "gene": gene,
                    }
                )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    df["variant_id"] = df.apply(
        lambda r: make_variant_id(r["chrom"], r["pos"], r["ref"], r["alt"]),
        axis=1,
    )
    return df


def _extract_gene_from_record(rec) -> str:
    """Best-effort gene extraction from a pysam VariantRecord."""
    info = rec.info
    for key in ("GENE", "SYMBOL", "GENE_NAME"):
        if key in info:
            val = info[key]
            if isinstance(val, tuple):
                val = val[0]
            if val:
                return str(val)

    # Try VEP CSQ field
    if "CSQ" in info:
        csq = info["CSQ"]
        if isinstance(csq, tuple):
            csq = csq[0]
        parts = str(csq).split("|")
        if len(parts) > 3:
            return parts[3]

    if "vep" in info:
        vep = info["vep"]
        if isinstance(vep, tuple):
            vep = vep[0]
        parts = str(vep).split("|")
        if len(parts) > 3:
            return parts[3]

    return ""


def read_variant_list(path: str | Path) -> pd.DataFrame:
    """Read a tab-separated variant list (no header) into a DataFrame.

    Expected columns: chrom, pos, ref, alt, gene.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "gene"])
    df["chrom"] = df["chrom"].astype(str).apply(lambda c: c if c.startswith("chr") else f"chr{c}")
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()
    df["variant_id"] = df.apply(
        lambda r: make_variant_id(r["chrom"], r["pos"], r["ref"], r["alt"]),
        axis=1,
    )
    return df


def ensure_dir(path: str | Path) -> Path:
    """Create a directory if it does not already exist."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_gnomad_vcf_path(override: str | None = None) -> str:
    """Resolve the gnomAD VCF path from an override or the GNOMAD_VCF env var."""
    path = override or os.environ.get("GNOMAD_VCF")
    if not path:
        raise FileNotFoundError(
            "gnomAD VCF path not provided. Set the GNOMAD_VCF environment "
            "variable or pass --gnomad-vcf."
        )
    if not os.path.exists(path):
        raise FileNotFoundError(f"gnomAD VCF not found: {path}")
    return str(path)
