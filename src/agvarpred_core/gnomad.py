"""gnomAD VCF querying for allele frequency and VEP annotations."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import parse_variant_id
from .vep import parse_vep


def _pysam_available() -> bool:
    try:
        import pysam  # noqa: F401
        return True
    except ImportError:
        return False


def _tabix_available() -> bool:
    return subprocess.run(["which", "tabix"], capture_output=True).returncode == 0


def _read_vcf_header(vcf_path: str | Path) -> list[str]:
    """Read the VCF header and return contig names that start with 'chr'."""
    try:
        import pysam
        with pysam.VariantFile(str(vcf_path)) as vcf:
            return [str(c) for c in vcf.header.contigs]
    except Exception:
        return []


def _normalize_chrom(chrom: str, header_chroms: list[str] | None = None) -> str:
    """Normalize chromosome name to match the gnomAD VCF."""
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    if header_chroms and chrom not in header_chroms:
        # Mitochondrial naming
        if chrom == "chrMT" and "chrM" in header_chroms:
            return "chrM"
        if chrom == "chrM" and "chrMT" in header_chroms:
            return "chrMT"
    return chrom


def query_gnomad_by_variants(
    variant_ids: list[str],
    vcf_path: str | Path,
    method: str = "auto",
) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
    """Query gnomAD for allele frequency and VEP annotations by variant ID.

    Parameters
    ----------
    variant_ids:
        List of canonical variant IDs (``chr_pos_ref_alt``).
    vcf_path:
        Path to the gnomAD VCF (bgzipped and tabixed).
    method:
        ``auto`` (default) chooses the fastest available method.
        Other options: ``tabix`` or ``pysam``.

    Returns
    -------
    af_map:
        Mapping from variant_id to gnomAD AF (0.0 if absent).
    vep_map:
        Mapping from variant_id to parsed VEP dictionary.
    """
    if method == "auto":
        if _tabix_available():
            method = "tabix"
        elif _pysam_available():
            method = "pysam"
        else:
            raise RuntimeError(
                "Neither tabix nor pysam is available. Install pysam or tabix to query gnomAD."
            )

    if method == "tabix":
        return _query_gnomad_tabix(variant_ids, vcf_path)
    if method == "pysam":
        return _query_gnomad_pysam(variant_ids, vcf_path)
    raise ValueError(f"Unknown method: {method}")


def _query_gnomad_pysam(
    variant_ids: list[str],
    vcf_path: str | Path,
) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
    import pysam

    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"gnomAD VCF not found: {vcf_path}")

    target = {vid: parse_variant_id(vid) for vid in variant_ids}
    # Index by (chrom, pos) for efficient fetching
    pos_index: dict[tuple[str, int], list[str]] = {}
    for vid, (chrom, pos, ref, alt) in target.items():
        pos_index.setdefault((chrom, pos), []).append(vid)

    header_chroms = _read_vcf_header(vcf_path)
    af_map: dict[str, float] = {}
    vep_map: dict[str, dict[str, Any]] = {}

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for (chrom, pos), vids in pos_index.items():
            norm_chrom = _normalize_chrom(chrom, header_chroms)
            try:
                for rec in vcf.fetch(norm_chrom, pos - 1, pos):
                    rec_chrom = str(rec.chrom)
                    if not rec_chrom.startswith("chr"):
                        rec_chrom = f"chr{rec_chrom}"
                    rec_pos = rec.pos
                    rec_ref = str(rec.ref).upper()
                    if not rec.alts:
                        continue
                    for i, alt in enumerate(rec.alts):
                        rec_alt = str(alt).upper()
                        for vid in vids:
                            _, _, want_ref, want_alt = target[vid]
                            if rec_ref == want_ref and rec_alt == want_alt:
                                af = rec.info.get("AF")
                                if af is None:
                                    af_value = 0.0
                                else:
                                    af_tuple = af if isinstance(af, tuple) else (af,)
                                    af_value = float(af_tuple[i]) if i < len(af_tuple) else 0.0
                                af_map[vid] = af_value

                                vep = rec.info.get("vep")
                                if vep:
                                    vep_raw = vep[0] if isinstance(vep, tuple) else vep
                                    vep_map[vid] = parse_vep(str(vep_raw))
            except Exception:
                # Position may not be in the VCF (e.g., alt contigs)
                continue

    return af_map, vep_map


def _query_gnomad_tabix(
    variant_ids: list[str],
    vcf_path: str | Path,
) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"gnomAD VCF not found: {vcf_path}")

    target = {vid: parse_variant_id(vid) for vid in variant_ids}
    bed_entries = []
    vid_lookup: dict[tuple[str, int, str, str], str] = {}
    for vid, (chrom, pos, ref, alt) in target.items():
        bed_entries.append((chrom, pos - 1, pos, vid))
        vid_lookup[(chrom, pos, ref, alt)] = vid

    bed_entries.sort(key=lambda x: (x[0], x[1]))

    tmp_dir = Path(os.environ.get("TMPDIR", "/tmp"))
    bed_path = tmp_dir / "agvarpred_gnomad_query.bed"
    with open(bed_path, "w") as fh:
        for chrom, start, end, _ in bed_entries:
            fh.write(f"{chrom}\t{start}\t{end}\n")

    cmd = ["tabix", "-R", str(bed_path), str(vcf_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"tabix failed: {proc.stderr}")

    af_map: dict[str, float] = {}
    vep_map: dict[str, dict[str, Any]] = {}

    for line in proc.stdout.strip().split("\n"):
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        chrom, pos, _, ref, alt = cols[0], int(cols[1]), cols[2], cols[3].upper(), cols[4]
        info = cols[7]
        alts = alt.split(",")
        for i, a in enumerate(alts):
            a = a.upper()
            key = (chrom, pos, ref, a)
            if key not in vid_lookup:
                continue
            vid = vid_lookup[key]
            m = re.search(r"AF=([^;]+)", info)
            if m:
                af_vals = m.group(1).split(",")
                af_map[vid] = float(af_vals[i]) if i < len(af_vals) else 0.0
            else:
                af_map[vid] = 0.0

            vep_m = re.search(r"vep=([^;]+)", info)
            if vep_m:
                vep_map[vid] = parse_vep(vep_m.group(1))

    return af_map, vep_map
