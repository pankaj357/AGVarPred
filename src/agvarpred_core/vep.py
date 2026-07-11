"""VEP annotation parsing utilities."""

from __future__ import annotations

import re
from typing import Any


VEP_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND",
    "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID",
    "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL",
    "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET",
    "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
    "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF",
    "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF",
    "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
    "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter",
    "LoF_flags", "LoF_info",
]

VEP_IDX = {f: i for i, f in enumerate(VEP_FIELDS)}


def _extract_score_pred(raw: str | None) -> tuple[float | None, str | None]:
    """Parse strings like 'deleterious(0.01)' or 'probably_damaging(0.999)'."""
    if not raw:
        return None, None
    m = re.search(r"([^(]+)\(([\d.]+)\)", raw)
    if not m:
        return None, None
    try:
        score = float(m.group(2))
    except ValueError:
        score = None
    return score, m.group(1)


def parse_vep(vep_string: str | None) -> dict[str, Any]:
    """Parse a gnomAD/VEP CSQ annotation string.

    Uses the first transcript annotation only. Returns a dictionary of
    raw VEP-derived features suitable for further encoding by
    :func:`agvarpred_core.preprocessing.encode_vep_features`.
    """
    if not vep_string:
        return {}

    first_ann = vep_string.split(",")[0]
    parts = first_ann.split("|")
    if len(parts) < 65:
        return {}

    def get(idx: int) -> str | None:
        val = parts[idx] if idx < len(parts) else ""
        return val if val != "" else None

    sift_raw = get(VEP_IDX["SIFT"])
    sift_score, sift_pred = _extract_score_pred(sift_raw)

    polyphen_raw = get(VEP_IDX["PolyPhen"])
    polyphen_score, polyphen_pred = _extract_score_pred(polyphen_raw)

    consequence = get(VEP_IDX["Consequence"])
    impact = get(VEP_IDX["IMPACT"])
    lof = get(VEP_IDX["LoF"])
    protein_pos = get(VEP_IDX["Protein_position"])
    amino_acids = get(VEP_IDX["Amino_acids"])

    impact_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    impact_score = impact_map.get(impact, 0) if impact else 0

    protein_pos_num = None
    if protein_pos:
        m = re.search(r"^(\d+)", protein_pos)
        if m:
            try:
                protein_pos_num = int(m.group(1))
            except ValueError:
                pass

    cons_list = consequence.split("&") if consequence else []

    return {
        "vep_SIFT_score": sift_score,
        "vep_SIFT_pred": sift_pred,
        "vep_PolyPhen_score": polyphen_score,
        "vep_PolyPhen_pred": polyphen_pred,
        "vep_IMPACT": impact,
        "vep_IMPACT_score": impact_score,
        "vep_Consequence": consequence,
        "vep_is_missense": 1 if "missense_variant" in cons_list else 0,
        "vep_is_synonymous": 1 if "synonymous_variant" in cons_list else 0,
        "vep_is_stop_gained": 1 if "stop_gained" in cons_list else 0,
        "vep_is_frameshift": 1 if "frameshift_variant" in cons_list else 0,
        "vep_is_splice": 1 if any(
            c in cons_list for c in ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]
        ) else 0,
        "vep_LoF": lof,
        "vep_is_LoF_HC": 1 if lof == "HC" else 0,
        "vep_Protein_position": protein_pos_num,
        "vep_Amino_acids": amino_acids,
        "vep_has_SIFT": 1 if sift_score is not None else 0,
        "vep_has_PolyPhen": 1 if polyphen_score is not None else 0,
    }


def parse_vep_api(response: dict) -> dict[str, Any]:
    """Parse an Ensembl VEP REST API response for a single variant.

    ``response`` is the dictionary for one element of the VEP API response list.
    """
    consequence = response.get("most_severe_consequence", "")
    cons_list = consequence.split(",") if consequence else []

    transcripts = response.get("transcript_consequences", [])
    chosen = None
    for tx in transcripts:
        if tx.get("canonical") == 1:
            chosen = tx
            break
    if not chosen and transcripts:
        # Prefer a transcript with SIFT or PolyPhen annotations
        for tx in transcripts:
            if tx.get("sift_score") is not None or tx.get("polyphen_score") is not None:
                chosen = tx
                break
        if not chosen:
            chosen = transcripts[0]

    sift_score = chosen.get("sift_score") if chosen else None
    sift_pred = chosen.get("sift_prediction") if chosen else None
    polyphen_score = chosen.get("polyphen_score") if chosen else None
    polyphen_pred = chosen.get("polyphen_prediction") if chosen else None
    protein_pos = chosen.get("protein_start") if chosen else None

    impact = chosen.get("impact") if chosen else None
    impact_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    impact_score = impact_map.get(impact, 0) if impact else 0

    return {
        "vep_SIFT_score": sift_score,
        "vep_SIFT_pred": sift_pred,
        "vep_PolyPhen_score": polyphen_score,
        "vep_PolyPhen_pred": polyphen_pred,
        "vep_IMPACT": impact,
        "vep_IMPACT_score": impact_score,
        "vep_Consequence": consequence,
        "vep_is_missense": 1 if "missense_variant" in cons_list else 0,
        "vep_is_synonymous": 1 if "synonymous_variant" in cons_list else 0,
        "vep_is_stop_gained": 1 if "stop_gained" in cons_list else 0,
        "vep_is_frameshift": 1 if "frameshift_variant" in cons_list else 0,
        "vep_is_splice": 1 if any(
            c in cons_list for c in ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]
        ) else 0,
        "vep_LoF": None,
        "vep_is_LoF_HC": 0,
        "vep_Protein_position": protein_pos,
        "vep_Amino_acids": None,
        "vep_has_SIFT": 1 if sift_score is not None else 0,
        "vep_has_PolyPhen": 1 if polyphen_score is not None else 0,
    }
