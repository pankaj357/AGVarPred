#!/usr/bin/env python3
"""
Prepare DECIPHER bulk export request specifications.

Run this script to generate:
  1. A ready-to-send email template
  2. Gene statistics for your request

Then copy the email text and send to: decipher@sanger.ac.uk
"""

import json
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[3]
TRAIN_GENES_PATH = ROOT / "external_validation/train_genes_upper.json"

def main():
    with open(TRAIN_GENES_PATH) as f:
        train_genes = set(json.load(f))

    # DECIPHER covers primarily developmental disorder (DD) and congenital anomaly genes
    # We want genes NOT in training to maintain independence
    email_template = f"""Subject: Bulk Variant Export Request for Academic Research — Pathogenic SNVs (GRCh38)

Dear DECIPHER Team,

I am writing to request a bulk export of variant data from DECIPHER for use in external validation of a machine learning model for variant pathogenicity prediction.

PROJECT CONTEXT
---------------
- Institution: [YOUR INSTITUTION]
- Principal Investigator: [YOUR NAME]
- Purpose: Independent benchmark for a germline variant pathogenicity classifier
- The model was trained on ClinVar-derived data; we need fully independent expert-curated variants for validation

REQUESTED DATA
--------------
1. Variant type: Single nucleotide variants (SNVs) and small indels only
2. Assembly: GRCh38 (chromosome, position, ref, alt)
3. Gene constraint: Only variants in genes NOT present in our training set ({len(train_genes):,} genes — list attached as train_genes.txt)
4. Classification: Pathogenic (P) and Likely Pathogenic (LP) only
5. Minimum review: Variants with at least one supporting publication or clinical assertion
6. Required columns:
   - chrom, pos, ref, alt
   - gene_symbol
   - decipher_patient_id (or aggregate count if patient-level restricted)
   - pathogenicity assertion
   - inheritance_pattern (if available)
   - phenotype_HPO_terms (if available)

INDEPENDENCE VERIFICATION
-------------------------
We will apply strict filtering:
- Zero variant overlap with ClinVar training data (134,002 variants)
- Zero gene overlap with training gene set ({len(train_genes):,} genes)
- Results will be used solely for academic publication and model benchmarking

DATA HANDLING
-------------
- Data will be stored on a secure institutional server
- No re-distribution or commercial use
- DECIPHER will be cited as the data source in all publications
- We are happy to sign a data sharing agreement if required

Please let me know if you need any additional information or documentation.

Thank you for maintaining this invaluable resource.

Best regards,
[YOUR NAME]
[YOUR EMAIL]
[YOUR INSTITUTION]
[DATE]
"""

    out_path = ROOT / "external_validation/DECIPHER_export_request.txt"
    with open(out_path, "w") as f:
        f.write(email_template)

    # Also save the training gene list as a simple text file for attachment
    genes_path = ROOT / "external_validation/train_genes_for_attachment.txt"
    with open(genes_path, "w") as f:
        for g in sorted(train_genes):
            f.write(g + "\n")

    print("=" * 70)
    print("DECIPHER EXPORT REQUEST PREPARED")
    print("=" * 70)
    print(f"\n📧 Email template saved to:")
    print(f"   {out_path}")
    print(f"\n📎 Training gene list attachment saved to:")
    print(f"   {genes_path} ({len(train_genes):,} genes)")
    print(f"\n📝 NEXT STEPS:")
    print(f"   1. Edit the template with your name/institution")
    print(f"   2. Send to: decipher@sanger.ac.uk")
    print(f"   3. Attach: train_genes_for_attachment.txt")
    print(f"   4. Typical response time: 1–2 weeks")
    print("=" * 70)


if __name__ == "__main__":
    main()
