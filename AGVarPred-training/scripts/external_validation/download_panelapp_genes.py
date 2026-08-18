#!/usr/bin/env python3
"""
Download Genomics England PanelApp panels and identify expert-curated genes
NOT present in the training set.

PanelApp API docs: https://panelapp.genomicsengland.co.uk/api/docs/
"""

import json
import requests
import pandas as pd
from pathlib import Path
from time import sleep

ROOT = Path(__file__).resolve().parents[3]
TRAIN_GENES_PATH = ROOT / "external_validation/train_genes_upper.json"
OUTDIR = ROOT / "external_validation/panelapp_data"
OUTDIR.mkdir(exist_ok=True, parents=True)

PANELAPP_API = "https://panelapp.genomicsengland.co.uk/api/v1/panels"
DELAY = 0.2  # seconds between requests


def load_train_genes():
    with open(TRAIN_GENES_PATH) as f:
        return set(json.load(f))


def fetch_all_panels():
    """Fetch all PanelApp panels (paginated)."""
    print("🔄 Fetching PanelApp panels...")
    panels = []
    url = PANELAPP_API + "/"
    while url:
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            panels.extend(data.get("results", []))
            url = data.get("next")
            print(f"   Fetched {len(panels)} panels so far...")
            sleep(DELAY)
        except Exception as e:
            print(f"   ⚠ Error: {e}")
            break
    print(f"✅ Total panels: {len(panels)}")
    return panels


def fetch_panel_genes(panel_id):
    """Fetch genes for a specific panel."""
    url = f"{PANELAPP_API}/{panel_id}/"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        genes = []
        for g in data.get("genes", []):
            gene_data = g.get("gene_data", {})
            symbol = gene_data.get("gene_symbol", "")
            if symbol:
                genes.append({
                    "panel_id": panel_id,
                    "panel_name": data.get("name", ""),
                    "gene_symbol": symbol.upper(),
                    "entity_type": g.get("entity_type", ""),
                    "confidence_level": g.get("confidence_level", ""),
                    "mode_of_inheritance": g.get("mode_of_inheritance", ""),
                })
        return genes
    except Exception as e:
        print(f"   ⚠ Panel {panel_id} error: {e}")
        return []


def main():
    train_genes = load_train_genes()
    print(f"Training genes: {len(train_genes):,}")

    panels = fetch_all_panels()
    if not panels:
        print("❌ No panels fetched. Exiting.")
        return

    # Save raw panel list
    pd.DataFrame(panels).to_csv(OUTDIR / "panels_list.csv", index=False)

    print("\n🔄 Fetching genes per panel (this may take a few minutes)...")
    all_genes = []
    for i, panel in enumerate(panels):
        panel_id = panel.get("id")
        if not panel_id:
            continue
        genes = fetch_panel_genes(panel_id)
        all_genes.extend(genes)
        if (i + 1) % 50 == 0:
            print(f"   Processed {i+1}/{len(panels)} panels | Genes collected: {len(all_genes)}")
        sleep(DELAY)

    df = pd.DataFrame(all_genes)
    df.to_csv(OUTDIR / "panelapp_all_genes.csv", index=False)
    print(f"✅ Saved: {OUTDIR / 'panelapp_all_genes.csv'}")

    # Identify novel genes
    df["in_training"] = df["gene_symbol"].isin(train_genes)
    novel = df[~df["in_training"]].copy()
    novel_genes = set(novel["gene_symbol"].unique())

    print(f"\n📊 PANELAPP GENE ANALYSIS")
    print(f"   Total gene-panel entries: {len(df):,}")
    print(f"   Unique genes: {df['gene_symbol'].nunique():,}")
    print(f"   Genes in training: {df['in_training'].sum():,}")
    print(f"   Genes NOT in training: {len(novel):,}")
    print(f"   Unique novel genes: {len(novel_genes):,}")

    novel.to_csv(OUTDIR / "panelapp_novel_genes.csv", index=False)
    print(f"✅ Saved: {OUTDIR / 'panelapp_novel_genes.csv'}")

    # Summary by panel
    panel_summary = []
    for panel_name, group in novel.groupby("panel_name"):
        panel_summary.append({
            "panel": panel_name,
            "novel_genes": group["gene_symbol"].nunique(),
            "entries": len(group),
        })
    summary_df = pd.DataFrame(panel_summary).sort_values("novel_genes", ascending=False)
    summary_df.to_csv(OUTDIR / "panelapp_novel_panels_summary.csv", index=False)
    print(f"✅ Saved: {OUTDIR / 'panelapp_novel_panels_summary.csv'}")

    print(f"\n📝 NEXT STEPS:")
    print(f"   1. Review {OUTDIR / 'panelapp_novel_panels_summary.csv'} for top panels")
    print(f"   2. Apply for Genomics England 100KGP Research Environment")
    print(f"      URL: https://www.genomicsengland.co.uk/research/academic-research")
    print(f"   3. Query 100KGP for pathogenic variants in the novel genes")


if __name__ == "__main__":
    main()
