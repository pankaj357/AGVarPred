#!/usr/bin/env python3
"""
Create a publication-quality multi-panel benchmark comparison figure.
Panel 1: Coverage comparison (% of variants scored by each tool)
Panel 2: AUC on full benchmark (our model) vs. same-subset (all tools)
Panel 3: Coverage-adjusted AUC (only variants where ALL tools have scores)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score

# Use a publication-ready style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
})

ROOT = Path(os.environ.get("PROJECT_ROOT", "."))
RESULTS_ROOT = ROOT / "external_validation/results"
TOOL_DIR = RESULTS_ROOT / "benchmark_tools"
OUTDIR = ROOT / "external_validation/results/benchmark_tools"
OUTDIR.mkdir(parents=True, exist_ok=True)

BENCHMARKS = [
    ("humsavar", "Humsavar"),
    ("mave_independent", "MAVE Independent"),
    ("gnomad_benign", "gnomAD Benign"),
    ("vip", "VIP"),
]

# Colors
COLOR_OUR = '#2E86AB'      # Blue
COLOR_CADD = '#A23B72'     # Magenta
COLOR_REVEL = '#F18F01'    # Orange
COLOR_AM = '#C73E1D'       # Red


def load_our_predictions(bench_name):
    """Load our model predictions."""
    pred_file = RESULTS_ROOT / bench_name / "regularized" / f"{bench_name}_predictions.csv"
    if not pred_file.exists():
        # Try alternative naming
        alt_names = {
            "mave_independent": "mave_independent_predictions.csv",
        }
        if bench_name in alt_names:
            pred_file = RESULTS_ROOT / bench_name / "regularized" / alt_names[bench_name]
    
    if pred_file.exists():
        df = pd.read_csv(pred_file)
        # Standardize column names
        if 'prob_calibrated' in df.columns:
            df['score'] = df['prob_calibrated']
        elif 'pathogenic_probability' in df.columns:
            df['score'] = df['pathogenic_probability']
        elif 'probability' in df.columns:
            df['score'] = df['probability']
        return df[['variant_id', 'score', 'true_label']].dropna()
    return None


def load_tool_scores(bench_name):
    """Load CADD/REVEL/AlphaMissense scores."""
    tool_file = TOOL_DIR / f"{bench_name}_tool_scores.csv"
    if tool_file.exists():
        df = pd.read_csv(tool_file)
        # Create variant_id to match predictions (add chr prefix if missing)
        chrom = df['chrom'].astype(str)
        chrom = chrom.apply(lambda x: x if x.startswith('chr') else f'chr{x}')
        df['variant_id'] = chrom + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']
        return df
    return None


def compute_auc(y_true, y_scores):
    """Compute ROC-AUC safely."""
    mask = ~np.isnan(y_scores)
    y_t = y_true[mask]
    y_s = y_scores[mask]
    if len(y_t) == 0 or len(np.unique(y_t)) < 2:
        return np.nan
    return roc_auc_score(y_t, y_s)


def analyze_benchmark(bench_key, bench_label):
    """Analyze a single benchmark."""
    our_df = load_our_predictions(bench_key)
    tool_df = load_tool_scores(bench_key)
    
    if our_df is None or tool_df is None:
        return None
    
    n_total = len(tool_df)
    
    # Coverage
    cov_our = len(our_df) / n_total * 100 if n_total > 0 else 0
    cov_cadd = tool_df['cadd'].notna().sum() / n_total * 100
    cov_revel = tool_df['revel'].notna().sum() / n_total * 100
    cov_am = tool_df['alphamissense'].notna().sum() / n_total * 100
    
    # Merge our predictions with tool scores on variant_id
    merged = tool_df.copy()
    our_scores = our_df.set_index('variant_id')['score']
    merged['our_score'] = merged['variant_id'].map(our_scores)
    
    # AUC on FULL benchmark (our model)
    auc_our_full = compute_auc(merged['label'].values, merged['our_score'].values)
    
    # AUC on same subset where each tool has scores
    auc_our_on_cadd = compute_auc(
        merged.loc[merged['cadd'].notna(), 'label'].values,
        merged.loc[merged['cadd'].notna(), 'our_score'].values
    )
    auc_cadd = compute_auc(merged['label'].values, merged['cadd'].values)
    
    auc_our_on_revel = compute_auc(
        merged.loc[merged['revel'].notna(), 'label'].values,
        merged.loc[merged['revel'].notna(), 'our_score'].values
    )
    auc_revel = compute_auc(merged['label'].values, merged['revel'].values)
    
    auc_our_on_am = compute_auc(
        merged.loc[merged['alphamissense'].notna(), 'label'].values,
        merged.loc[merged['alphamissense'].notna(), 'our_score'].values
    )
    auc_am = compute_auc(merged['label'].values, merged['alphamissense'].values)
    
    # Coverage-adjusted: only variants where ALL tools have scores
    all_mask = (
        merged['our_score'].notna() &
        merged['cadd'].notna() &
        merged['revel'].notna() &
        merged['alphamissense'].notna()
    )
    
    auc_our_adjusted = compute_auc(
        merged.loc[all_mask, 'label'].values,
        merged.loc[all_mask, 'our_score'].values
    )
    auc_cadd_adjusted = compute_auc(
        merged.loc[all_mask, 'label'].values,
        merged.loc[all_mask, 'cadd'].values
    )
    auc_revel_adjusted = compute_auc(
        merged.loc[all_mask, 'label'].values,
        merged.loc[all_mask, 'revel'].values
    )
    auc_am_adjusted = compute_auc(
        merged.loc[all_mask, 'label'].values,
        merged.loc[all_mask, 'alphamissense'].values
    )
    n_adjusted = all_mask.sum()
    
    return {
        'benchmark': bench_label,
        'n_total': n_total,
        'coverage': {
            'our': cov_our,
            'cadd': cov_cadd,
            'revel': cov_revel,
            'am': cov_am,
        },
        'auc_full': {
            'our': auc_our_full,
        },
        'auc_same_subset': {
            'our_on_cadd': auc_our_on_cadd,
            'cadd': auc_cadd,
            'our_on_revel': auc_our_on_revel,
            'revel': auc_revel,
            'our_on_am': auc_our_on_am,
            'am': auc_am,
        },
        'auc_adjusted': {
            'n': n_adjusted,
            'our': auc_our_adjusted,
            'cadd': auc_cadd_adjusted,
            'revel': auc_revel_adjusted,
            'am': auc_am_adjusted,
        },
    }


def create_figure(results):
    """Create publication-quality multi-panel figure."""
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
    
    benchmarks = [r['benchmark'] for r in results]
    x = np.arange(len(benchmarks))
    width = 0.18
    
    # === Panel A: Coverage Comparison ===
    ax1 = fig.add_subplot(gs[0, 0])
    
    cov_our = [r['coverage']['our'] for r in results]
    cov_cadd = [r['coverage']['cadd'] for r in results]
    cov_revel = [r['coverage']['revel'] for r in results]
    cov_am = [r['coverage']['am'] for r in results]
    
    ax1.bar(x - 1.5*width, cov_our, width, label='Our Model', color=COLOR_OUR, edgecolor='black', linewidth=0.5)
    ax1.bar(x - 0.5*width, cov_cadd, width, label='CADD', color=COLOR_CADD, edgecolor='black', linewidth=0.5)
    ax1.bar(x + 0.5*width, cov_revel, width, label='REVEL', color=COLOR_REVEL, edgecolor='black', linewidth=0.5)
    ax1.bar(x + 1.5*width, cov_am, width, label='AlphaMissense', color=COLOR_AM, edgecolor='black', linewidth=0.5)
    
    ax1.set_ylabel('Coverage (%)')
    ax1.set_title('A. Benchmark Coverage', fontweight='bold', loc='left')
    ax1.set_xticks(x)
    ax1.set_xticklabels(benchmarks, rotation=15, ha='right')
    ax1.set_ylim(0, 110)
    ax1.axhline(y=100, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax1.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Add value labels on bars
    for i, (v1, v2, v3, v4) in enumerate(zip(cov_our, cov_cadd, cov_revel, cov_am)):
        if v1 < 100:
            ax1.text(i - 1.5*width, v1 + 2, f'{v1:.1f}%', ha='center', va='bottom', fontsize=7)
        if v2 > 0:
            ax1.text(i - 0.5*width, v2 + 2, f'{v2:.1f}%', ha='center', va='bottom', fontsize=7)
        if v3 > 0:
            ax1.text(i + 0.5*width, v3 + 2, f'{v3:.1f}%', ha='center', va='bottom', fontsize=7)
        if v4 > 0:
            ax1.text(i + 1.5*width, v4 + 2, f'{v4:.1f}%', ha='center', va='bottom', fontsize=7)
    
    # === Panel B: AUC on Tool's Subset ===
    ax2 = fig.add_subplot(gs[0, 1])
    
    # For each tool, show our model's AUC on that tool's subset vs. the tool's AUC
    tool_names = ['CADD\nSubset', 'REVEL\nSubset', 'AlphaMissense\nSubset']
    
    for i, bench in enumerate(benchmarks):
        r = results[i]
        subset_aucs = [
            (r['auc_same_subset']['our_on_cadd'], r['auc_same_subset']['cadd']),
            (r['auc_same_subset']['our_on_revel'], r['auc_same_subset']['revel']),
            (r['auc_same_subset']['our_on_am'], r['auc_same_subset']['am']),
        ]
        
        x_pos = np.array([0, 1, 2]) + i * 4
        for j, (our_auc, tool_auc) in enumerate(subset_aucs):
            if not np.isnan(our_auc) and not np.isnan(tool_auc):
                ax2.bar(x_pos[j] - width/2, our_auc, width, color=COLOR_OUR, edgecolor='black', linewidth=0.5, alpha=0.9)
                ax2.bar(x_pos[j] + width/2, tool_auc, width, color=[COLOR_CADD, COLOR_REVEL, COLOR_AM][j], edgecolor='black', linewidth=0.5, alpha=0.9)
    
    ax2.set_ylabel('ROC-AUC')
    ax2.set_title('B. AUC on Tool\'s Scored Subset', fontweight='bold', loc='left')
    ax2.set_xticks([0.5, 1.5, 2.5, 4.5, 5.5, 6.5, 8.5, 9.5, 10.5])
    ax2.set_xticklabels(['CADD', 'REVEL', 'αMiss', 'CADD', 'REVEL', 'αMiss', 'CADD', 'REVEL', 'αMiss'], fontsize=8)
    ax2.set_ylim(0, 1.05)
    ax2.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Add benchmark labels
    for i, bench in enumerate(benchmarks):
        ax2.text(1.5 + i * 4, -0.08, bench, ha='center', va='top', fontsize=9, fontweight='bold', transform=ax2.get_xaxis_transform())
    
    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLOR_OUR, edgecolor='black', label='Our Model'),
        Patch(facecolor=COLOR_CADD, edgecolor='black', label='CADD'),
        Patch(facecolor=COLOR_REVEL, edgecolor='black', label='REVEL'),
        Patch(facecolor=COLOR_AM, edgecolor='black', label='AlphaMissense'),
    ]
    ax2.legend(handles=legend_elements, loc='upper right', frameon=True, fancybox=False, edgecolor='black')
    
    # === Panel C: Coverage-Adjusted AUC ===
    ax3 = fig.add_subplot(gs[1, :])
    
    adj_our = [r['auc_adjusted']['our'] for r in results]
    adj_cadd = [r['auc_adjusted']['cadd'] for r in results]
    adj_revel = [r['auc_adjusted']['revel'] for r in results]
    adj_am = [r['auc_adjusted']['am'] for r in results]
    adj_n = [r['auc_adjusted']['n'] for r in results]
    
    x = np.arange(len(benchmarks))
    ax3.bar(x - 1.5*width, adj_our, width, label='Our Model', color=COLOR_OUR, edgecolor='black', linewidth=0.5)
    ax3.bar(x - 0.5*width, adj_cadd, width, label='CADD', color=COLOR_CADD, edgecolor='black', linewidth=0.5)
    ax3.bar(x + 0.5*width, adj_revel, width, label='REVEL', color=COLOR_REVEL, edgecolor='black', linewidth=0.5)
    ax3.bar(x + 1.5*width, adj_am, width, label='AlphaMissense', color=COLOR_AM, edgecolor='black', linewidth=0.5)
    
    ax3.set_ylabel('ROC-AUC')
    ax3.set_title('C. Coverage-Adjusted AUC (Variants Scored by ALL Tools)', fontweight='bold', loc='left')
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{b}\n(n={n})" for b, n in zip(benchmarks, adj_n)])
    ax3.set_ylim(0, 1.05)
    ax3.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax3.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    # Add value labels
    for i, (v1, v2, v3, v4) in enumerate(zip(adj_our, adj_cadd, adj_revel, adj_am)):
        if not np.isnan(v1):
            ax3.text(i - 1.5*width, v1 + 0.02, f'{v1:.3f}', ha='center', va='bottom', fontsize=8)
        if not np.isnan(v2):
            ax3.text(i - 0.5*width, v2 + 0.02, f'{v2:.3f}', ha='center', va='bottom', fontsize=8)
        if not np.isnan(v3):
            ax3.text(i + 0.5*width, v3 + 0.02, f'{v3:.3f}', ha='center', va='bottom', fontsize=8)
        if not np.isnan(v4):
            ax3.text(i + 1.5*width, v4 + 0.02, f'{v4:.3f}', ha='center', va='bottom', fontsize=8)
    
    plt.suptitle('External Validation Benchmark Comparison', fontsize=14, fontweight='bold', y=0.98)
    
    out_path = OUTDIR / 'benchmark_comparison_figure.png'
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.savefig(out_path.with_suffix('.pdf'), bbox_inches='tight', dpi=300)
    print(f"Figure saved to {out_path}")
    
    return fig


def create_summary_table(results):
    """Create a summary table of all metrics."""
    rows = []
    for r in results:
        rows.append({
            'Benchmark': r['benchmark'],
            'N Total': r['n_total'],
            'Our Coverage %': f"{r['coverage']['our']:.1f}",
            'CADD Coverage %': f"{r['coverage']['cadd']:.1f}",
            'REVEL Coverage %': f"{r['coverage']['revel']:.1f}",
            'AM Coverage %': f"{r['coverage']['am']:.1f}",
            'Our AUC (Full)': f"{r['auc_full']['our']:.4f}" if not np.isnan(r['auc_full']['our']) else 'N/A',
            'Our AUC (CADD Subset)': f"{r['auc_same_subset']['our_on_cadd']:.4f}" if not np.isnan(r['auc_same_subset']['our_on_cadd']) else 'N/A',
            'CADD AUC': f"{r['auc_same_subset']['cadd']:.4f}" if not np.isnan(r['auc_same_subset']['cadd']) else 'N/A',
            'Our AUC (REVEL Subset)': f"{r['auc_same_subset']['our_on_revel']:.4f}" if not np.isnan(r['auc_same_subset']['our_on_revel']) else 'N/A',
            'REVEL AUC': f"{r['auc_same_subset']['revel']:.4f}" if not np.isnan(r['auc_same_subset']['revel']) else 'N/A',
            'Our AUC (AM Subset)': f"{r['auc_same_subset']['our_on_am']:.4f}" if not np.isnan(r['auc_same_subset']['our_on_am']) else 'N/A',
            'AM AUC': f"{r['auc_same_subset']['am']:.4f}" if not np.isnan(r['auc_same_subset']['am']) else 'N/A',
            'N (All Tools)': r['auc_adjusted']['n'],
            'Our AUC (Adjusted)': f"{r['auc_adjusted']['our']:.4f}" if not np.isnan(r['auc_adjusted']['our']) else 'N/A',
            'CADD AUC (Adjusted)': f"{r['auc_adjusted']['cadd']:.4f}" if not np.isnan(r['auc_adjusted']['cadd']) else 'N/A',
            'REVEL AUC (Adjusted)': f"{r['auc_adjusted']['revel']:.4f}" if not np.isnan(r['auc_adjusted']['revel']) else 'N/A',
            'AM AUC (Adjusted)': f"{r['auc_adjusted']['am']:.4f}" if not np.isnan(r['auc_adjusted']['am']) else 'N/A',
        })
    
    df = pd.DataFrame(rows)
    csv_path = OUTDIR / 'coverage_adjusted_comparison.csv'
    df.to_csv(csv_path, index=False)
    print(f"Summary table saved to {csv_path}")
    return df


def main():
    print("Analyzing benchmarks...")
    results = []
    for bench_key, bench_label in BENCHMARKS:
        print(f"\n{bench_label}:")
        r = analyze_benchmark(bench_key, bench_label)
        if r:
            results.append(r)
            print(f"  Coverage: Our={r['coverage']['our']:.1f}%, CADD={r['coverage']['cadd']:.1f}%, REVEL={r['coverage']['revel']:.1f}%, AM={r['coverage']['am']:.1f}%")
            print(f"  Adjusted N: {r['auc_adjusted']['n']}")
            print(f"  Adjusted AUC: Our={r['auc_adjusted']['our']:.4f}, CADD={r['auc_adjusted']['cadd']:.4f}, REVEL={r['auc_adjusted']['revel']:.4f}, AM={r['auc_adjusted']['am']:.4f}")
    
    print("\nCreating figure...")
    create_figure(results)
    
    print("\nCreating summary table...")
    df = create_summary_table(results)
    print("\n", df.to_string(index=False))
    
    print("\nDone!")


if __name__ == "__main__":
    main()
