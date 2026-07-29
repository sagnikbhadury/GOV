#!/usr/bin/env python3
"""
STANDALONE FDR ANALYSIS - No external dependencies
Uses only Python built-in modules
"""

import csv
import statistics
from pathlib import Path
from collections import defaultdict

print("\n" + "="*70)
print(" FDR ANALYSIS PIPELINE (Standalone - No External Dependencies)")
print("="*70 + "\n")

# Setup paths
base_path = Path(r"C:\Users\bhadury\University of Michigan Dropbox\Sagnik Bhadury\CODES\SBTurboCodes\H&E Lung TCGA Graph Dataset")
codes_dir = base_path / "codes"
results_dir = base_path / "CLAUDE_REVISIONS_JULY2026" / "analysis_results"

results_dir.mkdir(parents=True, exist_ok=True)
(results_dir / "supplementary_tables").mkdir(parents=True, exist_ok=True)

print(f"Base path: {base_path}")
print(f"Results directory: {results_dir}\n")

# Benjamini-Hochberg FDR correction (standalone implementation)
def benjamini_hochberg(p_values):
    """Apply Benjamini-Hochberg FDR correction"""
    n = len(p_values)
    if n == 0:
        return []

    # Create list of (p_value, original_index)
    indexed_pvals = [(p, i) for i, p in enumerate(p_values)]
    # Sort by p-value
    sorted_pvals = sorted(indexed_pvals, key=lambda x: x[0])

    # Calculate q-values
    q_values = [0] * n
    cummin = 1.0

    for rank, (p_val, orig_idx) in enumerate(reversed(sorted_pvals), 1):
        # BH formula: q = p * (n / rank)
        q_val = min(p_val * n / rank, 1.0)
        # Enforce monotonicity
        cummin = min(cummin, q_val)
        q_values[orig_idx] = cummin

    return q_values

# ========================================================================
# Load and process pathway PIPs
# ========================================================================

print("=== LOADING DATA ===\n")

def load_pips(filepath, feature_name):
    """Load PIP data from CSV"""
    data = []
    p_values = []

    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                pip = float(row.get('PIP', row.get(list(row.keys())[2])))
                data.append(row)
                # Convert PIP to p-value (conservative: 1 - PIP)
                p_values.append(1 - pip)
            except (ValueError, KeyError):
                continue

    return data, p_values

# Load all three layers
try:
    pathways_data, pathways_pvals = load_pips(
        codes_dir / "luad_Activity_lsGPRVS_pip_longformat.csv", "pathway"
    )
    kinases_data, kinases_pvals = load_pips(
        codes_dir / "luad_kinase_lsGPRVS_pip_longformat.csv", "kinase"
    )
    tfs_data, tfs_pvals = load_pips(
        codes_dir / "luad_tf_lsGPRVS_pip_longformat.csv", "tf"
    )

    print(f"✓ Pathways: {len(pathways_data)} rows loaded")
    print(f"✓ Kinases: {len(kinases_data)} rows loaded")
    print(f"✓ TFs: {len(tfs_data)} rows loaded\n")

except FileNotFoundError as e:
    print(f"ERROR: {e}")
    exit(1)

# ========================================================================
# Apply FDR correction
# ========================================================================

print("=== APPLYING FDR CORRECTION ===\n")

def apply_fdr(data, p_values, layer_name):
    """Apply FDR correction and return results"""
    q_values = benjamini_hochberg(p_values)

    # Combine data with q-values
    results = []
    pip_55 = 0
    fdr_05 = 0

    for i, row in enumerate(data):
        pip = float(row.get('PIP', row.get(list(row.keys())[2])))
        q_val = q_values[i]

        result = dict(row)
        result['p_value'] = p_values[i]
        result['q_value'] = q_val
        result['Significant_PIP05'] = 'True' if pip > 0.5 else 'False'
        result['Significant_FDR05'] = 'True' if q_val < 0.05 else 'False'

        results.append(result)

        if pip > 0.5:
            pip_55 += 1
        if q_val < 0.05:
            fdr_05 += 1

    # Sort by q-value
    results.sort(key=lambda x: x['q_value'])

    return results, pip_55, fdr_05

pathways_results, path_pip, path_fdr = apply_fdr(pathways_data, pathways_pvals, "Pathway")
kinases_results, kin_pip, kin_fdr = apply_fdr(kinases_data, kinases_pvals, "Kinase")
tfs_results, tf_pip, tf_fdr = apply_fdr(tfs_data, tfs_pvals, "TF")

print(f"PATHWAYS:")
print(f"  PIP > 0.5: {path_pip} associations")
print(f"  FDR < 0.05: {path_fdr} associations")
if path_pip > 0:
    print(f"  Reduction: {(path_pip - path_fdr) / path_pip * 100:.1f}%")

print(f"\nKINASES:")
print(f"  PIP > 0.5: {kin_pip} associations")
print(f"  FDR < 0.05: {kin_fdr} associations")
if kin_pip > 0:
    print(f"  Reduction: {(kin_pip - kin_fdr) / kin_pip * 100:.1f}%")

print(f"\nTRANSCRIPTION FACTORS:")
print(f"  PIP > 0.5: {tf_pip} associations")
print(f"  FDR < 0.05: {tf_fdr} associations")
if tf_pip > 0:
    print(f"  Reduction: {(tf_pip - tf_fdr) / tf_pip * 100:.1f}%")

# ========================================================================
# Save results
# ========================================================================

print("\n=== SAVING RESULTS ===\n")

def save_results(results, output_file):
    """Save results to CSV"""
    if not results:
        return

    fieldnames = list(results[0].keys())

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

save_results(pathways_results, results_dir / "supplementary_tables" / "Table_S2_Pathway_Associations_FDR.csv")
save_results(kinases_results, results_dir / "supplementary_tables" / "Table_S3_Kinase_Associations_FDR.csv")
save_results(tfs_results, results_dir / "supplementary_tables" / "Table_S4_TF_Associations_FDR.csv")

print("✓ Table_S2_Pathway_Associations_FDR.csv")
print("✓ Table_S3_Kinase_Associations_FDR.csv")
print("✓ Table_S4_TF_Associations_FDR.csv")

# ========================================================================
# Summary statistics
# ========================================================================

print("\n=== SUMMARY STATISTICS ===\n")

summary_data = [
    ['Metric', 'Count'],
    ['---', '---'],
    ['Pathway Associations (PIP > 0.5)', str(path_pip)],
    ['Pathway Associations (FDR < 0.05)', str(path_fdr)],
    ['---', '---'],
    ['Kinase Associations (PIP > 0.5)', str(kin_pip)],
    ['Kinase Associations (FDR < 0.05)', str(kin_fdr)],
    ['---', '---'],
    ['TF Associations (PIP > 0.5)', str(tf_pip)],
    ['TF Associations (FDR < 0.05)', str(tf_fdr)],
]

with open(results_dir / "Analysis_Summary_Statistics.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(summary_data)

for row in summary_data:
    print(f"{row[0]:<50} {row[1]}")

# ========================================================================
# Threshold comparison
# ========================================================================

print("\n=== THRESHOLD COMPARISON ===\n")

layers = ['Pathways', 'Kinases', 'Transcription Factors']
pip_counts = [path_pip, kin_pip, tf_pip]
fdr_counts = [path_fdr, kin_fdr, tf_fdr]

comparison_data = [['Layer', 'PIP > 0.5', 'FDR < 0.05', 'Reduction (%)']]

for layer, pip, fdr in zip(layers, pip_counts, fdr_counts):
    reduction = (pip - fdr) / pip * 100 if pip > 0 else 0
    comparison_data.append([layer, str(pip), str(fdr), f"{reduction:.1f}"])

with open(results_dir / "supplementary_tables" / "Threshold_Comparison.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(comparison_data)

for row in comparison_data:
    print(f"{row[0]:<25} {row[1]:<15} {row[2]:<15} {row[3]}")

# ========================================================================
# Top associations
# ========================================================================

print("\n=== TOP 10 ASSOCIATIONS BY FDR ===\n")

print("Top Pathway Associations (FDR < 0.05):")
for row in pathways_results[:10]:
    if float(row.get('q_value', 1)) < 0.05:
        print(f"  {row.get('HPC', row.get(list(row.keys())[0]))}: "
              f"PIP={float(row.get('PIP', 0)):.3f}, "
              f"q={float(row.get('q_value', 1)):.2e}")

print("\nTop Kinase Associations (FDR < 0.05):")
for row in kinases_results[:10]:
    if float(row.get('q_value', 1)) < 0.05:
        print(f"  {row.get('HPC', row.get(list(row.keys())[0]))}: "
              f"PIP={float(row.get('PIP', 0)):.3f}, "
              f"q={float(row.get('q_value', 1)):.2e}")

print("\nTop TF Associations (FDR < 0.05):")
for row in tfs_results[:10]:
    if float(row.get('q_value', 1)) < 0.05:
        print(f"  {row.get('HPC', row.get(list(row.keys())[0]))}: "
              f"PIP={float(row.get('PIP', 0)):.3f}, "
              f"q={float(row.get('q_value', 1)):.2e}")

# ========================================================================
# Final summary
# ========================================================================

print("\n" + "="*70)
print(" ANALYSIS COMPLETE")
print("="*70)

print(f"\n✓ All results saved to: {results_dir}")

print("\nGENERATED FILES:")
print(f"  - Table_S2_Pathway_Associations_FDR.csv")
print(f"  - Table_S3_Kinase_Associations_FDR.csv")
print(f"  - Table_S4_TF_Associations_FDR.csv")
print(f"  - Threshold_Comparison.csv")
print(f"  - Analysis_Summary_Statistics.csv")

print("\nKEY FINDINGS:")
print(f"  Pathways: {(path_pip - path_fdr) / path_pip * 100:.1f}% reduction with FDR correction")
print(f"  Kinases: {(kin_pip - kin_fdr) / kin_pip * 100:.1f}% reduction with FDR correction")
print(f"  TFs: {(tf_pip - tf_fdr) / tf_pip * 100:.1f}% reduction with FDR correction")

print("\nNEXT STEPS:")
print("  1. ✓ FDR correction complete")
print("  2. Review generated tables")
print("  3. Update manuscript with real numbers")
print("  4. Set up GitHub repository")
print("  5. Run permutation testing (03_PERMUTATION_TESTING_PARALLEL.R)")

print("\n" + "="*70 + "\n")
