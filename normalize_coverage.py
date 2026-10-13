#!/usr/bin/env python3
"""
normalize_coverage.py
Normalisation de la matrice de couverture mosdepth par la somme des coverages
sur les régions BED uniquement (per-million scaling).
coverage_norm = coverage_brute * 1_000_000 / somme_coverage_sample
Usage: python3 normalize_coverage.py <coverage_matrix.tsv>
"""

import os, sys

if len(sys.argv) == 2:
    FINAL = sys.argv[1]
else:
    print("Usage: python3 normalize_coverage.py <coverage_matrix.tsv>")
    sys.exit(1)

FINAL_RAW  = FINAL.replace(".tsv", "_raw.tsv")
FINAL_NORM = FINAL.replace(".tsv", "_normalized.tsv")

# ── Lecture ────────────────────────────────────────────────────────────────────
print(f"Lecture : {FINAL}")
with open(FINAL) as f:
    lines = f.readlines()

header = lines[0].strip().split('\t')
# Les 4 premières colonnes sont chr, start, end, gene
meta_cols = header[:4]
samples   = header[4:]
print(f"  {len(lines)-1} régions x {len(samples)} samples")

regions = []
matrix  = []
for line in lines[1:]:
    parts = line.strip().split('\t')
    regions.append(parts[:4])  # chr, start, end, gene
    values = []
    for v in parts[4:]:
        try:
            values.append(float(v))
        except ValueError:
            values.append(None)
    matrix.append(values)

# ── Sommes par sample ──────────────────────────────────────────────────────────
sums = []
for j in range(len(samples)):
    total = sum(matrix[i][j] for i in range(len(regions)) if matrix[i][j] is not None)
    sums.append(total)

print("\nFacteurs de normalisation (somme / 1M) :")
for j, s in enumerate(samples):
    print(f"  {s} : {sums[j]:.2f}")

# ── Normalisation ──────────────────────────────────────────────────────────────
matrix_norm = []
for i, row in enumerate(matrix):
    row_norm = []
    for j, val in enumerate(row):
        if val is None or sums[j] == 0:
            row_norm.append("NA")
        else:
            row_norm.append(f"{val * 1_000_000 / sums[j]:.4f}")
    matrix_norm.append(row_norm)

# ── Sauvegarde brute ───────────────────────────────────────────────────────────
os.rename(FINAL, FINAL_RAW)
print(f"\nValeurs brutes : {FINAL_RAW}")

# ── Écriture normalisée ────────────────────────────────────────────────────────
with open(FINAL_NORM, 'w') as out:
    out.write("\t".join(header) + "\n")
    for i, meta in enumerate(regions):
        out.write("\t".join(meta + matrix_norm[i]) + "\n")

print(f"Matrice normalisée : {FINAL_NORM}")
