#!/usr/bin/env python3
"""
merge_coverage.py
Merge des fichiers .regions.bed.gz de mosdepth en une matrice de couverture TSV.
- La couverture est toujours la dernière colonne (index -1)
- Le nom de gène est extrait de la colonne 4 du BED si présente
- Le tissu est extrait du nom du sample (tag avant Aligned, sauf Run\d+)
- Le nom du sample est coupé avant Aligned
Usage: python3 merge_coverage.py <mosdepth_results_dir> <bed_file> <output_matrix.tsv>
"""

import gzip, os, glob, sys, re

if len(sys.argv) == 4:
    OUTDIR = sys.argv[1]
    BED    = sys.argv[2]
    FINAL  = sys.argv[3]
else:
    print("Usage: python3 merge_coverage.py <mosdepth_results_dir> <bed_file> <output_matrix.tsv>")
    sys.exit(1)

# ── Nom court du sample ────────────────────────────────────────────────────────
def get_short_name(sample_name):
    # Cherche le pattern: base-TAGAligned
    m = re.search(r'^(.+)-([^-]+)Aligned', sample_name)
    if m:
        base = m.group(1)   # tout avant -TAGAligned
        tag  = m.group(2)
        if re.match(r'^Run\d+$', tag, re.IGNORECASE):
            return base             # ex: 6622CY000272  (Run = technique, ignoré)
        return f"{base}_{tag}"      # ex: 6624CY000712_PAX
    # Pas de pattern Aligned reconnu -> on garde le nom brut
    return sample_name

# ── Lecture du BED pour les noms de gènes ─────────────────────────────────────
gene_names = {}
print(f"Lecture BED : {BED}")
with open(BED) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        key = (parts[0], parts[1], parts[2])
        gene_names[key] = parts[3] if len(parts) >= 4 else "NA"
print(f"  {len(gene_names)} régions avec noms de gènes")

# ── Lister les fichiers mosdepth ───────────────────────────────────────────────
files = sorted(glob.glob(os.path.join(OUTDIR, "*.regions.bed.gz")))
print(f"Fichiers trouvés : {len(files)}")
if not files:
    print("ERREUR : aucun fichier .regions.bed.gz trouvé.")
    sys.exit(1)

raw_samples = [os.path.basename(f).replace(".regions.bed.gz", "") for f in files]
samples     = [get_short_name(s) for s in raw_samples]

print("Samples :")
for r, s in zip(raw_samples, samples):
    print(f"  {r} -> {s}")

# ── Lecture des fichiers mosdepth ──────────────────────────────────────────────
data = {}
all_regions = None

for f, sample in zip(files, samples):
    sample_data = {}
    region_order = []
    with gzip.open(f, 'rt') as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom, start, end = parts[0], parts[1], parts[2]
            coverage = parts[-1]
            key = (chrom, start, end)
            sample_data[key] = coverage
            region_order.append(key)
    data[sample] = sample_data
    if all_regions is None:
        all_regions = region_order

print(f"Régions : {len(all_regions)}")

# ── Écriture ───────────────────────────────────────────────────────────────────
print(f"Écriture : {FINAL}")
with open(FINAL, 'w') as out:
    header = ["chr", "start", "end", "gene"] + samples
    out.write("\t".join(header) + "\n")
    for key in all_regions:
        chrom, start, end = key
        gene = gene_names.get(key, "NA")
        coverages = [data[s].get(key, "NA") for s in samples]
        out.write("\t".join([chrom, start, end, gene] + coverages) + "\n")

print(f"Terminé : {len(all_regions)} régions x {len(samples)} samples")
print(f"Matrice : {FINAL}")
