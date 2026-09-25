import pandas as pd
from pathlib import Path

GWAS = Path("data/gwas/raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip")
PVAR = Path("data/processed/1000g_prs_genomewide.pvar")

OUTPUT = Path("results/t2d_1000g_harmonised.tsv")

CHUNK_SIZE = 500_000

print("Loading 1000G variant IDs...")

target_ids = set()

with open(PVAR) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.split("\t")

        variant_id = fields[2]
        target_ids.add(variant_id)

print(f"1000G variants loaded: {len(target_ids):,}")

matched_parts = []
total_gwas = 0
matched = 0
ambiguous_removed = 0

reader = pd.read_csv(
    GWAS,
    sep=r"\s+",
    compression="zip",
    chunksize=CHUNK_SIZE)

for i, chunk in enumerate(reader, start=1):

    total_gwas += len(chunk)

    # Remove strand-ambiguous SNPs
    allele_pair = (
        chunk["EffectAllele"] +
        chunk["NonEffectAllele"])

    ambiguous = allele_pair.isin(
        ["AT", "TA", "CG", "GC"])

    ambiguous_removed += ambiguous.sum()

    chunk = chunk[~ambiguous].copy()

    # Two possible IDs depending on REF/ALT orientation
    id1 = (
        chunk["Chromsome"].astype(str) + ":" +
        chunk["Position"].astype(str) + ":" +
        chunk["EffectAllele"] + ":" +
        chunk["NonEffectAllele"])

    id2 = (
        chunk["Chromsome"].astype(str) + ":" +
        chunk["Position"].astype(str) + ":" +
        chunk["NonEffectAllele"] + ":" +
        chunk["EffectAllele"])

    match1 = id1.isin(target_ids)
    match2 = id2.isin(target_ids)

    keep = match1 | match2

    matched_chunk = chunk.loc[keep].copy()

    matched_chunk["ID"] = id1.loc[keep]

    # If the reversed orientation matched, use that target ID instead
    reverse_indices = match2[keep]

    matched_chunk.loc[reverse_indices, "ID"] = id2.loc[keep][reverse_indices]

    matched_parts.append(matched_chunk)

    matched += len(matched_chunk)

    print(
        f"Chunk {i}: matched {len(matched_chunk):,} | "
        f"total matched {matched:,}")

harmonised = pd.concat(
    matched_parts,
    ignore_index=True)

harmonised = harmonised[
    [
        "ID",
        "Chromsome",
        "Position",
        "EffectAllele",
        "NonEffectAllele",
        "Beta",
        "SE",
        "EAF",
        "Pval"]
]

harmonised.to_csv(
    OUTPUT,
    sep="\t",
    index=False)

print()
print(f"Total GWAS variants: {total_gwas:,}")
print(f"Ambiguous SNPs removed: {ambiguous_removed:,}")
print(f"Matched variants: {len(harmonised):,}")
print(f"Saved to: {OUTPUT}")
