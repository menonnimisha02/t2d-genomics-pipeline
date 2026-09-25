import pandas as pd
from pathlib import Path

CLUMPS = Path("results/t2d_clumped.clumps")
HARMONISED = Path("results/t2d_1000g_harmonised.tsv")
OUTPUT = Path("results/t2d_prs_weights.tsv")

# Get the 3,650 lead variants from the clumping result
clumps = pd.read_csv(
    CLUMPS,
    sep=r"\s+",
    usecols=["ID"])

lead_ids = set(clumps["ID"])

print(f"Lead clump variants: {len(lead_ids):,}")

matched_parts = []

# Search the large harmonised GWAS file in chunks
for chunk in pd.read_csv(
    HARMONISED,
    sep="\t",
    chunksize=500_000):
    found = chunk[chunk["ID"].isin(lead_ids)]

    if not found.empty:
        matched_parts.append(found)

weights = pd.concat(matched_parts, ignore_index=True)

# If an ID somehow occurs more than once, keep the strongest association
weights = (
    weights
    .sort_values("Pval")
    .drop_duplicates("ID", keep="first"))

score_file = weights[
    ["ID", "EffectAllele", "Beta"]
].copy()

score_file.to_csv(
    OUTPUT,
    sep="\t",
    index=False)

print(f"Variants found for scoring: {len(score_file):,}")
print(f"Lead variants missing: {len(lead_ids) - len(score_file):,}")
print(f"Saved PRS weights to: {OUTPUT}")
