import pandas as pd
from pathlib import Path

INPUT = Path("data/gwas/raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip")

SIG_OUTPUT = Path("results/t2d_significant_variants.tsv.gz")
TOP_OUTPUT = Path("results/t2d_top20_variants.tsv")
ZERO_OUTPUT = Path("results/t2d_zero_pvalue_variants.tsv")

CHUNK_SIZE = 1_000_000

significant_parts = []
zero_parts = []

reader = pd.read_csv(
    INPUT,
    sep=r"\s+",
    compression="zip",
    chunksize=CHUNK_SIZE)

for i, chunk in enumerate(reader, start=1):

    # Keep genome-wide significant variants
    significant = chunk[chunk["Pval"] < 5e-8].copy()

    if not significant.empty:
        significant_parts.append(significant)

    # Keep variants whose P-value is recorded as exactly zero
    zero_rows = chunk[chunk["Pval"] == 0].copy()

    if not zero_rows.empty:
        zero_parts.append(zero_rows)

    print(f"Processed chunk {i}")

# Combine all significant variants
sig_df = pd.concat(significant_parts, ignore_index=True)

# Rename the typo from the original dataset for cleaner output
sig_df = sig_df.rename(columns={"Chromsome": "Chromosome"})

# Make an easy-to-read variant identifier
sig_df["Variant"] = (
    sig_df["Chromosome"].astype(str)
    + ":"
    + sig_df["Position"].astype(str)
    + ":"
    + sig_df["EffectAllele"]
    + ":"
    + sig_df["NonEffectAllele"])

# Sort strongest associations first
sig_df = sig_df.sort_values("Pval")

# Save all significant variants
sig_df.to_csv(
    SIG_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip")

# Top 20 variants with positive P-values
top20 = sig_df[sig_df["Pval"] > 0].head(20)

top20.to_csv(
    TOP_OUTPUT,
    sep="\t",
    index=False)

# Save P=0 variants separately
if zero_parts:
    zero_df = pd.concat(zero_parts, ignore_index=True)
    zero_df = zero_df.rename(columns={"Chromsome": "Chromosome"})

    zero_df["Variant"] = (
        zero_df["Chromosome"].astype(str)
        + ":"
        + zero_df["Position"].astype(str)
        + ":"
        + zero_df["EffectAllele"]
        + ":"
        + zero_df["NonEffectAllele"])

    zero_df.to_csv(
        ZERO_OUTPUT,
        sep="\t",
        index=False)

print()
print(f"Genome-wide significant variants: {len(sig_df):,}")
print(f"Saved significant variants to: {SIG_OUTPUT}")
print(f"Saved top 20 positive-P variants to: {TOP_OUTPUT}")
print(f"Saved P=0 variants to: {ZERO_OUTPUT}")
