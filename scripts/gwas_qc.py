import pandas as pd
from pathlib import Path

INPUT = Path("data/gwas/raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip")
OUTPUT = Path("results/gwas_summary.txt")

CHUNK_SIZE = 1_000_000

total_rows = 0
significant_rows = 0
missing_rows = 0
invalid_p = 0
zero_p = 0
invalid_eaf = 0
invalid_alleles = 0
same_alleles = 0
invalid_se = 0

valid_bases = {"A", "C", "G", "T"}

reader = pd.read_csv(
    INPUT,
    sep=r"\s+",
    compression="zip",
    chunksize=CHUNK_SIZE)

for i, chunk in enumerate(reader, start=1):
    total_rows += len(chunk)

    significant_rows += (chunk["Pval"] < 5e-8).sum()
    missing_rows += chunk.isna().any(axis=1).sum()

    invalid_p += (
        (chunk["Pval"] < 0) |
        (chunk["Pval"] > 1)).sum()

    zero_p += (chunk["Pval"] == 0).sum()

    invalid_eaf += (
        (chunk["EAF"] < 0) |
        (chunk["EAF"] > 1)).sum()

    invalid_alleles += (
        ~chunk["EffectAllele"].isin(valid_bases) |
        ~chunk["NonEffectAllele"].isin(valid_bases)).sum()

    same_alleles += (
        chunk["EffectAllele"] ==
        chunk["NonEffectAllele"]).sum()

    invalid_se += (chunk["SE"] <= 0).sum()

    print(
        f"Processed chunk {i}: "
        f"{len(chunk):,} rows | "
        f"total {total_rows:,}")


summary = f"""T2D GWAS Summary Statistics QC

Total variants processed: {total_rows:,}
Genome-wide significant variant associations (P < 5e-8): {significant_rows:,}

QC results:
- Missing-value rows: {missing_rows:,}
- Invalid P-values: {invalid_p:,}
- P-values recorded as zero: {zero_p:,}
- Invalid EAF values: {invalid_eaf:,}
- Invalid allele rows: {invalid_alleles:,}
- Same effect/non-effect allele rows: {same_alleles:,}
- Invalid SE values: {invalid_se:,}
"""

print()
print(summary)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(summary)

print(f"Saved QC summary to: {OUTPUT}")
