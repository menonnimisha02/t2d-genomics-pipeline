import pandas as pd
from pathlib import Path

INPUT = Path("data/gwas/raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip")
OUTPUT_DIR = Path("data/gwas/processed")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CHUNK_SIZE = 1_000_000

reader = pd.read_csv(
    INPUT,
    sep=r"\s+",
    compression="zip",
    chunksize=CHUNK_SIZE
)

total_rows = 0
significant_rows = 0
missing_rows = 0
invalid_p = 0
invalid_eaf = 0
invalid_alleles = 0
same_alleles = 0
invalid_se = 0

for i, chunk in enumerate(reader, start=1):
    total_rows += len(chunk)
    significant_rows += (chunk["Pval"] < 5e-8).sum()
    missing_rows += chunk.isna().any(axis=1).sum()
    invalid_p += ((chunk["Pval"] < 0) | (chunk["Pval"] > 1)).sum()
    invalid_eaf += ((chunk["EAF"] < 0) | (chunk["EAF"] > 1)).sum()
    valid_bases = ["A", "C", "G", "T"]
    invalid_alleles += (
    ~chunk["EffectAllele"].isin(valid_bases) |
    ~chunk["NonEffectAllele"].isin(valid_bases)).sum()
    same_alleles += (
    chunk["EffectAllele"] == chunk["NonEffectAllele"]).sum()
    invalid_se += (chunk["SE"] <= 0).sum()
    print(f"Processed chunk {i}: {len(chunk):,} rows | total: {total_rows:,}")

print(f"Total variant rows: {total_rows:,}")
print(f"Genome-wide significant variants: {significant_rows:,}")
print(f"Rows with missing values: {missing_rows:,}")
print(f"Invalid P-values: {invalid_p:,}")
print(f"Invalid EAF values: {invalid_eaf:,}")
print(f"Invalid allele rows: {invalid_alleles:,}")
print(f"Same effect/non-effect allele: {same_alleles:,}")
print(f"Invalid SE values: {invalid_se:,}")
