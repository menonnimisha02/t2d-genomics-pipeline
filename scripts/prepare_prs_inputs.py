import pandas as pd
from pathlib import Path

META = Path("data/metadata/1000G_samples.panel")
HARMONISED = Path("results/t2d_1000g_harmonised.tsv")

EUR_OUTPUT = Path("results/eur_samples.txt")
CLUMP_OUTPUT = Path("results/t2d_clump_input.tsv")

CHUNK_SIZE = 500_000


# Create European reference sample list
meta = pd.read_csv(META, sep=r"\s+")

eur_samples = meta.loc[
    meta["super_pop"] == "EUR",
    ["sample"]].copy()

eur_samples.columns = ["#IID"]

eur_samples.to_csv(
    EUR_OUTPUT,
    sep="\t",
    index=False)

print(f"EUR samples: {len(eur_samples):,}")


# Create clumping input from harmonised GWAS variants
clump_parts = []

for chunk in pd.read_csv(
    HARMONISED,
    sep="\t",
    usecols=["ID", "Pval"],
    chunksize=CHUNK_SIZE):
    selected = chunk.loc[
        chunk["Pval"] < 1e-5,
        ["ID", "Pval"]].copy()

    if not selected.empty:
        clump_parts.append(selected)

clump_input = pd.concat(
    clump_parts,
    ignore_index=True)

clump_input = (
    clump_input
    .sort_values("Pval")
    .drop_duplicates("ID", keep="first"))

clump_input = clump_input.rename(
    columns={"Pval": "P"})

clump_input.to_csv(
    CLUMP_OUTPUT,
    sep="\t",
    index=False)

print(f"Clumping candidates: {len(clump_input):,}")
print(f"Saved EUR samples to: {EUR_OUTPUT}")
print(f"Saved clumping input to: {CLUMP_OUTPUT}")
