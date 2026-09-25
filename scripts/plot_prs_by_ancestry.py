import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

PRS_FILE = Path("results/t2d_prs.sscore")
META_FILE = Path("data/metadata/1000G_samples.panel")

OUT_TABLE = Path("results/t2d_prs_with_ancestry.tsv")
OUT_FIG = Path("figures/t2d_prs_by_ancestry.png")

# Read PRS results
prs = pd.read_csv(
    PRS_FILE,
    sep=r"\s+")

# Read 1000G ancestry metadata
meta = pd.read_csv(
    META_FILE,
    sep=r"\s+")

# Join ancestry labels
merged = prs.merge(
    meta,
    left_on="#IID",
    right_on="sample",
    how="left")

# Standardise PRS across all 1000G individuals
eur_scores = merged.loc[
    merged["super_pop"] == "EUR",
    "Beta_SUM"
]

eur_mean = eur_scores.mean()
eur_std = eur_scores.std()

merged["PRS_Z"] = (
    merged["Beta_SUM"] - eur_mean
) / eur_std

print(f"EUR reference mean raw PRS: {eur_mean:.4f}")
print(f"EUR reference SD raw PRS: {eur_std:.4f}")

merged.to_csv(
    OUT_TABLE,
    sep="\t",
    index=False)

print(f"Individuals scored: {len(merged):,}")
print(f"Missing ancestry labels: {merged['super_pop'].isna().sum()}")

print()
print("Mean standardised PRS by superpopulation:")
print(
    merged.groupby("super_pop")["PRS_Z"]
    .agg(["count", "mean", "std"]))

# Boxplot
order = ["AFR", "AMR", "EAS", "EUR", "SAS"]

data = [
    merged.loc[
        merged["super_pop"] == pop,
        "PRS_Z"]
    for pop in order
]

plt.figure(figsize=(9, 6))

plt.boxplot(
    data,
    tick_labels=order)

plt.axhline(0, linestyle="--")

plt.xlabel("1000 Genomes superpopulation")
plt.ylabel("T2D PRS (EUR-reference SD units)")
plt.title("European-derived T2D PRS across 1000 Genomes populations")

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=300)

print()
print(f"Saved plot to {OUT_FIG}")
