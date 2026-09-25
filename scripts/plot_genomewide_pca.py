import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

PCA_FILE = Path("results/genomewide_pca.eigenvec")
META_FILE = Path("data/metadata/1000G_samples.panel")

OUT_TABLE = Path("results/genomewide_pca_with_population.tsv")
OUT_FIG = Path("figures/genomewide_pca.png")

# Read PCA results
pca = pd.read_csv(
    PCA_FILE,
    sep=r"\s+")

# Read 1000 Genomes population metadata
meta = pd.read_csv(
    META_FILE,
    sep=r"\s+")

# PLINK may call the sample column #IID
iid_col = "#IID" if "#IID" in pca.columns else "IID"

# Join PCA coordinates to population labels
merged = pca.merge(
    meta,
    left_on=iid_col,
    right_on="sample",
    how="left")

merged.to_csv(
    OUT_TABLE,
    sep="\t",
    index=False)

print(f"Individuals in PCA: {len(merged):,}")
print(f"Missing population labels: {merged['super_pop'].isna().sum()}")

# Plot PC1 vs PC2
plt.figure(figsize=(9, 7))

for population, group in merged.groupby("super_pop"):
    plt.scatter(
        group["PC1"],
        group["PC2"],
        s=18,
        alpha=0.7,
        label=population)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("1000 Genomes Genome-wide PCA")
plt.legend(title="Superpopulation")

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=300)

print(f"Saved PCA plot to {OUT_FIG}")
