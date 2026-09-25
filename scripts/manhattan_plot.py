import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

INPUT = Path("data/gwas/raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip")
OUTPUT = Path("figures/t2d_manhattan.png")

CHUNK_SIZE = 1_000_000

reader = pd.read_csv(
    INPUT,
    sep=r"\s+",
    compression="zip",
    usecols=["Chromsome", "Position", "Pval"],
    chunksize=CHUNK_SIZE
)

chrom_max = {}

for chunk in reader:
    chunk = chunk[chunk["Chromsome"].between(1, 22)]

    for chrom, max_pos in chunk.groupby("Chromsome")["Position"].max().items():
        chrom = int(chrom)
        chrom_max[chrom] = max(chrom_max.get(chrom, 0), int(max_pos))

print(chrom_max)

offsets = {}
chrom_centers = {}
running_position = 0

for chrom in range(1, 23):
    offsets[chrom] = running_position
    chrom_centers[chrom] = running_position + chrom_max[chrom] / 2
    running_position += chrom_max[chrom]

print(offsets)

plot_parts = []

reader = pd.read_csv(
    INPUT,
    sep=r"\s+",
    compression="zip",
    usecols=["Chromsome", "Position", "Pval"],
    chunksize=CHUNK_SIZE
)

for chunk in reader:
    chunk = chunk[chunk["Chromsome"].between(1, 22)].copy()

    significant = chunk[chunk["Pval"] < 5e-8].copy()

    background = chunk[chunk["Pval"] >= 5e-8].sample(
        frac=0.01,
        random_state=42)

    plot_parts.append(
        pd.concat([background, significant]))

plot_df = pd.concat(plot_parts, ignore_index=True)

plot_df["PlotPos"] = (
    plot_df["Position"]
    + plot_df["Chromsome"].map(offsets))

plot_df["MinusLog10P"] = -np.log10(
    plot_df["Pval"].clip(lower=1e-300))

print(f"Points prepared for plotting: {len(plot_df):,}")

plot_df["Chromsome"] = plot_df["Chromsome"].astype(int)

plt.figure(figsize=(16, 6))

for chrom in range(1, 23):
    chrom_data = plot_df[plot_df["Chromsome"] == chrom]
    plt.scatter(
        chrom_data["PlotPos"],
        chrom_data["MinusLog10P"],
        s=4,
        alpha=0.6,
        label=str(chrom) if chrom <= 2 else None)

plt.axhline(-np.log10(5e-8), linestyle="--")

plt.xticks(
    [chrom_centers[c] for c in range(1, 23)],
    [str(c) for c in range(1, 23)],
    rotation=0)

plt.xlabel("Chromosome")
plt.ylabel("-log10(P-value)")
plt.title("T2D GWAS Manhattan Plot")
plt.tight_layout()
plt.savefig(OUTPUT, dpi=300)

print(f"Saved plot to {OUTPUT}")
