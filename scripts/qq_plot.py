import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

INPUT = Path("data/gwas/raw/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip")
OUTPUT = Path("figures/t2d_qq.png")

CHUNK_SIZE = 1_000_000

pvalue_parts = []
zero_pvalues = 0

reader = pd.read_csv(
    INPUT,
    sep=r"\s+",
    compression="zip",
    usecols=["Pval"],
    chunksize=CHUNK_SIZE)

for i, chunk in enumerate(reader, start=1):
    p = chunk["Pval"].to_numpy(dtype=float)
    zero_pvalues += (p == 0).sum()
    # QQ plot cannot take log10(0), so keep positive P-values
    p = p[(p > 0) & (p <= 1)]
    pvalue_parts.append(p)
    print(f"Processed chunk {i}")

pvalues = np.concatenate(pvalue_parts)

print(f"Positive P-values used: {len(pvalues):,}")
print(f"P-values equal to zero: {zero_pvalues:,}")

# Sort from smallest P-value to largest
pvalues.sort()

n = len(pvalues)

observed = -np.log10(pvalues)

expected_p = (np.arange(1, n + 1) - 0.5) / n
expected = -np.log10(expected_p)

# Plot 100,000 evenly spaced points instead of all ~19 million
indices = np.linspace(
    0,
    n - 1,
    min(100000, n),
    dtype=int)

x = expected[indices]
y = observed[indices]

plt.figure(figsize=(7, 7))

plt.scatter(
    x,
    y,
    s=5,
    alpha=0.5)

max_expected = x.max()

plt.plot(
    [0, max_expected],
    [0, max_expected],
    linestyle="--")

plt.xlim(0, max_expected * 1.05)
plt.ylim(bottom=0)

plt.xlabel("Expected -log10(P-value)")
plt.ylabel("Observed -log10(P-value)")
plt.title("T2D GWAS QQ Plot")

plt.tight_layout()
plt.savefig(OUTPUT, dpi=300)

print(f"Saved QQ plot to {OUTPUT}")
