#!/usr/bin/env bash
set -euo pipefail

mkdir -p results figures

# Merge chromosomes 1-22 into one genome-wide PLINK dataset
plink2 \
    --pfile data/processed/chr1_prs \
    --pmerge-list config/autosomes_merge_list.txt \
    --make-pgen \
    --out data/processed/1000g_prs_genomewide

# Select common SNPs and perform LD pruning for PCA
plink2 \
    --pfile data/processed/1000g_prs_genomewide \
    --maf 0.05 \
    --indep-pairwise 500kb 0.2 \
    --out results/genomewide_ld

# Genome-wide PCA
plink2 \
    --pfile data/processed/1000g_prs_genomewide \
    --extract results/genomewide_ld.prune.in \
    --pca 10 \
    --out results/genomewide_pca

python3 scripts/plot_genomewide_pca.py

# Harmonise GWAS variants with the 1000 Genomes target dataset
python3 scripts/harmonise_gwas_1000g.py

# Prepare EUR LD reference samples and GWAS clumping input
python3 scripts/prepare_prs_inputs.py

# LD clumping using European 1000 Genomes samples
plink2 \
    --pfile data/processed/1000g_prs_genomewide \
    --keep results/eur_samples.txt \
    --clump results/t2d_clump_input.tsv \
    --clump-p1 1e-5 \
    --clump-p2 1e-5 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --out results/t2d_clumped

# Extract GWAS effect sizes for the clumped variants
python3 scripts/make_prs_weights.py

# Calculate PRS for all 1000 Genomes individuals
plink2 \
    --pfile data/processed/1000g_prs_genomewide \
    --score results/t2d_prs_weights.tsv 1 2 3 header-read cols=+scoresums \
    --out results/t2d_prs

python3 scripts/plot_prs_by_ancestry.py

echo "Genome-wide PCA and PRS analysis complete"
