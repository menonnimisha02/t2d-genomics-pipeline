# Type 2 Diabetes Genomics & Polygenic Risk Score Pipeline

An end-to-end bioinformatics project integrating Type 2 Diabetes GWAS summary statistics with 1000 Genomes genotype data using Python, PLINK2, Linux and Bash.

## Project Summary

This project:

- processed 19.3 million T2D GWAS variant associations
- processed 12.1 million genome-wide variants from 2,504 individuals
- performed genome-wide PCA using 286,036 LD-pruned SNPs
- harmonised 9.25 million GWAS and genotype variants
- generated a clumping-and-thresholding PRS using 3,650 lead variants
- compared PRS distributions across AFR, AMR, EAS, EUR and SAS populations

## Main Results

- GWAS variants processed: **19,303,703**
- Genome-wide significant associations: **64,948**
- 1000 Genomes individuals: **2,504**
- Genome-wide SNPs: **12,057,349**
- SNPs used for PCA: **286,036**
- Harmonised GWAS–genotype variants: **9,251,622**
- PRS variants after LD clumping: **3,650**

## Figures

### Manhattan Plot
![Manhattan Plot](figures/t2d_manhattan.png)

### QQ Plot
![QQ Plot](figures/t2d_qq.png)

### Genome-wide PCA
![Genome-wide PCA](figures/genomewide_pca.png)

### PRS by Ancestry
![PRS by Ancestry](figures/t2d_prs_by_ancestry.png)

## Tools

- Python
- pandas
- NumPy
- Matplotlib
- PLINK2
- Linux / WSL
- Bash
- Git

## Important Limitation

1000 Genomes does not provide suitable Type 2 Diabetes phenotype data for prediction testing. Therefore, this project evaluates PRS distributions and cross-population portability rather than clinical prediction accuracy.

## Author

Nimisha Menon
