#!/bin/bash

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

for chr in {2..22}
do
    echo "======================================"
    echo "Starting chromosome $chr"
    echo "======================================"

    VCF="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

    # Download chromosome
    wget -c \
      "${BASE_URL}/${VCF}" \
      -P data/raw/

    # Convert and filter with PLINK2
    plink2 \
      --vcf "data/raw/${VCF}" \
      --snps-only just-acgt \
      --max-alleles 2 \
      --maf 0.01 \
      --set-all-var-ids '@:#:$r:$a' \
      --make-pgen \
      --out "data/processed/chr${chr}_prs"

    echo "Finished chromosome $chr"
done

echo "======================================"
echo "All chromosomes 2-22 finished"
echo "======================================"
