#!/usr/bin/env bash
set -euo pipefail

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

mkdir -p data/raw data/processed

for chr in {1..22}; do
    VCF="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

    echo "Processing chromosome ${chr}"

    wget -c "${BASE_URL}/${VCF}" -P data/raw/

    gzip -t "data/raw/${VCF}"

    plink2 \
        --vcf "data/raw/${VCF}" \
        --snps-only just-acgt \
        --max-alleles 2 \
        --maf 0.01 \
        --set-all-var-ids '@:#:$r:$a' \
        --make-pgen \
        --out "data/processed/chr${chr}_prs"
done

echo "Finished processing chromosomes 1-22"
