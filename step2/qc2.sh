#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#$ -l mem_free=4G
#

# your code goes here

VCFTOOLS=/share/apps/vcftools/vcftools-0.1.15

$VCFTOOLS --gzvcf sri750_2894720snps_genotype_imputed.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.01 \
    --recode --stdout | gzip -c > sri750_imputed_final.vcf.gz
