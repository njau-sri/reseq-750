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

$VCFTOOLS --gzvcf ../sri750_2894720snps_genotype.vcf.gz --gzdiff sri750_genotype_qc_2894720snps_imputed.vcf.gz --diff-site-discordance
