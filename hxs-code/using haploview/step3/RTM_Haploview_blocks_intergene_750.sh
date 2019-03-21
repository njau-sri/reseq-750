#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#$ -l mem_free=4G
#


# your code goes here
export OMP_NUM_THREADS=16

./rtm-gwas-snpldb --openmp --vcf /ds3512/home/haoxsh/haodata/SNPLDB/sri750_genotype_final_2745637snps.vcf --block intergeneGAM_all.block --maf 0.01 --out GAMETES_intergene_750
./rtm-gwas-snpldb --openmp --vcf /ds3512/home/haoxsh/haodata/SNPLDB/sri750_genotype_final_2745637snps.vcf --block intergeneGABR_all.block --maf 0.01 --out GABRIEL_intergene_750
./rtm-gwas-snpldb --openmp --vcf /ds3512/home/haoxsh/haodata/SNPLDB/sri750_genotype_final_2745637snps.vcf --block intergeneSPI_all.block --maf 0.01 --out SPINE_intergene_750

