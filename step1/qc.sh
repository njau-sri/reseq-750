#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#$ -l mem_free=4G
#

# your code goes here

# sri reseq 5x 750

PYTHON2=/opt/python/bin/python2
VCFTOOLS=/share/apps/vcftools/vcftools-0.1.15

# MAF >=0.01 and MISSING <=0.1
$VCFTOOLS --vcf GATK.final.snp_change.vcf --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --max-missing 0.9 --recode --out filtered.GATK  || exit 1
$VCFTOOLS --vcf mpileup.final.snp_change.vcf --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --max-missing 0.9 --recode --out filtered.mpileup  || exit 1

# calculates discordance on a site by site basis
$VCFTOOLS --vcf filtered.GATK.recode.vcf --diff filtered.mpileup.recode.vcf --diff-site-discordance  || exit 1

# common SNP and alleles match
awk '{if(NR==1||($3=="B"&&$4=="1"))print}' out.diff.sites > out.diff.sites.common  || exit 1

# discordance histogram
# cut -f 7 out.diff.sites.common | tail -n +2 > z.txt
# x <- scan("z.txt")
# h <- hist(x, breaks=seq(0,1,0.01))
# write.table(cbind(h$breaks,c(h$counts,NA)),"freq.txt",quote=F,sep="\t")

# mumber of genotypes called in both files >=675 =750*0.9
# discordance rate <= 0.1
awk '{if(NR==1||($5>=675&&$7<=0.1))print}' out.diff.sites.common > out.diff.sites.common.good  || exit 1

# extract CHROM, POS
cut -f 1,2 out.diff.sites.common.good | tail -n +2 > pos.txt  || exit 1

# filter common SNP
$VCFTOOLS --vcf filtered.GATK.recode.vcf --positions pos.txt --recode --out gatk.common.good  || exit 1
$VCFTOOLS --vcf filtered.mpileup.recode.vcf --positions pos.txt --recode --out mpileup.common.good  || exit 1

# treat discordant gentoype as missing
$PYTHON2 drop_discord.py gatk.common.good.recode.vcf mpileup.common.good.recode.vcf common.good.vcf  || exit 1

# again MAF >=0.01 and MISSING <=0.1
$VCFTOOLS --vcf common.good.vcf --maf 0.01 --max-missing 0.9 --recode --out common.good.filtered  || exit 1

# filter heterozygous genotype: HET <=0.1
$PYTHON2 filter_het.py common.good.filtered.recode.vcf 0.1 > final.vcf
