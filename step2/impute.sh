#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#$ -l mem_free=4G
#

PYTHON=/opt/python/bin/python2

# your code goes here

$PYTHON impute_soybean_v3.py ../../sri750_2894720snps_genotype.vcf.gz
