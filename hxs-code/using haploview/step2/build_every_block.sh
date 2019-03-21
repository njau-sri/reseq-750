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

for i in $(seq 0 6)
do
	$PYTHON combine_haploview_result.py $i.txt  intergene_$i.GAMsnp.txt intergene_$i.x.ped.4GAMblocks  intergene_$i.GAMblock.txt
	$PYTHON combine_haploview_result.py $i.txt  intergene_$i.GABRsnp.txt intergene_$i.x.ped.GABRIELblocks  intergene_$i.GABRblock.txt
	$PYTHON combine_haploview_result.py $i.txt  intergene_$i.SPIsnp.txt intergene_$i.x.ped.SPINEblocks  intergene_$i.SPIblock.txt
	for j in $(seq 0 5)#实际计算时一个intergene区域可以最多得到30个block
	do
	$PYTHON test.py intergene_$i.GAMblock.txt $j intergene_$i.GAMsnp.txt intergeneGAM_$i.$j.txt intergeneGAM_$i.$j.block
	$PYTHON test.py intergene_$i.GABRblock.txt $j intergene_$i.GABRsnp.txt intergeneGABR_$i.$j.txt intergeneGABR_$i.$j.block
	$PYTHON test.py intergene_$i.SPIblock.txt $j intergene_$i.SPIsnp.txt intergeneSPI_$i.$j.txt intergeneSPI_$i.$j.block
	done
done

for i in $(seq 0 6)
do
	for j in $(seq 0 5)#实际计算时一个intergene区域可以最多得到30个block
	do
	$PYTHON combine_all_block_together.py intergeneGAM_$i.$j.block intergeneGAM_all.block
	$PYTHON combine_all_block_together.py intergeneGABR_$i.$j.block intergeneGABR_all.block
	$PYTHON combine_all_block_together.py intergeneSPI_$i.$j.block intergeneSPI_all.block
	done
done
