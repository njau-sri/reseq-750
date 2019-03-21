#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#$ -l mem_free=4G
#

PYTHON=/opt/python/bin/python2
VCFTOOLS=/share/apps/vcftools/vcftools-0.1.15
PLINK=/share/apps/plink/plink-1.90b4.3
HAPLOVIEW="java -jar /share/apps/haploview/4.2/haploview.jar -n"

# your code goes here

#按照soybase.org上提供的大豆基因组基因信息，得到基因组上共计基因间区53645个(删除重复后)
for i in $(seq 0 6)
do
	#block.txt文件为基因间区在基因组上物理位置文件，格式为 “染色体 起始位置 结束为止 区间长度”
	#snp.txt文件为全基因组上SNP的信息文件，格式示例“Chr01_4485”,文件中每行一个SNP的位点信息
	#生成文件$i.txt是每个基因间区中包含的SNP信息，编号代表基因组上第几个基因间区
	$PYTHON match_intergene_and_snp.py example_block.txt $i example_snp.txt $i.txt
	#sri750_genotype_final_2745637snps.vcf文件是测序完成之后，通过MAF过滤，去杂合，填补之后的SNP基因型VCF文件
	#$i.txt文件是上一步脚本生成的每个基因间区中SNP信息，通过VCFTOOLS软件的--snps命令将每个基因间区单独call出一个SNP文件，
	#文件中包含每个基因间区的全部SNP信息
	$VCFTOOLS --vcf /ds3512/home/haoxsh/haodata/SNPLDB/sri750_genotype_final_2745637snps.vcf --snps $i.txt --recode --out intergene_snp_$i
	#利用PLINK软件将上一步生成的基因间区vcf文件转化成Haploview需要的.PED的格式，最终生成
	$PLINK --vcf  intergene_snp_$i.recode.vcf --recode --out intergene_$i
	#在生成的.ped文件中的第六列为phnotype信息，plink以-9表示，但是haploview只识别0，所以要把.ped文件第六列改为0，
	#并重新写入到intergene_$i.x.ped中。
    $PYTHON change_nine_to_zero.py intergene_$i.ped intergene_$i.x.ped
	#plink生成的map文件也不可直接使用，haploview要求格式为“SNPID position”,自写Python脚本将haploview所需map文件内容写入intergene_$i.x.map
	$PYTHON make_infofile.py $i.txt intergene_$i.x.map
	$HAPLOVIEW -memory 20480 -pedfile intergene_$i.x.ped -info intergene_$i.x.map -log intergene_$i.log\
	-hwcutoff 0 -minMAF 0.01 -maxDistance 100 -blockMAFThresh 0.01 -blockoutput ALL -dprime -check
done
