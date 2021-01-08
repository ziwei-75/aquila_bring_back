#!/bin/sh

mkdir bam_index
for i in {1..22}
do 
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr$i.fa.gz -P ./bam_index
	gunzip bam_index/chr$i.fa.gz
	minimap2 -d bam_index/chr$i.mmi bam_index/chr$i.fa
done

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz -P ./bam_index
gunzip bam_index/chrX.fa.gz
minimap2 -d bam_index/chrX.mmi bam_index/chrX.fa
rm ./bam_index/*.fa

