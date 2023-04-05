#!/bin/bash

## Edit path for each directory 

RAW=/Path/9631-RAW_CTCF

TRIM=/Path/9631-TRIM_CTCF

ALIGN=/Path/9631-ALIGN_CTCF

GENOME=/Path/Genomes/hg19_ecoli/hg19_ecK12MG1655	

BAM=/Path/9631-BAM_CTCF

## Edit project number and initials 

for i in 0001 0002 0003 0004 0005 0006 0007 0008 0009 0010 0011 0012

do 

echo "Creating bam file ${i}"

samtools view -S -b -@ 14  ${ALIGN}/9631-MB-${i}.hg19ec.sam -o ${BAM}/9631-MB-${i}.hg19ec.bam

echo "F4q10 bam file ${i}"	
samtools view -b -F 4 -q 10 -@ 14 ${BAM}/9631-MB-${i}.hg19ec.bam -o ${BAM}/9631-MB-${i}.hg19ec.F4q10.bam

echo "Sorting bam file ${i}"	

samtools sort -@ 12 ${BAM}/9631-MB-${i}.hg19ec.F4q10.bam -o ${BAM}/9631-MB-${i}.hg19ec.F4q10.sorted.bam

samtools index ${BAM}/9631-MB-${i}.hg19ec.F4q10.sorted.bam

echo "Creating hg19 bam file for ${i}"	

samtools view -@ 12 -bh ${BAM}/9631-MB-${i}.hg19ec.F4q10.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${BAM}/9631-MB-${i}.hg19.F4q10.sorted.bam

samtools index ${BAM}/9631-MB-${i}.hg19.F4q10.sorted.bam

done