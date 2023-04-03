#!/bin/bash

### Edit the path based on your own directories 

RAW=/Volumes/Path/8823-JE-RAW

TRIM=/Volumes/Path/8823-JE-TRIM

ADAPTER=/Volumes/Path/TruSeq_CD_adapter.txt

TRIMSOFT=/Users/Path/software/Trimmomatic-0.39/trimmomatic-0.39.jar

ALIGN=/Volumes/Path/8823-JE-bowtie

GENOME=/Volumes/Path/Genomes/hg19_ecoli/hg19_ecK12MG1655

for i in 0001 0002 0003 0004 0005 0006 0007 0008 0009 0010 0011 0012

do 

### Edit the Project number and initials based on your own files

bowtie2 -p 16 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-overlap --no-dovetail \
 -x ${GENOME} -1 ${TRIM}/8823-JE-${i}_R1.noadap.paired.txt -2 ${TRIM}/8823-JE-${i}_R2.noadap.paired.txt \
 -S ${ALIGN}/8823-JE-${i}.hg19ec.sam 

done