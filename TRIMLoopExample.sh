#!/bin/bash

### Edit the path based on your own directories 

RAW=/Volumes/Path/8823-JE-RAW

TRIM=/Volumes/Path/8823-JE-TRIM

ADAPTER=/Volumes/Path/TruSeq_CD_adapter.txt

TRIMSOFT=/Users/Path/software/Trimmomatic-0.39/trimmomatic-0.39.jar

for i in 0001 0002 0003 0004 0005 0006 0007 0008 0009

do 

echo "Trimming ${i}"

### Edit the Project number and initials based on your own files
java -classpath $TRIMSOFT org.usadellab.trimmomatic.TrimmomaticPE -phred33 -thread 16 \
	${RAW}/8823-JE-${i}_S1_L005_R1_001.fastq.gz ${RAW}/8823-JE-${i}_S1_L005_R2_001.fastq.gz \
		${TRIM}/8823-JE-${i}_R1.noadap.paired.txt ${TRIM}/8823-JE-${i}_R1.noadap.unpaired.txt \
			${TRIM}/8823-JE-${i}_R2.noadap.paired.txt ${TRIM}/8823-JE-${i}_R2.noadap.unpaired.txt ILLUMINACLIP:$ADAPTER:2:30:7 LEADING:15 TRAILING:15 MINLEN:15

done 
