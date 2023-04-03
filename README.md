# CUT_RUN: Hiebert Lab How to
 A vignette of how to analyze CUT&RUN data for the Hiebert Lab 
## Trimmomatic 
### manual: http://www.usadellab.org/cms/?page=trimmomatic
#### phred33: specifies the base quality encoding
#### LEADING: Cut bases off the start of a read, if below a threshold quality
#### TRAILING: Cut bases off the end of a read, if below a threshold quality
#### MINLEN: Drop the read if it is below a specified length
#### threads: number of processors or threads to use

##### Option 1: 

```{shell}
java -classpath {SOFTWARE_DIR}/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 -threads 8 \
5176-MB-1-TCGGATTC-CGCAACTA_S01_L005_R1_001.fastq.gz 5176-MB-1-TCGGATTC-CGCAACTA_S01_L005_R2_001.fastq.gz \
5176-MB-1_R1.nodap.paired.txt 5176-MB-1_R1.noadap.unpaired.txt 5176-MB-1_R2.noadap.paired.txt  5176-MB-1_R2.noadap.unpaired.txt \
ILLUMINACLIP:TruSeq_CD_adapter.txt:2:30:7 LEADING:15 TRAILING:15 MINLEN:15 
```
##### Option 2: run in the background using nohup 
```{shell}
nohup java -classpath {SOFTWARE_DIR}/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 -threads 8 \
5176-MB-1-TCGGATTC-CGCAACTA_S01_L005_R1_001.fastq.gz 5176-MB-1-TCGGATTC-CGCAACTA_S01_L005_R2_001.fastq.gz \
5176-MB-1_R1.nodap.paired.txt 5176-MB-1_R1.noadap.unpaired.txt 5176-MB-1_R2.noadap.paired.txt  5176-MB-1_R2.noadap.unpaired.txt \
ILLUMINACLIP:TruSeq_CD_adapter.txt:2:30:7 LEADING:15 TRAILING:15 MINLEN:15 > 5176-trim-1.out &
```
##### Option 3: Refer to shell script loop example 

## Bowtie2
### manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
#### From CUT&RUN protocol.io For mapping spike-in fragments, we also use the --no-overlap --no-dovetail options to avoid cross-mapping of the experimental genome to that of the spike-in DNA.(https://www.protocols.io/view/cut-amp-run-targeted-in-situ-genome-wide-profiling-14egnr4ql5dy/v3?step=113)
#### local: Local alignment searches for the best alignment of a substring of the input sequence. While it can find an alignment for the entire sequence, if another, shorter, alignment has a higher score, it will be chosen. End-to-end will compute the score over the entire matching of the input sequence and its alignment with the reference. If there are adapters/long mismatches/indels etc. the local will work best. If you have a good reason to believe that the input sequence should be fully matched to the reference, then select end-to-end
#### very-sensitive-local: Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 
#### no-unal: Suppress SAM records for reads that failed to align.
#### no-mixed: By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.
#### no-discordant: A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints
#### no-overlap: If one mate alignment overlaps the other at all, consider that to be non-concordant
#### no-dovetail
#### I: The minimum fragment length for valid paired-end alignments. 
#### X: The maximum fragment length for valid paired-end alignments. 
#### p: threads or processors to use will running 

##### Option 1: 
```{shell} 
bowtie2 -p 8 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-overlap --no-dovetail -x {GENOME_DIR}/hg19_ec \
-1 4617-MB-1_R1.noadap.paired.txt -2 4617-MB-1_R2.noadap.paired.txt -S 4617-MB-1.hg19scer.sam 
```
##### Option 2: Run with nohup  

##### Option 3: Refer to shell script loop example 
