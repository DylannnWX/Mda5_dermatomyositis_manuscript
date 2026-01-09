#!/bin/bash
#SBATCH -c 6                               # Request one core
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)


set -e

Group=$1

module load gcc/14.2.0 samtools/1.21 bedtools/2.31.0 trimmomatic/0.39 bowtie2/2.5.4 bcftools/1.21 python/3.13.1 star/2.7.11b

java -jar $TRIMMOMATIC/trimmomatic-0.39.jar PE -threads 6 -phred33 ${Group}*R1_001.fastq.gz ${Group}*R2_001.fastq.gz ${Group}_R1_paired.fq ${Group}_R1_unpaired.fq ${Group}_R2_paired.fq ${Group}_R2_unpaired.fq ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:10
rm ${Group}_R1_unpaired.fq ${Group}_R2_unpaired.fq
STAR --runThreadN 6 --genomeDir your/STAR/reference/DIR/human --readFilesIn ${Group}_R1_paired.fq ${Group}_R2_paired.fq --outFileNamePrefix ${Group}_Hg38STAR --outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outSAMunmapped Within --outSAMattributes NH HI NM MD AS --outReadsUnmapped Fastx
samtools index ${Group}_Hg38STARAligned.sortedByCoord.out.bam


