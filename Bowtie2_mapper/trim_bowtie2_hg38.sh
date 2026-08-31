#!/bin/bash
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH --mem=40G
#SBATCH -o trim_bowtie2_%j.out
#SBATCH -e trim_bowtie2_%j.err

set -euo pipefail

# ============================================================
# Trimmomatic -> Bowtie2 hg38
#
# Usage:
#   sbatch trim_bowtie2_hg38.sh SAMPLE R1.fastq.gz R2.fastq.gz
#
# Edit HG38_INDEX and ADAPTERS below before running.
# ============================================================

module purge
module load gcc/14.2.0
module load trimmomatic/0.39
module load bowtie2/2.5.4
module load samtools/1.21

THREADS=${SLURM_CPUS_PER_TASK:-8}

# Bowtie2 index prefix, NOT a fasta filename.
HG38_INDEX="/path/to/Human_Gchr/hg38"

# Use the adapter FASTA appropriate for the library.
ADAPTERS="/path/to/adapters.fa"

if [ "$#" -ne 3 ]; then
    echo "Usage:"
    echo "sbatch $0 SAMPLE R1.fastq.gz R2.fastq.gz"
    exit 1
fi

SAMPLE=$1
R1=$2
R2=$3

if [ ! -f "${R1}" ]; then
    echo "ERROR: R1 not found: ${R1}"
    exit 1
fi

if [ ! -f "${R2}" ]; then
    echo "ERROR: R2 not found: ${R2}"
    exit 1
fi

if [ ! -f "${ADAPTERS}" ]; then
    echo "ERROR: adapter FASTA not found: ${ADAPTERS}"
    exit 1
fi

OUTDIR="${SAMPLE}_Hg38_BT2"
mkdir -p "${OUTDIR}/trimmed"

R1_PAIRED="${OUTDIR}/trimmed/${SAMPLE}_R1_paired.fastq.gz"
R1_UNPAIRED="${OUTDIR}/trimmed/${SAMPLE}_R1_unpaired.fastq.gz"
R2_PAIRED="${OUTDIR}/trimmed/${SAMPLE}_R2_paired.fastq.gz"
R2_UNPAIRED="${OUTDIR}/trimmed/${SAMPLE}_R2_unpaired.fastq.gz"

BAM="${OUTDIR}/${SAMPLE}_All_wholeHg38_BT2.bam"

echo "============================================================"
echo "Trimming ${SAMPLE}"
echo "============================================================"

trimmomatic PE \
    -threads "${THREADS}" \
    -phred33 \
    "${R1}" \
    "${R2}" \
    "${R1_PAIRED}" \
    "${R1_UNPAIRED}" \
    "${R2_PAIRED}" \
    "${R2_UNPAIRED}" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

echo
echo "============================================================"
echo "Bowtie2 mapping ${SAMPLE} to hg38"
echo "============================================================"

bowtie2 \
    --very-sensitive \
    -p "${THREADS}" \
    -x "${HG38_INDEX}" \
    -1 "${R1_PAIRED}" \
    -2 "${R2_PAIRED}" \
    2> "${OUTDIR}/${SAMPLE}.bowtie2.log" \
| samtools view \
    -@ "${THREADS}" \
    -b \
    - \
| samtools sort \
    -@ "${THREADS}" \
    -o "${BAM}" \
    -

samtools index \
    -@ "${THREADS}" \
    "${BAM}"

samtools flagstat \
    -@ "${THREADS}" \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.flagstat.txt"

samtools idxstats \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.idxstats.txt"

echo
echo "============================================================"
echo "DONE"
echo "============================================================"
echo "BAM:"
echo "  ${BAM}"
echo "Index:"
echo "  ${BAM}.bai"
echo "Bowtie2 log:"
echo "  ${OUTDIR}/${SAMPLE}.bowtie2.log"
echo "Flagstat:"
echo "  ${OUTDIR}/${SAMPLE}.flagstat.txt"
