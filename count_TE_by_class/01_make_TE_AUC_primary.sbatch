#!/bin/bash
#SBATCH -c 4
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH --mem=16G
#SBATCH -o logs_primary/auc_%x_%j.out
#SBATCH -e logs_primary/auc_%x_%j.err

set -euo pipefail

# ============================================================
# One-sample PRIMARY-mapping TE AUC calculation
#
# Expected annotation:
#   Repeatmasker_IGV_magic_L1.bed
#
# BED format:
#   chr  start  end  RepeatMasker_label  score  strand
#
# Primary mapped alignments are defined as:
#
#   samtools view -F 2308
#
# This excludes:
#   4     unmapped
#   256   secondary
#   2048  supplementary
#
# Usage:
#   sbatch 01_make_TE_AUC_primary.sbatch \
#       SAMPLE.bam SAMPLE Repeatmasker_IGV_magic_L1.bed
# ============================================================

module purge
module load gcc/14.2.0
module load samtools/1.21
module load bedtools/2.31.0

THREADS=${SLURM_CPUS_PER_TASK:-4}

if [ "$#" -ne 3 ]; then
    echo "Usage:"
    echo "  sbatch $0 SAMPLE.bam SAMPLE Repeatmasker_IGV_magic_L1.bed"
    exit 1
fi

BAM=$1
SAMPLE=$2
REPEAT_BED=$3

WORKDIR="$(pwd)"

mkdir -p "${WORKDIR}/primary_bam"
mkdir -p "${WORKDIR}/normalized_bedgraph_primary"
mkdir -p "${WORKDIR}/repeat_AUC_primary"
mkdir -p "${WORKDIR}/stats_primary"
mkdir -p "${WORKDIR}/logs_primary"

PRIMARY_BAM="${WORKDIR}/primary_bam/${SAMPLE}.primary.bam"
BG="${WORKDIR}/normalized_bedgraph_primary/${SAMPLE}.normalized.bedGraph"
AUC_OUT="${WORKDIR}/repeat_AUC_primary/${SAMPLE}.repeat_AUC.tsv"
COUNT_OUT="${WORKDIR}/stats_primary/${SAMPLE}.primary_count.tsv"
CLEAN_BED="${WORKDIR}/stats_primary/${SAMPLE}.Repeatmasker_IGV_magic_L1.cleaned.bed"

if [ ! -f "${BAM}" ]; then
    echo "ERROR: BAM not found: ${BAM}"
    exit 1
fi

if [ ! -f "${REPEAT_BED}" ]; then
    echo "ERROR: repeat BED not found: ${REPEAT_BED}"
    exit 1
fi

# Clean possible CRLF line endings.
awk 'BEGIN{OFS="\t"} { sub(/\r$/, "", $0); print }' \
    "${REPEAT_BED}" \
    > "${CLEAN_BED}"

# Validate RepeatMasker BED.
awk '
BEGIN { bad=0 }

NF < 6 {
    print "ERROR: RepeatMasker BED must have at least 6 columns." > "/dev/stderr"
    bad=1
    exit
}

$6 != "+" && $6 != "-" {
    print "ERROR: BED column 6 must be + or -." > "/dev/stderr"
    print "Problem line: " $0 > "/dev/stderr"
    print "Parsed strand field: [" $6 "]" > "/dev/stderr"
    bad=1
    exit
}

END {
    if (bad) exit 1
}
' "${CLEAN_BED}"


# ============================================================
# Count primary mapped alignments
# ============================================================

N_PRIMARY=$(samtools view \
    -@ "${THREADS}" \
    -c \
    -F 2308 \
    "${BAM}"
)

if [ "${N_PRIMARY}" -eq 0 ]; then
    echo "ERROR: ${SAMPLE} has zero primary mapped alignments."
    exit 1
fi

echo -e "sample\tprimary_mapped_alignments" > "${COUNT_OUT}"
echo -e "${SAMPLE}\t${N_PRIMARY}" >> "${COUNT_OUT}"

echo "============================================================"
echo "Sample: ${SAMPLE}"
echo "Primary mapped alignments: ${N_PRIMARY}"
echo "============================================================"


# ============================================================
# Primary-only BAM
# ============================================================

samtools view \
    -@ "${THREADS}" \
    -bh \
    -F 2308 \
    "${BAM}" \
    > "${PRIMARY_BAM}"

samtools index \
    -@ "${THREADS}" \
    "${PRIMARY_BAM}"


# ============================================================
# Coverage
# ============================================================

bedtools genomecov \
    -ibam "${PRIMARY_BAM}" \
    -bg \
| awk \
    -v N="${N_PRIMARY}" \
    'BEGIN { OFS="\t" }
     {
         print $1,$2,$3,$4*1000000/N
     }' \
> "${BG}"


# ============================================================
# Parse RepeatMasker annotation and calculate TE-region AUC
# ============================================================

bedtools intersect \
    -a <(
        awk '
        BEGIN {
            OFS="\t"
        }

        function first_piece(s,    a) {
            split(s,a,"_")
            return a[1]
        }

        {
            label=$4
            strand=$6

            if (label ~ /^Simple_repeat_/) {

                te_class="Simple_repeat"
                rest=substr(label,15)
                te_name=first_piece(rest)

            } else if (label ~ /^Low_complexity_/) {

                te_class="Low_complexity"
                rest=substr(label,16)
                te_name=first_piece(rest)

            } else if (label ~ /^Unknown_/) {

                te_class="Unknown"
                rest=substr(label,9)
                te_name=first_piece(rest)

            } else {

                n=split(label,a,"_")

                class_family=a[1]

                split(class_family,b,"/")
                te_class=b[1]

                if (n >= 2) {
                    te_name=a[2]
                } else {
                    te_name=class_family
                }
            }

            print $1,$2,$3,te_name,te_class,strand,label
        }
        ' "${CLEAN_BED}"
    ) \
    -b "${BG}" \
    -wao \
| awk '
    BEGIN {
        OFS="\t"

        print \
            "chr", \
            "start", \
            "end", \
            "TE_name", \
            "TE_class", \
            "TE_strand", \
            "RepeatMasker_label", \
            "AUC"
    }

    function output_previous() {

        if (have_region) {

            print \
                prev_chr, \
                prev_start, \
                prev_end, \
                prev_name, \
                prev_class, \
                prev_strand, \
                prev_label, \
                auc
        }
    }

    {
        chr=$1
        start=$2
        end=$3
        name=$4
        class=$5
        strand=$6
        label=$7

        key=chr FS start FS end FS name FS class FS strand FS label

        if (!have_region) {

            previous_key=key

            prev_chr=chr
            prev_start=start
            prev_end=end
            prev_name=name
            prev_class=class
            prev_strand=strand
            prev_label=label

            auc=0
            have_region=1
        }

        if (key != previous_key) {

            output_previous()

            previous_key=key

            prev_chr=chr
            prev_start=start
            prev_end=end
            prev_name=name
            prev_class=class
            prev_strand=strand
            prev_label=label

            auc=0
        }

        # A = 7 columns
        # bedGraph = 4 columns
        # $11 = coverage
        # $12 = overlap length
        if ($12 > 0 && $11 != ".") {
            auc += $11 * $12
        }
    }

    END {
        output_previous()
    }
' > "${AUC_OUT}"


echo
echo "Finished ${SAMPLE}"
echo "AUC table:"
echo "  ${AUC_OUT}"
