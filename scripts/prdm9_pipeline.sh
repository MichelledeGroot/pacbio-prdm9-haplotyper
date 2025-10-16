#!/usr/bin/env bash
set -euo pipefail

# Usage check
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <sample_id> <bam_file> <run_name>"
  exit 1
fi

PROBAND="$1"
PROBAND_BAM="$2"
RUN="$3"

# Base folder setup
BASE_DIR="/ifs/data/research/projects/michelle/17_prdm9_haplotypecheck"
INPUT_DIR="${BASE_DIR}/data"
SCRIPTS_DIR="${BASE_DIR}/scripts"
RESULTS_DIR="${BASE_DIR}/results/${RUN}/${PROBAND}"
mkdir -p "$INPUT_DIR" "$SCRIPTS_DIR" "$RESULTS_DIR"

# Log file
LOG_FILE="${RESULTS_DIR}/${PROBAND}_pipeline.log"

# Logger function
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

log "Processing ${PROBAND}..."
log "BAM: ${PROBAND_BAM}"

# load modules
module load bioinf/samtools
module load bioinf/minimap2
module load bioinf/bedtools

# PRDM9 region
REGION="chr5:23526573-23527864"

log "Extracting PRDM9 region..."
samtools view -b -h "${PROBAND_BAM}" "${REGION}" > "${RESULTS_DIR}/${PROBAND}.PRDM9.bam"

log "Splitting into haplotypes..."
samtools view -d HP:1 -o "${RESULTS_DIR}/${PROBAND}.HP1.PRDM9.bam" "${RESULTS_DIR}/${PROBAND}.PRDM9.bam"
samtools view -d HP:2 -o "${RESULTS_DIR}/${PROBAND}.HP2.PRDM9.bam" "${RESULTS_DIR}/${PROBAND}.PRDM9.bam"

log "Finding spanning reads..."
for HAP in HP1 HP2; do
  bedtools intersect -f 1.0 -wb \
    -a "${INPUT_DIR}/prdm9_target.bed" \
    -b "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.bam" \
    | cut -f7 | head -n 5 > "${RESULTS_DIR}/spanning_readnames_${HAP}.txt"
done

# Check haplotype read counts and determine mode
HP1_COUNT=$(wc -l < "${RESULTS_DIR}/spanning_readnames_HP1.txt")
HP2_COUNT=$(wc -l < "${RESULTS_DIR}/spanning_readnames_HP2.txt")

log "HP1 has ${HP1_COUNT} spanning reads."
log "HP2 has ${HP2_COUNT} spanning reads."

if [ "$HP1_COUNT" -lt 5 ] || [ "$HP2_COUNT" -lt 5 ]; then
  MODE="unphased"
  log "Detected homozygous region (one or both HP < 5 reads). Using unphased reads."

  # Create unphased spanning read list from all reads in the region
  bedtools intersect -f 1.0 -wb \
    -a "${INPUT_DIR}/prdm9_target.bed" \
    -b "${RESULTS_DIR}/${PROBAND}.PRDM9.bam" \
    | cut -f7 | sort -u | head -n 5 > "${RESULTS_DIR}/spanning_readnames_unphased.txt"

  # Duplicate unphased read list for HP1 and HP2
  cp "${RESULTS_DIR}/spanning_readnames_unphased.txt" "${RESULTS_DIR}/spanning_readnames_HP1.txt"
  cp "${RESULTS_DIR}/spanning_readnames_unphased.txt" "${RESULTS_DIR}/spanning_readnames_HP2.txt"

  # Duplicate BAMs so both HPs use the same reads
  cp "${RESULTS_DIR}/${PROBAND}.PRDM9.bam" "${RESULTS_DIR}/${PROBAND}.HP1.PRDM9.bam"
  cp "${RESULTS_DIR}/${PROBAND}.PRDM9.bam" "${RESULTS_DIR}/${PROBAND}.HP2.PRDM9.bam"

  log "Set HP1 and HP2 to identical unphased reads."
  log "Mode: ${MODE}"
  log "Unphased reads used:"
  awk '{print "  " $0}' "${RESULTS_DIR}/spanning_readnames_unphased.txt" | tee -a "$LOG_FILE"

else
  MODE="phased"
  log "Mode: ${MODE}"
  log "HP1 reads used:"
  awk '{print "  " $0}' "${RESULTS_DIR}/spanning_readnames_HP1.txt" | tee -a "$LOG_FILE"
  log "HP2 reads used:"
  awk '{print "  " $0}' "${RESULTS_DIR}/spanning_readnames_HP2.txt" | tee -a "$LOG_FILE"
fi

## Get fasta and find best haplotype match
log "Mapping reads and generating match files..."
for HAP in HP1 HP2; do
  log "Processing ${HAP}..."
  samtools fasta "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.bam" > "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.fasta"
  grep -f "${RESULTS_DIR}/spanning_readnames_${HAP}.txt" "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.fasta" -A 1 | grep -v "\-\-" > "${RESULTS_DIR}/targetreads_${HAP}.fasta"

  for i in {1..5}; do
    head -n $((i * 2)) "${RESULTS_DIR}/targetreads_${HAP}.fasta" | tail -n 2 > "${RESULTS_DIR}/read${i}_${HAP}.fasta"

    minimap2 -c -x map-hifi --seed 123 --secondary=no --paf-no-hit \
      "${RESULTS_DIR}/read${i}_${HAP}.fasta" "${INPUT_DIR}/prdm9_ref.fa" \
      > "${RESULTS_DIR}/read${i}.${HAP}_prdm9motifs.paf"

    awk '!seen[$1]++' "${RESULTS_DIR}/read${i}.${HAP}_prdm9motifs.paf" \
      > "${RESULTS_DIR}/read${i}.${HAP}_prdm9motifs_unique.paf"

    # ±10 bp tolerance for YES/NO
    awk 'BEGIN {FS=OFS="\t"} {
      diff = ($2 > $11) ? ($2 - $11) : ($11 - $2);
      if (diff <= 10) print $1, "YES", $13;
      else print $1, "NO", $13;
    }' "${RESULTS_DIR}/read${i}.${HAP}_prdm9motifs_unique.paf" \
      | sed 's/NM:i://g' | sort \
      > "${RESULTS_DIR}/read${i}.${HAP}.prdm9matches.tsv"
  done
done

log "Combining match files per haplotype..."
for HAP in HP1 HP2; do
  FIRST_FILE=$(ls "${RESULTS_DIR}"/read*.${HAP}.prdm9matches.tsv | head -n1)
  cut -f1 "$FIRST_FILE" > "${RESULTS_DIR}/tmp_col1.txt"
  paste "${RESULTS_DIR}/tmp_col1.txt" "${RESULTS_DIR}"/read*.${HAP}.prdm9matches.tsv | cut -f1,3,4,6,7,9,10,12,13,15,16 \
    > "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9_combined.tsv"
  rm -f "${RESULTS_DIR}/tmp_col1.txt"
done

log "Summarizing best haplotype(s) per haplotype and adding phasing info..."
for HAP in HP1 HP2; do
  HAPLO_MATCH="${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9_haplotypematch.tsv"
  awk -v SAMPLE="$PROBAND" -v HAP="$HAP" -v MODE="$MODE" -v HP1_COUNT="$HP1_COUNT" -v HP2_COUNT="$HP2_COUNT" '
  function median(a, n) {
    asort(a)
    if (n % 2) return a[(n+1)/2]
    else return (a[n/2] + a[n/2+1])/2
  }
  {
    yes_count=0
    delete mis_arr
    n=0
    for(i=2;i<=NF;i+=2){
      if($i=="YES") yes_count++
      n++; mis_arr[n]=$(i+1)
    }
    if(yes_count==5){
      med=median(mis_arr, n)
      mismatches=mis_arr[1]
      for(j=2;j<=n;j++) mismatches=mismatches","mis_arr[j]
      if(min_med=="" || med<min_med){
        min_med=med
        lines=SAMPLE "\t" HAP "\t" $1 "\t" yes_count "\t" med "\t" mismatches "\t" MODE "\t" HP1_COUNT "\t" HP2_COUNT
      } else if(med==min_med){
        lines=lines "\n" SAMPLE "\t" HAP "\t" $1 "\t" yes_count "\t" med "\t" mismatches "\t" MODE "\t" HP1_COUNT "\t" HP2_COUNT
      }
    }
  }
  END {
    if(lines!="") print "sample\thap\tbest_haplotype\t#YES\tmedian_mismatches\tmismatches_per_read\tmode\tHP1_reads\tHP2_reads"
    if(lines!="") print lines
  }' "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9_combined.tsv" > "$HAPLO_MATCH"

  log "Haplotypematch written to $HAPLO_MATCH with mode info."
done

log "Pipeline finished for ${PROBAND}."