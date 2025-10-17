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

# Default number of reads
number_of_reads=5

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
  HAP_BAM="${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.bam"
  OUT_FILE="${RESULTS_DIR}/spanning_readnames_${HAP}.txt"

  if [ ! -s "$HAP_BAM" ]; then
    log "Warning: ${HAP_BAM} missing or empty. Writing empty file."
    : > "$OUT_FILE"
    continue
  fi

  READ_COUNT=$(samtools view -c "$HAP_BAM" 2>/dev/null || echo 0)

  if [ "$READ_COUNT" -eq 0 ]; then
    log "No reads found in ${HAP_BAM}. Writing empty file."
    : > "$OUT_FILE"
    continue
  fi

  if ! bedtools intersect -f 1.0 -wb \
    -a "${INPUT_DIR}/prdm9_target.bed" \
    -b "$HAP_BAM" \
    | cut -f7 | head -n "$number_of_reads" > "$OUT_FILE"; then
      log "bedtools failed for ${HAP}. Writing empty file."
      : > "$OUT_FILE"
      continue
  fi

  log "${HAP}: wrote $(wc -l < "$OUT_FILE") read names."
done

# Check haplotype read counts and determine mode
HP1_COUNT=$(wc -l < "${RESULTS_DIR}/spanning_readnames_HP1.txt")
HP2_COUNT=$(wc -l < "${RESULTS_DIR}/spanning_readnames_HP2.txt")

log "HP1 has ${HP1_COUNT} spanning reads."
log "HP2 has ${HP2_COUNT} spanning reads."

if [ "$HP1_COUNT" -lt "$number_of_reads" ] || [ "$HP2_COUNT" -lt "$number_of_reads" ]; then
  MODE="unphased"
  log "Detected homozygous region (one or both HP < $number_of_reads reads). Using unphased reads."

  # Create unphased spanning read list from all reads in the region
  bedtools intersect -f 1.0 -wb \
    -a "${INPUT_DIR}/prdm9_target.bed" \
    -b "${RESULTS_DIR}/${PROBAND}.PRDM9.bam" \
    | cut -f7 | sort -u | head -n "$number_of_reads" > "${RESULTS_DIR}/spanning_readnames_unphased.txt"

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

## Mapping reads and generating match files
log "Mapping reads and generating match files..."
for HAP in HP1 HP2; do
  log "Processing ${HAP}..."
  samtools fasta "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.bam" > "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.fasta"
  grep -f "${RESULTS_DIR}/spanning_readnames_${HAP}.txt" "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9.fasta" -A 1 | grep -v "\-\-" > "${RESULTS_DIR}/targetreads_${HAP}.fasta"

  for i in $(seq 1 $number_of_reads); do
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
  FILES=(${RESULTS_DIR}/read*.${HAP}.prdm9matches.tsv)
  if [ ${#FILES[@]} -eq 0 ]; then
    log "No match files found for ${HAP}, skipping."
    continue
  fi

  paste "${FILES[@]}" | awk 'BEGIN{FS=OFS="\t"}{
      # Print first column once (haplotype ID)
      printf "%s", $1
      for(i=2; i<=NF; i+=3){
          # Print YES/NO and mismatch only, skip repeated hap name
          printf "\t%s\t%s", $(i), $(i+1)
      }
      print ""
  }' > "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9_combined.tsv"

done

# Summarize best haplotype(s) per haplotype and add phasing info
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
    if(yes_count==n){   # Only haplotypes with all YES
      med=median(mis_arr, n)
      mismatches=mis_arr[1]
      for(j=2;j<=n;j++) mismatches=mismatches","mis_arr[j]
      if(min_med=="" || med<min_med){
        min_med=med
        lines=SAMPLE "\t" HAP "\t" $1 "\t" yes_count "\t" med "\t" mismatches "\t" MODE "\t" HP1_COUNT "\t" HP2_COUNT
      } else if(med==min_med){
        lines=lines "\n"SAMPLE "\t" HAP "\t" $1 "\t" yes_count "\t" med "\t" mismatches "\t" MODE "\t" HP1_COUNT "\t" HP2_COUNT
      }
    }
  }
  END {
    if(lines!="") print lines
  }' "${RESULTS_DIR}/${PROBAND}.${HAP}.PRDM9_combined.tsv" > "$HAPLO_MATCH"

  log "Haplotypematch written to $HAPLO_MATCH with mode info."
done

log "Pipeline finished for ${PROBAND}."