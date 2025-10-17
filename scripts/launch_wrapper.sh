#!/bin/bash

# === Config ===
BASE_DIR="/ifs/data/research/projects/michelle/17_prdm9_haplotypecheck"
TRIO_FILE="${BASE_DIR}/data/trios_bam_vcf.tsv"
PIPELINE_SCRIPT="${BASE_DIR}/scripts/prdm9_pipeline.sh"
RESULTS_BASE="${BASE_DIR}/results"

MAX_SAMPLES=0  # 0 = all
RUN_NAME="use_5_reads"

# === Slurm submission ===
count=0
while IFS=$'\t' read -r child father mother child_bam father_bam mother_bam vcf_raw; do
    ((count++))
    if [[ $MAX_SAMPLES -gt 0 && $count -gt $MAX_SAMPLES ]]; then
        echo "[INFO] Reached MAX_SAMPLES=$MAX_SAMPLES. Stopping."
        break
    fi

    for SAMPLE_NAME in "$child" "$father" "$mother"; do
        BAM_FILE=""
        case "$SAMPLE_NAME" in
            "$child") BAM_FILE="$child_bam" ;;
            "$father") BAM_FILE="$father_bam" ;;
            "$mother") BAM_FILE="$mother_bam" ;;
        esac

        echo "[INFO] Submitting sample: $SAMPLE_NAME"
        SAMPLE_RESULTS="${RESULTS_BASE}/${RUN_NAME}/${SAMPLE_NAME}"
        mkdir -p "$SAMPLE_RESULTS"
        LOG_FILE="${SAMPLE_RESULTS}/prdm9_${SAMPLE_NAME}.log"

        sbatch --job-name="prdm9_${SAMPLE_NAME}" \
               --output="${LOG_FILE}" \
               --ntasks=1 \
               --cpus-per-task=1 \
               --mem-per-cpu=16100M \
               --partition=develop \
               --wrap="bash $PIPELINE_SCRIPT ${SAMPLE_NAME} ${BAM_FILE} ${RUN_NAME}"
    done
done < "$TRIO_FILE"

# === Local execution (optional) ===
# echo "[INFO] Starting prdm9 analysis locally..."
# count=0
# while IFS=$'\t' read -r child father mother child_bam father_bam mother_bam vcf_raw; do
#     ((count++))
#     if [[ $MAX_SAMPLES -gt 0 && $count -gt $MAX_SAMPLES ]]; then
#         echo "[INFO] Reached MAX_SAMPLES=$MAX_SAMPLES. Stopping local execution."
#         break
#     fi

#     for SAMPLE_NAME in "$child" "$father" "$mother"; do
#         BAM_FILE=""
#         case "$SAMPLE_NAME" in
#             "$child") BAM_FILE="$child_bam" ;;
#             "$father") BAM_FILE="$father_bam" ;;
#             "$mother") BAM_FILE="$mother_bam" ;;
#         esac

#         echo "[INFO] Running in terminal for sample: $SAMPLE_NAME"
#         SAMPLE_RESULTS="${RESULTS_BASE}/${SAMPLE_NAME}"
#         mkdir -p "$SAMPLE_RESULTS"
#         LOG_FILE="${SAMPLE_RESULTS}/prdm9_${SAMPLE_NAME}_local.log"

#         bash "$PIPELINE_SCRIPT" "$SAMPLE_NAME" "$BAM_FILE" > "$LOG_FILE" 2>&1
#         echo "[INFO] Finished $SAMPLE_NAME. Log: $LOG_FILE"
#     done
# done < "$TRIO_FILE"