#!/bin/bash
set -o errexit
set -o xtrace

export JOB_NAME="score"
export INPUT_CHR=${1}

export BILLING_PROJECT="${GOOGLE_PROJECT}"
export PROJECT_BUCKET="${WORKSPACE_BUCKET}"
export SHARED_PROJECT_PATH="${PROJECT_BUCKET}/zeyun/aou"
export SUSHIE_TMP_FOLDER="${PROJECT_BUCKET}/zeyun/sushie_tmp"
export SUSHIE_TMP_SUBFOLDER="${SUSHIE_TMP_FOLDER}"
export SUSHIE_WEIGHT_INPUT_FOLDER="${SHARED_PROJECT_PATH}/weights/chr${INPUT_CHR}/"
export SUSHIE_SNP_INPUT_FOLDER="${SHARED_PROJECT_PATH}/snps/chr${INPUT_CHR}"
export SUSHIE_WEIGHT_OUTPUT_FOLDER="${SHARED_PROJECT_PATH}/outputs/chr${INPUT_CHR}"

export METADATA_FILE="${SHARED_PROJECT_PATH}/metadata/all.chr${INPUT_CHR}.metadata.tsv"

export AOU_NETWORK=network
export AOU_SUBNETWORK=subnetwork
export DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"
export MAX_PREEMPTION=4

PLINK2="gs://fc-secure-b522bed9-b978-46f7-ad16-99da7c1f5490/bin/plink2"
SCRIPT="gs://fc-secure-b522bed9-b978-46f7-ad16-99da7c1f5490/zeyun/scripts/score_run.sh"
BFILE_PREFIX="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed"

# Localize your tasks file.
# gsutil cp ${TASKS_FILE} ./
# Note: we are not using the tasks file, but instead using the metadata file

# input list
# mesa.proteins.chr22.ENSG00000015475_SL003704.sushie.weights.tsv
# mesa.proteins.chr22.ENSG00000015475_SL003704.sushie.snps.tsv

dsub \
    --provider google-cls-v2 \
    --regions us-central1 \
    --user-project "${GOOGLE_PROJECT}" \
    --project "${GOOGLE_PROJECT}" \
    --image 'marketplace.gcr.io/google/ubuntu2004:latest' \
    --network "${AOU_NETWORK}" \
    --subnetwork "${AOU_SUBNETWORK}" \
    --service-account "$(gcloud config get-value account)" \
    --user "${DSUB_USER_NAME}" \
    --disk-size 512 \
    --boot-disk-size 20 \
    --min-ram 16 \
    --preemptible \
    --input BEDFILE=${BFILE_PREFIX}/acaf_threshold.chr${INPUT_CHR}.bed \
    --input BIMFILE=${BFILE_PREFIX}/acaf_threshold.chr${INPUT_CHR}.bim \
    --input FAMFILE=${BFILE_PREFIX}/acaf_threshold.chr${INPUT_CHR}.fam \
    --input METADATA_FILE=${METADATA_FILE} \
    --input PLINK2=$PLINK2 \
    --input-recursive SUSHIE_WEIGHT_INPUT_FOLDER=${SUSHIE_WEIGHT_INPUT_FOLDER} \
    --input-recursive SUSHIE_SNP_INPUT_FOLDER=${SUSHIE_SNP_INPUT_FOLDER} \
    --output-recursive SUSHIE_WEIGHT_OUTPUT_FOLDER=${SUSHIE_WEIGHT_OUTPUT_FOLDER} \
    --env SUSHIE_TMP_FOLDER=${SUSHIE_TMP_FOLDER} \
    --env INPUT_CHR=${INPUT_CHR} \
    --logging "${SUSHIE_TMP_FOLDER}" \
    --name "TEST.${JOB_NAME}.${INPUT_CHR}" \
    --script "${SCRIPT}"

# # we have 2 input recursive folders (weights and snps), im worried that this will cause problems, so I added a line in score.sh to extract the snps.tsv file
