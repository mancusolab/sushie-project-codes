#!/bin/bash
set -o errexit
set -o xtrace

export JOB_NAME="assocTest"
export AOU_NETWORK=network
export AOU_SUBNETWORK=subnetwork
export DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"
export MAX_PREEMPTION=4

# export INPUT_CHR=$1
export BILLING_PROJECT="${GOOGLE_PROJECT}"
export PROJECT_BUCKET="${WORKSPACE_BUCKET}"
export SHARED_PROJECT_PATH="${PROJECT_BUCKET}/zeyun/aou"

export INPUT_FILE_METADATA="${SHARED_PROJECT_PATH}/metadata_twas2/twas.group${1}.meta.tsv"
export INPUT_FOLDER_SUSHIE_WEIGHT="${SHARED_PROJECT_PATH}/outputs/"
export INPUT_FOLDER_ANCESTRY_LIST="${PROJECT_BUCKET}/data/aou/ancestry_pred/"
export EXCLUSION_PT_FILE="${PROJECT_BUCKET}/data/aou/pheno/20240927/"
export INPUT_FOLDER_PHENO_LIST="${PROJECT_BUCKET}/data/aou/pheno/20240927/"

export OUTPUT_FOLDER_SUSHIE_TWAS="${SHARED_PROJECT_PATH}/twas_outputs/"
export LOGGING_FOLDER_SUSHIE_TWAS="${SHARED_PROJECT_PATH}/twas_tmp"

INPUT_SCRIPT="gs://fc-secure-b522bed9-b978-46f7-ad16-99da7c1f5490/zeyun/scripts/twas_run.sh"
INPUT_RSCRIPT="gs://fc-secure-b522bed9-b978-46f7-ad16-99da7c1f5490/zeyun/scripts/twas.R"

dsub \
    --provider google-cls-v2 \
    --regions us-central1 \
    --user-project "${GOOGLE_PROJECT}" \
    --project "${GOOGLE_PROJECT}" \
    --image rocker/tidyverse:latest \
    --network "${AOU_NETWORK}" \
    --subnetwork "${AOU_SUBNETWORK}" \
    --service-account "$(gcloud config get-value account)" \
    --user "${DSUB_USER_NAME}" \
    --disk-size 512 \
    --boot-disk-size 20 \
    --min-ram 16 \
    --preemptible \
    --input INPUT_RSCRIPT=${INPUT_RSCRIPT} \
    --input INPUT_FILE_METADATA=${INPUT_FILE_METADATA} \
    --input-recursive INPUT_FOLDER_ANCESTRY_LIST=${INPUT_FOLDER_ANCESTRY_LIST} \
    --input-recursive INPUT_FOLDER_PHENO_LIST=${INPUT_FOLDER_PHENO_LIST} \
    --input-recursive INPUT_FOLDER_SUSHIE_WEIGHT=${INPUT_FOLDER_SUSHIE_WEIGHT} \
    --input-recursive EXCLUSION_PT=${EXCLUSION_PT_FILE} \
    --output-recursive OUTPUT_FOLDER_SUSHIE_TWAS=${OUTPUT_FOLDER_SUSHIE_TWAS} \
    --env INPUT_CHR=${INPUT_CHR} \
    --logging "${LOGGING_FOLDER_SUSHIE_TWAS}" \
    --name "group.${1}" \
    --script "${INPUT_SCRIPT}"
