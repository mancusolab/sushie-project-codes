#!/bin/bash

set -o errexit
set -o nounset

# Let's see what variables are really set - is everything set the way we
# thought?
env | sort

NGENES_METADATA=`wc -l ${INPUT_FILE_METADATA} | awk '{print $1}'`

for idx in `seq 1 $NGENES_METADATA`
do
    # parse the $GENES file
    # 1       ENSG00000000971_SL000415        mesa.proteins
    params=`sed "${idx}q;d" ${INPUT_FILE_METADATA}`
    set -- junk $params
    shift
    # get INPUT_STUDY from the --env flag
    CHROM=${1}
    INPUT_GENE=${3}     # ENSG00000000938_SL006912
    STUDY=${2}    # mesa.proteins

    # echo "analyzing ${INPUT_GENE}, ${INPUT_STUDY}"
    PREFIX=${STUDY}.chr${CHROM}.${INPUT_GENE}
    INPUT_FILE_SUSHIE_WEIGHT=${INPUT_FOLDER_SUSHIE_WEIGHT}/chr${CHROM}/${PREFIX}.ge.sscore.gz

    # echo "INPUT_FILE_SUSHIE_WEIGHT ${INPUT_FILE_SUSHIE_WEIGHT}"
    # echo "Running Rscript $(date)"

    # echo "INPUT_FILE_ANCESTRY ${INPUT_FILE_ANCESTRY}"

    Rscript ${INPUT_RSCRIPT} \
            ${CHROM} \
            ${INPUT_GENE} \
            ${STUDY} \
            ${INPUT_FILE_SUSHIE_WEIGHT} \
            ${INPUT_FOLDER_PHENO_LIST}/ \
            ${INPUT_FOLDER_ANCESTRY_LIST}/ancestry_preds.tsv \
            ${EXCLUSION_PT}/exclusion.criteria.tsv \
            ${OUTPUT_FOLDER_SUSHIE_TWAS}/${PREFIX}.twas.tsv

    gzip ${OUTPUT_FOLDER_SUSHIE_TWAS}/${PREFIX}.twas.tsv
done

