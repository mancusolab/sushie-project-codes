#!/bin/bash

set -o errexit
set -o nounset

chmod +x ${PLINK2}

# Let's see what variables are really set - is everything set the way we
# thought?
env | sort

NGENES_METADATA=`wc -l ${METADATA_FILE} | awk '{print $1}'`

for idx in `seq 1 $NGENES_METADATA`
do
    # parse the $GENES file
    # 1       ENSG00000000971_SL000415        mesa.proteins

    params=`sed "${idx}q;d" ${METADATA_FILE}`
    set -- junk $params
    shift
    CHR=${1}
    STUDY=${2}
    INPUT_GENE=${3}     

    echo "analyzing ${INPUT_GENE} in ${STUDY}"
    
    PREFIX=${STUDY}.chr${INPUT_CHR}.${INPUT_GENE}
    INPUT_WEIGHTS=${SUSHIE_WEIGHT_INPUT_FOLDER}/${PREFIX}.weights.tsv
    INPUT_SNPS=${SUSHIE_SNP_INPUT_FOLDER}/${PREFIX}.snps.tsv
    OUTPUT_SCORE=${SUSHIE_WEIGHT_OUTPUT_FOLDER}/${PREFIX}.ge

    echo "input_weights: ${INPUT_WEIGHTS}"
    echo "input_snps: ${INPUT_SNPS}"
    echo "output_score ${OUTPUT_SCORE}"
    echo "Running ${PLINK2} $(date)"
    

    ${PLINK2} \
        --bed ${BEDFILE} \
        --bim ${BIMFILE} \
        --fam ${FAMFILE} \
        --extract ${INPUT_SNPS} \
        --memory 12000 \
        --score ${INPUT_WEIGHTS} 1 2 header-read \
        --score-col-nums 3-18 \
        --out ${OUTPUT_SCORE}

    echo "gzipping $(date)"
    gzip ${OUTPUT_SCORE}.sscore

done
