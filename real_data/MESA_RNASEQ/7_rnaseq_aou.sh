#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=4Gb
#SBATCH --array=701-4740
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 14218
# 14220 = 4740 * 3

source /home1/zeyunlu/init.sh
conda activate jax2

# SCRIPTF=/home1/zeyunlu/github

SCRIPTF=/project/nmancuso_8/zeyunlu/projects

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

EUR_COVAR=${SCRATCH}/mesa_rnaseq_covar_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_rnaseq_covar_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_rnaseq_covar_HIS.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_pred_rnaseq/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 3 *int(int($NR-1)))"`
stop=$((start + 2))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_rnaseq_cs_trait_list.tsv`
  set -- junk $params
  shift
  NAME=$2
  CHR=$4
  P0=$8
  P1=$9
  ID=${12}

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}

  # get genotype data
  bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis
  wfile=/scratch1/zeyunlu/sushie_rnaseq/sushie/weights/${ID}.normal.sushie.weights.tsv
  OUT=/project/nmancuso_8/data/sushie/aou_reg/

  if [ -f ${wfile} ]; then
    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

    for pop in EUR AFR HIS
    do
      echo "Getting geno and pheno data for ${pop}"

      KEEP=${SCRATCH}/mesa_rnaseq_pt_${pop}.tsv
      PHENO=${SCRATCH}/mesa_rnaseq_v1_${pop}.tsv.gz

      ${PLINK} --silent --bfile ${bfile} --extract $TMPDIR/${ID}.snp --rm-dup force-first \
        --make-bed --keep $KEEP --out $TMPDIR/${ID}.${pop}.geno

      python ${SCRIPTF}/sushie-data-codes/real_data/utils/extract_pheno.py \
      --file $PHENO \
      --gene $NAME \
      --subject $KEEP \
      --out $TMPDIR/${ID}.${pop}.pheno
    done

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict_aou.py \
      --pheno $TMPDIR/${ID}.EUR.pheno $TMPDIR/${ID}.AFR.pheno $TMPDIR/${ID}.HIS.pheno \
      --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
      --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
      --trait ${ID} \
      --w_files /scratch1/zeyunlu/sushie_rnaseq/ \
      --seed $((IDX + NR*50000)) \
      --study mesa.mrna \
      --out $OUT

    rm -rf ${TMPDIR}
  fi
done
rm -rf ${bigTMP}
