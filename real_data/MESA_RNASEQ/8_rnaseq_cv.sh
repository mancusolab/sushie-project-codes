#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=48Gb
#SBATCH --array=2320
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 21855
# 21855 = 4371 * 5

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

# SCRIPTF=/home1/zeyunlu/github
OUT=/scratch1/zeyunlu/sushie_rnaseq
SCRIPTF=/project/nmancuso_8/zeyunlu/projects

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

EUR_COVAR=${SCRATCH}/mesa_rnaseq_covar_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_rnaseq_covar_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_rnaseq_covar_HIS.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_cv_rnaseq/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 5 *int(int($NR-1)))"`
stop=$((start + 4))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_rnaseq_gene_list_noMHC.tsv`
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
  if [ -f ${wfile} ]; then
    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

    for pop in EUR AFR HIS
    do
      # split people
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

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict_old.py \
    --pheno $TMPDIR/${ID}.EUR.pheno $TMPDIR/${ID}.AFR.pheno $TMPDIR/${ID}.HIS.pheno \
    --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
    --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
    --trait ${ID} \
    --seed $IDX \
    --out $TMPDIR \
    --out_r2 $OUT/r2/sushie.${ID}.cvr2.tsv

    Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/old_pred_mesusie_3pop.R \
      $TMPDIR $ID $OUT/r2/mesusie.${ID}.cvr2.tsv

  rm -rf ${TMPDIR}
  fi
done

rm -rf ${bigTMP}
