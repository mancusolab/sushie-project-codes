#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --mem=8Gb
#SBATCH --array=1-98
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1274
# 1274 = 98 * 13

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/scratch1/zeyunlu/sushie

EUR_COVAR=${SCRATCH}/mesa_proteins_covar_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_proteins_covar_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_proteins_covar_HIS.tsv.gz

start=`python -c "print( 1 + 13 *int(int($NR-1)))"`
stop=$((start + 12))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_proteins_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  ID=${15}
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  TMPDIR=/scratch1/zeyunlu/temp_pred_proteins/${ID}

  # get genotype data
  bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis
  wfile=/scratch1/zeyunlu/sushie_proteins/weights/${ID}.normal.sushie.weights.tsv

  if [ -f ${wfile} ]; then
    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

    for pop in EUR AFR HIS
    do
      echo "Getting geno and pheno data for ${pop}"

      KEEP=${SCRATCH}/mesa_proteins_pt_${pop}.tsv
      PHENO=${SCRATCH}/mesa_proteins_v1_${pop}.tsv.gz

      ${PLINK} --silent --bfile ${bfile} --extract $TMPDIR/${ID}.snp --rm-dup force-first \
        --make-bed --keep $KEEP --out $TMPDIR/${ID}.${pop}.geno

      python ${SCRATCH}/extract_pheno.py \
      --file $PHENO \
      --gene $ID \
      --subject $KEEP \
      --out $TMPDIR/${ID}.${pop}.pheno
    done

    OUT=/project/nmancuso_8/data/sushie/aou_reg/mesa.proteins/chr${CHR}/mesa.proteins.chr${CHR}.${ID}.sushie.weights.tsv
    OUT2=/scratch1/zeyunlu/sushie_proteins/r2/$ID.cvr2.tsv

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict.py \
      --pheno $TMPDIR/${ID}.EUR.pheno.sushie.pheno $TMPDIR/${ID}.AFR.pheno.sushie.pheno $TMPDIR/${ID}.HIS.pheno.sushie.pheno \
      --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
      --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
      --trait ${ID} \
      --w_files /scratch1/zeyunlu/sushie_proteins/weights/${ID} \
      --seed $IDX \
      --out $OUT \
      --out_r2 $OUT2
  fi
done
