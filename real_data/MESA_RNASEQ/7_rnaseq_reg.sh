#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --mem=8Gb
#SBATCH --array=1041,1043,1075,1137,1176,1179,1181,1210,1221,1244,1246,1247,1319,1346,902,998,999
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 21855
# 21855 = 1457 * 15

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/scratch1/zeyunlu/sushie
EUR_COVAR=${SCRATCH}/mesa_rnaseq_covar_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_rnaseq_covar_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_rnaseq_covar_HIS.tsv.gz

start=`python -c "print( 1 + 15 *int(int($NR-1)))"`
stop=$((start + 14))

bigTMP=/scratch1/zeyunlu/temp_pred_rnaseq/tempf_${NR}

mkdir ${bigTMP}

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
  wfile=/scratch1/zeyunlu/sushie_rnaseq/weights/${ID}.normal.sushie.weights.tsv
  OUT=/project/nmancuso_8/data/sushie/aou_reg/mesa.mrna/chr${CHR}/mesa.mrna.chr${CHR}.${ID}.sushie.weights.tsv
  OUT2=/scratch1/zeyunlu/sushie_rnaseq/r2/$ID.cvr2.tsv

  if [ ! -f ${OUT2} ]; then
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

        python ${SCRATCH}/extract_pheno.py \
        --file $PHENO \
        --gene $NAME \
        --subject $KEEP \
        --out $TMPDIR/${ID}.${pop}.pheno
      done

      python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict.py \
        --pheno $TMPDIR/${ID}.EUR.pheno.sushie.pheno $TMPDIR/${ID}.AFR.pheno.sushie.pheno $TMPDIR/${ID}.HIS.pheno.sushie.pheno \
        --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
        --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
        --trait ${ID} \
        --w_files /scratch1/zeyunlu/sushie_rnaseq/weights/${ID} \
        --seed $IDX \
        --out $OUT \
        --out_r2 $OUT2
    fi
  fi
done
