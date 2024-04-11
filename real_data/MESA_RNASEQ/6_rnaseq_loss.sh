#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --mem=4Gb
#SBATCH --array=1-547
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
# 21855 < 21880 = 537 * 40

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/scratch1/zeyunlu/sushie

start=`python -c "print( 1 + 40 *int(int($NR-1)))"`
stop=$((start + 39))

bigTMP=/scratch1/zeyunlu/loss_rnaseq/tempf_${NR}

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

  wfile=/scratch1/zeyunlu/sushie_rnaseq/weights/${ID}.normal.sushie.weights.tsv

  if [ -f ${wfile} ]; then
    TMPDIR=${bigTMP}/${ID}

    mkdir -p $TMPDIR

    # get genotype data
    bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis

    for pop in EUR AFR HIS EAS
    do
      echo "Getting geno and pheno data for ${pop}"

      KEEP=${SCRATCH}/mesa_rnaseq_pt_${pop}.tsv

      ${PLINK} --silent --bfile ${bfile} --chr $CHR --from-bp $P0 --to-bp $P1 \
        --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 \
        --make-bed --keep $KEEP --out $TMPDIR/${ID}.${pop}.geno
    done

    OUT=/scratch1/zeyunlu/sushie_rnaseq/loss/$ID.loss.tsv

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/loss.py \
        --eurfile $TMPDIR/${ID}.EUR.geno.bim \
        --afrfile $TMPDIR/${ID}.AFR.geno.bim \
        --hisfile $TMPDIR/${ID}.HIS.geno.bim \
        --easfile $TMPDIR/${ID}.EAS.geno.bim \
        --gene $ID \
        --out $OUT
  fi
done
