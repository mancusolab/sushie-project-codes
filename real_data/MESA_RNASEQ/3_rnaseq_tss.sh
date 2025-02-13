#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=2Gb
#SBATCH --array=1-625
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
# 21875 = 625 * 35

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 35 *int(int($NR-1)))"`
stop=$((start + 34))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_rnaseq_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  STRAND=$5
  TSS=$6
  TES=$7
  RNASEQ_ID=${12}

  echo SUMMARY ${IDX}, ${RNASEQ_ID}, ${STRAND}, ${TSS}, ${TES}

  FILE=/scratch1/zeyunlu/sushie_rnaseq/sushie/alphas/${RNASEQ_ID}.normal.sushie.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/tss/tss.rnaseq.${RNASEQ_ID}.tsv

  if [ -e $FILE ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/calc_dist.py \
      --alpha-path $FILE \
      --gene ${RNASEQ_ID} \
      --strand ${STRAND} \
      --tss $TSS \
      --tes $TES \
      --tss-out $OUT
  fi
done
