#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2Gb
#SBATCH --array=1-116
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1274
# 1276 = 116 * 11

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 11 *int(int($NR-1)))"`
stop=$((start + 10))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_proteins_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  PROTEINS_ID=${15}
  STRAND=$5
  TSS=$6
  TES=$7

  echo SUMMARY ${IDX}, ${PROTEINS_ID}, ${STRAND}, ${TSS}, ${TES}

  FILE=/scratch1/zeyunlu/sushie_proteins/sushie/alphas/${PROTEINS_ID}.normal.sushie.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_proteins/tss/tss.proteins.${PROTEINS_ID}.tsv

  if [ -e $FILE ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/calc_dist.py \
      --alpha-path $FILE \
      --gene ${PROTEINS_ID} \
      --strand ${STRAND} \
      --tss $TSS \
      --tes $TES \
      --tss-out $OUT
  fi
done
