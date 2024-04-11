#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=2Gb
#SBATCH --array=1-927
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 13905
# 13905 = 15 * 927

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/scratch1/zeyunlu/sushie

start=`python -c "print( 1 + 15 *int(int($NR-1)))"`
stop=$((start + 14))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d" ${SCRATCH}/genoa_sushie_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  GENOA_ID=$2
  STRAND=$5
  TSS=$6
  TES=$7

  echo SUMMARY ${NR}, ${IDX}, ${GENOA_ID}, ${STRAND}, ${TSS}, ${TES}

  FILE=/scratch1/zeyunlu/sushie_genoa/alphas/${GENOA_ID}.normal.sushie.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_genoa/tss/tss.genoa.${GENOA_ID}.tsv

  if [ -e $FILE ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/calc_dist.py \
      --alpha-path $FILE \
      --gene ${GENOA_ID} \
      --strand $STRAND \
      --tss $TSS \
      --tes $TES \
      --tss-out $OUT
  fi
done
