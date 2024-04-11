#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=1Gb
#SBATCH --array=1-531
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 21240
# 21240 = 531 * 40

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/scratch1/zeyunlu/sushie

start=`python -c "print( 1 + 40 *int(int($NR-1)))"`
stop=$((start + 39))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d" ${SCRATCH}/v5_overlap_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  ID=$1

  echo SUMMARY $NR, $IDX, $ID

  # sushie
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/alphas/${ID}.normal.sushie.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/alphas/${ID}.normal.sushie.alphas.tsv

  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid2/valid.normal.${ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${ID} \
      --gene2 ${ID} \
      --seed $((IDX + NR*50000)) \
      --out $OUT
  else
    echo "both sushie files do not exist."
  fi

  # indep
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/alphas/${ID}.indep.sushie.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/alphas/${ID}.indep.sushie.alphas.tsv

  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid2/valid.indep.${ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${ID} \
      --gene2 ${ID} \
      --seed $((IDX + NR*50000)) \
      --out $OUT
  else
    echo "both indep files do not exist."
  fi

  # meta
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/alphas/${ID}.normal.meta.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/alphas/${ID}.normal.meta.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid2/valid.meta.${ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${ID} \
      --gene2 ${ID} \
      --seed $((IDX + NR*50000)) \
      --ans 3 \
      --out $OUT
  else
    echo "both meta files do not exist."
  fi

  # mega
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/alphas/${ID}.normal.mega.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/alphas/${ID}.normal.mega.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid2/valid.mega.${ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${ID} \
      --gene2 ${ID} \
      --seed $((IDX + NR*50000)) \
      --out $OUT
  else
    echo "both mega files do not exist."
  fi

done
