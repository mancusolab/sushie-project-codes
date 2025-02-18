#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=3Gb
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

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

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
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/sushie/alphas/${ID}.normal.sushie.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/sushie/alphas/${ID}.normal.sushie.alphas.tsv

  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid/valid.normal.${ID}.tsv

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
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/sushie/alphas/${ID}.indep.sushie.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/sushie/alphas/${ID}.indep.sushie.alphas.tsv

  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid/valid.indep.${ID}.tsv

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
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/sushie/alphas/${ID}.normal.meta.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/sushie/alphas/${ID}.normal.meta.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid/valid.meta.${ID}.tsv

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
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/sushie/alphas/${ID}.normal.mega.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/sushie/alphas/${ID}.normal.mega.alphas.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid/valid.mega.${ID}.tsv

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

  # mesusie
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/mesusie/weights/mesusie.${ID}.weights.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/mesusie/weights/mesusie.${ID}.weights.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid/valid.mesusie.${ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp2.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${ID} \
      --gene2 ${ID} \
      --seed $((IDX + NR*50000)) \
      --out $OUT
  else
    echo "both mesusie files do not exist."
  fi

  # susiex
  FILE1=/scratch1/zeyunlu/sushie_rnaseq/susiex/weights/susiex.${ID}.weights.tsv
  FILE2=/scratch1/zeyunlu/sushie_v5/susiex/weights/susiex.${ID}.weights.tsv
  OUT=/scratch1/zeyunlu/sushie_rnaseq/valid/valid.susiex.${ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp2.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${ID} \
      --gene2 ${ID} \
      --seed $((IDX + NR*50000)) \
      --out $OUT
  else
    echo "both susiex files do not exist."
  fi

done
