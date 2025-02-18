#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=1Gb
#SBATCH --array=1-146
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1313
# 1314 = 9 * 146

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 9 *int(int($NR-1)))"`
stop=$((start + 8))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d" ${SCRATCH}/interval_overlap_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  pro_ID=$3
  int_ID=$2

  echo SUMMARY $NR, $IDX, $pro_ID, $int_ID

  # sushie
  FILE1=/scratch1/zeyunlu/sushie_proteins/sushie/alphas/${pro_ID}.normal.sushie.alphas.tsv
  FILE2=/scratch1/zeyunlu/sushie_interval/sushie/alphas/${int_ID}.normal.sushie.alphas.tsv

  OUT=/scratch1/zeyunlu/sushie_proteins/valid/valid.normal.${pro_ID}.tsv

  if [ -e $FILE1 ] && [ -e $FILE2 ];
  then
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/valid_snp.py \
      --file_w1 $FILE1 \
      --file_w2 $FILE2 \
      --gene1 ${pro_ID} \
      --gene2 ${int_ID} \
      --seed $((IDX + NR*50000)) \
      --out $OUT
  else
    echo "both sushie files do not exist."
  fi

done
