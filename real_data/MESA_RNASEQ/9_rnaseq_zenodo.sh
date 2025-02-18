#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --mem=4Gb
#SBATCH --array=1-4371
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

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

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

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

  wfile=/scratch1/zeyunlu/sushie_rnaseq/sushie/weights/${ID}.normal.sushie.weights.tsv

  if [ -f ${wfile} ]; then

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict_zenodo.py \
      --trait ${ID} \
      --w_files /scratch1/zeyunlu/sushie_rnaseq/ \
      --study mesa.mrna \
      --out /project/nmancuso_8/data/sushie/zenodo_weights/mesa.mrna/chr${CHR}/mesa.mrna.${ID}.sushie.weights.tsv

  fi
done
