#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --mem=4Gb
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
# 13905 = 927 * 15

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 15 *int(int($NR-1)))"`
stop=$((start + 14))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d" ${SCRATCH}/genoa_sushie_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  ID=$2
  NAME=$3
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${NAME}, ${CHR}, ${P0}, ${P1}

  wfile=/scratch1/zeyunlu/sushie_genoa/sushie/weights/${ID}.normal.sushie.weights.tsv

  if [ -f ${wfile} ]; then

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict_zenodo.py \
      --trait ${ID} \
      --w_files /scratch1/zeyunlu/sushie_genoa/ \
      --study genoa.mrna \
      --out /project/nmancuso_8/data/sushie/zenodo_weights/genoa.mrna/chr${CHR}/genoa.mrna.${ID}.sushie.weights.tsv

  fi
done
