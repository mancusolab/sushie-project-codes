#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --mem=4Gb
#SBATCH --array=1-98
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1274
# 1274 = 98 * 13

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 13 *int(int($NR-1)))"`
stop=$((start + 12))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_proteins_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  ID=${15}
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  wfile=/scratch1/zeyunlu/sushie_proteins/sushie/weights/${ID}.normal.sushie.weights.tsv

  if [ -f ${wfile} ]; then

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict_zenodo.py \
      --trait ${ID} \
      --w_files /scratch1/zeyunlu/sushie_proteins/ \
      --study mesa.proteins \
      --out /project/nmancuso_8/data/sushie/zenodo_weights/mesa.proteins/chr${CHR}/mesa.proteins.${ID}.sushie.weights.tsv
  fi
done
