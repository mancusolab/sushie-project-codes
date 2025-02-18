#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu


if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/MESA_RNASEQ/prepare_ldsc_annot_genoa.R
