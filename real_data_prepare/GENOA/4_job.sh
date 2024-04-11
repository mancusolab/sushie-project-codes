#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=40GB
#SBATCH --array=2,4
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu


if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source /home1/zeyunlu/init.sh

Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/GENOA/4_split_r2.R $NR
