#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

plink=/project/nmancuso_8/zeyunlu/tools/plink2
out=/project/nmancuso_8/data/INTERVAL/genotype/ld_pruned/INTERVAL_SOMALOGIC_POSTQC_ALL_final_v1.dbsnp153.ld_pruned
pc_res=/project/nmancuso_8/data/INTERVAL/covariates/geno_pcs
$plink --bfile ${out} --pca 30 --out ${pc_res}/interval_genotype_all_pcs
