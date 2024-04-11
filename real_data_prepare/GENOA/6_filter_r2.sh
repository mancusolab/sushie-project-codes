#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --array=1-22
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu


if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source /home1/zeyunlu/init.sh
conda activate bcf


file=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/weighted_r2
out=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/filter_06

bcftools filter ${file}/ea_chr${NR}_anno_weighted_r2.vcf.gz -e "INFO/WEIGHTEDR2<0.6" -o ${out}/ea_chr${NR}_weighted_r2_06.vcf.gz

bcftools filter ${file}/aa_chr${NR}_anno_weighted_r2.vcf.gz -e "INFO/WEIGHTEDR2<0.6" -o ${out}/aa_chr${NR}_weighted_r2_06.vcf.gz
