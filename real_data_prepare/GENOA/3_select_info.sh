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

cd /project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/merged

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/R2\n" ea_chr${NR}.vcf.gz -o /project/nmancuso_8/zeyunlu/result/ea_${NR}_r2.tsv

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/R2\n" aa_chr${NR}.vcf.gz -o /project/nmancuso_8/zeyunlu/result/aa_${NR}_r2.tsv
