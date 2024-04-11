#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
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

out=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/merged

ea1=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/ea/ea_Affymetrix
ea2=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/ea/ea_IlluminaALL

bcftools merge $ea1/chr${NR}.vcf.gz $ea2/chr${NR}.vcf.gz --info-rules R2:join -o $out/ea_chr${NR}.vcf.gz

aa1=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/aa/aa_Affymetrix
aa2=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/aa/aa_Illumina

bcftools merge $aa1/chr${NR}.vcf.gz $aa2/chr${NR}.vcf.gz --info-rules R2:join -o $out/aa_chr${NR}.vcf.gz

