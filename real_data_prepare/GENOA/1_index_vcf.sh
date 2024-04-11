#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=20GB
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

# index file
cd /project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/ea/ea_Affymetrix
bcftools index chr${NR}.vcf.gz

cd /project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/ea/ea_IlluminaALL
bcftools index chr${NR}.vcf.gz

cd /project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/ea/ea_IlluminaALL_geno0.15
bcftools index chr${NR}.vcf.gz

cd /project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/aa/aa_Affymetrix
bcftools index chr${NR}.vcf.gz

cd /project/nmancuso_8/data/GENOA/processed/genotype/vcf/imputed/aa/aa_Illumina
bcftools index chr${NR}.vcf.gz
