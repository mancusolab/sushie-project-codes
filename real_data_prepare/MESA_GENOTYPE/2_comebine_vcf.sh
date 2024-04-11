#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1GB
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

geno=/project/nmancuso_8/data/MESA/processed/geno/vcf/filtered

geno_file1=${geno}/TOPMED.WGS.filtered.chr${NR}.hg38.c1.vcf.gz
geno_file2=${geno}/TOPMED.WGS.filtered.chr${NR}.hg38.c2.vcf.gz

out=/project/nmancuso_8/data/MESA/processed/geno/vcf/combined/TOPMED.WGS.chr${NR}.combined.vcf.gz
bcftools merge ${geno_file1} ${geno_file2} -Oz -o ${out}

