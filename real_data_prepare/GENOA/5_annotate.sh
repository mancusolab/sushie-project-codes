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

cd /project/nmancuso_8/zeyunlu/result

file=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/merged
out=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/weighted_r2

# EA
bgzip ea_${NR}_newr2.tsv
tabix -s1 -b2 -e2 ea_${NR}_newr2.tsv.gz
bcftools annotate -a ea_${NR}_newr2.tsv.gz -h header -c CHROM,POS,REF,ALT,WEIGHTEDR2 ${file}/ea_chr${NR}.vcf.gz -o ${out}/ea_chr${NR}_anno_weighted_r2.vcf.gz

# AA
bgzip aa_${NR}_newr2.tsv
tabix -s1 -b2 -e2 aa_${NR}_newr2.tsv.gz
bcftools annotate -a aa_${NR}_newr2.tsv.gz -h header -c CHROM,POS,REF,ALT,WEIGHTEDR2 ${file}/aa_chr${NR}.vcf.gz -o ${out}/aa_chr${NR}_anno_weighted_r2.vcf.gz

