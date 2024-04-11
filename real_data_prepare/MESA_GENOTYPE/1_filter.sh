#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --array=1-44
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

params=`sed "${NR}q;d" /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/MESA_GENOTYPE/1_file.tsv`
echo "NR=${NR}: ${params}"
set -- junk $params
shift

# PARAMETERS
file=$1
out=$2

path=/project/nmancuso_8/data/MESA/processed/geno/vcf/filtered

bcftools view -i "FILTER='PASS'" -m2 -M2 -v snps,indels ${file} -Ob \
 | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o ${path}/${out}

