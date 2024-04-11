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

plink=/project/nmancuso_8/zeyunlu/tools/plink/plink
plink2=/project/nmancuso_8/zeyunlu/tools/plink2
path=/project/nmancuso_8/data/GEUVADIS/GenotypeData/plink/ld_pruned/

LIST=${path}/grm_all_list.txt

if [ -f ${LIST} ]
then
  rm -rf ${LIST}
fi
touch ${LIST}

for NR in `seq 1 22`
do
  vcf=/project/nmancuso_8/data/GEUVADIS/GenotypeData/annotated_vcf_dbsnpv153/GEUVADIS.chr${NR}.annotated.dbSNP_v153.vcf.gz
  out=${path}/geuvadis_genotype_chr${NR}_ld_prune
  ${plink} --vcf ${vcf} -list-duplicate-vars suppress-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 --indep-pairwise 200 1 0.3 --make-bed --out ${out}
  echo "${out}" >> ${LIST}
done

pc_res=/project/nmancuso_8/data/GEUVADIS/ProcessedData/geno_pcs/
$plink2 --pmerge-list ${LIST} bfile --pca 30 --out ${pc_res}/geuvadis_genotype_all_pcs

