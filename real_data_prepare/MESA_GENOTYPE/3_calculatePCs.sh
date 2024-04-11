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
path=/project/nmancuso_8/data/MESA/processed/geno/plink/ld_pruned/all

LIST=${path}/grm_all_list.txt

if [ -f ${LIST} ]
then
  rm -rf ${LIST}
fi
touch ${LIST}

for NR in `seq 1 22`
do
  vcf=/project/nmancuso_8/data/MESA/processed/geno/vcf/combined/TOPMED.WGS.chr${NR}.combined.vcf.gz
  out=${path}/mesa_genotype_chr${NR}_ld_prune
  ${plink} --silent --vcf ${vcf} --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 --indep-pairwise 200 1 0.3 --make-bed --out ${out}

  echo "${out}" >> ${LIST}
done

pc_res=/project/nmancuso_8/data/MESA/processed/covariates/genotype_pcs/all
$plink --pmerge-list ${LIST} bfile --pca 30 --out ${pc_res}/mesa_genotype_all_pcs

