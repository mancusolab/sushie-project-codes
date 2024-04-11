#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --array=1-4
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

plink=/project/nmancuso_8/zeyunlu/tools/plink2
path=/project/nmancuso_8/data/MESA/processed/geno/plink/ld_pruned

params=`sed "${NR}q;d" ${path}/ancestry_list.tsv`
set -- junk $params
shift
pop=$1

echo "Running for ${pop}"

LIST=${path}/by_ancestry/${pop}/grm_${pop}_list.txt

if [ -f ${LIST} ]
then
  rm -rf ${LIST}
fi

touch ${LIST}

KEEP=${path}/${pop}_genoid.tsv

for chr in `seq 1 22`
do
  vcf=/project/nmancuso_8/data/MESA/processed/geno/vcf/combined/TOPMED.WGS.chr${chr}.combined.vcf.gz
  out=${path}/by_ancestry/${pop}/mesa_genotype_${pop}_chr${chr}_ld_prune
  ${plink} --silent --vcf ${vcf} --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 --keep $KEEP --indep-pairwise 200 1 0.3 --make-bed --out ${out}

  echo "${out}" >> ${LIST}
done

pc_res=/project/nmancuso_8/data/MESA/processed/covariates/genotype_pcs/by_ancestry
$plink --pmerge-list ${LIST} bfile --pca 30 --out ${pc_res}/mesa_genotype_${pop}_pcs

