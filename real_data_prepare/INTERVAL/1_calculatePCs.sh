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

# split
bfile=/project/nmancuso_8/data/INTERVAL/genotype/INTERVAL_SOMALOGIC_POSTQC_ALL_final_v1.dbsnp153
for chr in `seq 1 22`
do
  out=/project/nmancuso_8/data/INTERVAL/genotype/all_final_split/INTERVAL_SOMALOGIC_POSTQC_chr${chr}_final_v1.dbsnp153
  ${plink} --bfile ${bfile} --chr $chr --make-bed --out ${out}
done

# fusion split
bfile=/project/nmancuso_8/data/INTERVAL/genotype/INTERVAL_SOMALOGIC_POSTQC_ALL_final_v1.dbsnp153.FUSION
for chr in `seq 1 22`
do
  out=/project/nmancuso_8/data/INTERVAL/genotype/all_final_fusion_split/INTERVAL_SOMALOGIC_POSTQC_chr${chr}_final_v1.dbsnp153.FUSION
  ${plink} --bfile ${bfile} --chr $chr --make-bed --out ${out}
done

## calculate genotype PCs
# no need to provide genotype PCs as they already provided 3
#
#path=/project/nmancuso_8/data/INTERVAL/genotype/ld_pruned/
#bfile=/project/nmancuso_8/data/INTERVAL/genotype/INTERVAL_SOMALOGIC_POSTQC_ALL_final_v1.dbsnp153
#out=/project/nmancuso_8/data/INTERVAL/genotype/ld_pruned/INTERVAL_SOMALOGIC_POSTQC_ALL_final_v1.dbsnp153.ld_pruned
#${plink} --silent --bfile ${bfile} --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 --indep-pairwise 200 1 0.3 --make-bed --out ${out}
#
#pc_res=/project/nmancuso_8/data/INTERVAL/covariates/geno_pcs
#$plink --bfile ${out} --pca 30 --out ${pc_res}/interval_genotype_all_pcs
#
