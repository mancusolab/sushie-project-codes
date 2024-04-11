#!/bin/bash

out=/project/nmancuso_8/data/sushie/plink
KEEP=/scratch1/zeyunlu/sushie/mesa_all_pt.tsv
PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
for chr in `seq 1 22`
do
  vcf=/project/nmancuso_8/data/MESA/processed/geno/vcf/combined/TOPMED.WGS.chr${chr}.combined.vcf.gz
  ${PLINK} --silent --vcf ${vcf} --keep $KEEP --make-bed --out ${out}/TOPMED.WGS.chr${chr}.filtered.analysis
done
