#!/bin/bash

out=/project/nmancuso_8/data/GEUVADIS/GenotypeData/plink/raw
PLINK=/project/nmancuso_8/zeyunlu/tools/plink/plink
for chr in `seq 1 22`
do
  vcf=/project/nmancuso_8/data/GEUVADIS/GenotypeData/annotated_vcf_dbsnpv153/GEUVADIS.chr${chr}.annotated.dbSNP_v153.vcf.gz
  ${PLINK} --silent --vcf ${vcf} --make-bed --out ${out}/GEUVADIS.chr${chr}.annotated.dbSNP_v153
done
