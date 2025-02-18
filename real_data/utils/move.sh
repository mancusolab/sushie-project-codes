#!/bin/bash

# proteins
for IDX in `seq 1 21088`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  /project/nmancuso_8/data/sushie/aou_csonly/metadata.tsv`
  set -- junk $params
  shift
  CHR=$1
  NAME=$2
  STUDY=$3

  echo SUMMARY $IDX ${CHR}, ${NAME}, ${STUDY}

  wfolder=/project/nmancuso_8/data/sushie/aou_csonly/weights/chr${CHR}/
  sfolder=/project/nmancuso_8/data/sushie/aou_csonly/snps/chr${CHR}/

  cp /project/nmancuso_8/data/sushie/all_of_us/snps/${STUDY}/chr${CHR}/${STUDY}.chr${CHR}.${NAME}.sushie.snps.tsv ${sfolder}
  cp /project/nmancuso_8/data/sushie/all_of_us/weights/${STUDY}/chr${CHR}/${STUDY}.chr${CHR}.${NAME}.sushie.weights.tsv ${wfolder}
done



for study in mesa.mrna mesa.proteins genoa.mrna
do
  for chr in `seq 1 22`
  do
    echo $study, $chr
    mv snps/chr${chr}/${study}* snps_bp/chr${chr}/
    mv weights/chr${chr}/${study}* weights_bp/chr${chr}/
  done
done


for file in /scratch1/zeyunlu/sushie_rnaseq/weights/*.normal.meta.weights.tsv; do
  cp $file .
done


for file in /scratch1/zeyunlu/sushie_rnaseq/cs/*.indep.sushie.cs.tsv; do
  cp $file .
done

for file in /scratch1/zeyunlu/sushie_rnaseq/alphas/*.normal.meta.alphas.tsv; do
  cp $file .
done


