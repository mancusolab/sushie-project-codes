#!/bin/bash

liftover=/project/nmancuso_8/zeyunlu/tools/liftOver
out=/project/nmancuso_8/data/LDSC/baseline_annotation/bed_hg38
chain=/project/nmancuso_8/data/liftover_chains/hg19ToHg38.over.chain.gz

cd /project/nmancuso_8/data/LDSC/baseline_annotation/bed
for FILE in *
do
  echo $FILE
  $liftover ${FILE} $chain ${out}/hg38_${FILE} /scratch1/zeyunlu/trash2/unmapped_${FILE}
done


# chiou et al
liftover=/project/nmancuso_8/zeyunlu/tools/liftOver
out=~/trash/chieou_grch38
chain=/project/nmancuso_8/data/liftover_chains/hg19ToHg38.over.chain.gz
$liftover ~/trash/chiou_grch37.bed $chain ${out} /scratch1/zeyunlu/trash2/unmapped_chiou


# satpathy et al
liftover=/project/nmancuso_8/zeyunlu/tools/liftOver
chain=/project/nmancuso_8/data/liftover_chains/hg19ToHg38.over.chain.gz

cd ~/trash/
for file in fresh_peaks.bed frozen_peaks.bed frozensort_peaks.bed
do
  $liftover $file $chain ~/trash/grch38_${file} /scratch1/zeyunlu/trash2/unmapped_$file
done


# linker
liftover=/project/nmancuso_8/zeyunlu/tools/liftOver
out1=~/trash/linker1_grch38
out2=~/trash/linker2_grch38
chain=/project/nmancuso_8/data/liftover_chains/hg19ToHg38.over.chain.gz
$liftover ~/trash/linker1.tsv $chain ${out1} /scratch1/zeyunlu/trash/unmapped_linker1
$liftover ~/trash/linker2.tsv $chain ${out2} /scratch1/zeyunlu/trash/unmapped_linker2
