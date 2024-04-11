#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2Gb
#SBATCH --array=1-22
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  chr=$1
  SCRATCHDIR=`echo $PWD`
else
  chr=$SLURM_ARRAY_TASK_ID
fi

source /home1/zeyunlu/init.sh
conda activate jax2

mesa_geno=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${chr}.filtered.analysis.bim
genoa_geno=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/ea_chr${chr}.bim

for geno in ${mesa_geno} ${genoa_geno}
do
  if [ "$geno" == "$mesa_geno" ]
  then
    geno_name="mesa"
  else
    geno_name="genoa"
  fi

  cd /scratch1/zeyunlu/new_annot/original
  for file in *
  do
    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    ${file} \
    ${geno} \
    ${chr} \
    -f -l \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.ldsc.chr${chr}.${file}
  done

  # v2
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ENCODE/raw/GRCh38-cCREs.v2.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "," -l -d "CTCF-only" \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.encode.chr${chr}.ccre.v2.bed.gz

  # v3
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ENCODE/raw/GRCh38-cCREs.v3.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "," -l -d "CTCF-only" \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.encode.chr${chr}.ccre.v3.bed.gz

  # v4
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ENCODE/raw/GRCh38-cCREs.v4.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "-" -l \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.encode.chr${chr}.ccre.v4.bed.gz

  # chiou et al
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ANNOTATIONS/CHIOU_ET_AL/chiou_et_al_grch38.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "," -l \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.chiou.chr${chr}.ccre.bed.gz

  # satpathy et al
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ANNOTATIONS/SATPATHY_ET_AL/snATAC-seq-fresh-peaks_grch38.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "," -l \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.satpathy.chr${chr}.fresh.peaks.bed.gz

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ANNOTATIONS/SATPATHY_ET_AL/snATAC-seq-frozen-peaks_grch38.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "," -l \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.satpathy.chr${chr}.frozen.peaks.bed.gz

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data_prepare/annotation/make.annot.encode.py \
    /project/nmancuso_8/data/ANNOTATIONS/SATPATHY_ET_AL/snATAC-seq-frozen-sort-peaks_grch38.bed.gz \
    ${geno} \
    ${chr} \
    -f -s "," -l \
    -o /scratch1/zeyunlu/new_annot/processed/${geno_name}/${geno_name}.satpathy.chr${chr}.frozen.sort.peaks.bed.gz
done









