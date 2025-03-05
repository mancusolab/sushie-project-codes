#!/bin/bash

# output 1000G SNPs to a separate folder
cd /project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/EUR

for file in *.bim; do
    CHR=$(awk 'NR==1 {print $1}' "$file")
    awk '{print $2}' "$file" > "/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/mesa_snps/mesa_${CHR}_$file"
done

# select these SNPs in the TOPMed-MESA genotypes
cd /project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/mesa_snps
SCRATCH=/project/nmancuso_8/data/sushie/meta_data
PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
OUT=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes

for file in *.bim
do
  for pop in EUR AFR
  do
    echo "Getting geno data for ${file}...${pop}"
    CHR=$(echo "${file}" | grep -oP "(?<=mesa_)\d+(?=_EUR)")
    bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis
    KEEP=${SCRATCH}/mesa_rnaseq_pt_${pop}.tsv

    ${PLINK} --silent --bfile ${bfile} --extract ${file} \
      --make-bed --keep $KEEP --out $OUT/${pop}_mesa/${file}.${pop}.geno
  done
done

# make sure the simulation SNPs are the same for the two studies
cd /project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/EUR
for file in *.bim; do
    CHR=$(awk 'NR==1 {print $1}' "$file")

done

FOLDER=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes
for IDX in `seq 1 500`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_data_prepare/sim_gene_index_list.tsv`
  set -- junk $params
  shift
  INDEX=$1
  CHR=$2
  NAME=$3
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_data_prepare/common_snps_for_two_studies.py \
      --file1 $FOLDER/EUR/EUR_gene${INDEX}_${NAME}_snps.bim \
      --file2 $FOLDER/AFR/AFR_gene${INDEX}_${NAME}_snps.bim \
      --file3 $FOLDER/EUR_mesa/mesa_${CHR}_EUR_gene${INDEX}_${NAME}_snps.bim.EUR.geno.bim \
      --file4 $FOLDER/AFR_mesa/mesa_${CHR}_EUR_gene${INDEX}_${NAME}_snps.bim.AFR.geno.bim \
      --file5 $FOLDER/common_snps/common_snps_index${INDEX}_chr${CHR}_${NAME}.geno.bim
done

FOLDER=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes
for IDX in `seq 1 500`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_data_prepare/sim_gene_index_list.tsv`
  set -- junk $params
  shift
  INDEX=$1
  CHR=$2
  NAME=$3

  $PLINK --bfile $FOLDER/EUR/EUR_gene${INDEX}_${NAME}_snps \
    --extract $FOLDER/common_snps/common_snps_index${INDEX}_chr${CHR}_${NAME}.geno.bim \
    --make-bed --out $FOLDER/new/EUR_1000G/${NAME}_${CHR}_geno

  $PLINK --bfile $FOLDER/AFR/AFR_gene${INDEX}_${NAME}_snps \
    --extract $FOLDER/common_snps/common_snps_index${INDEX}_chr${CHR}_${NAME}.geno.bim \
    --make-bed --out $FOLDER/new/AFR_1000G/${NAME}_${CHR}_geno

  $PLINK --bfile $FOLDER/EUR_mesa/mesa_${CHR}_EUR_gene${INDEX}_${NAME}_snps.bim.EUR.geno \
    --extract $FOLDER/common_snps/common_snps_index${INDEX}_chr${CHR}_${NAME}.geno.bim \
    --make-bed --out $FOLDER/new/EUR_mesa/${NAME}_${CHR}_geno

  $PLINK --bfile $FOLDER/AFR_mesa/mesa_${CHR}_EUR_gene${INDEX}_${NAME}_snps.bim.AFR.geno \
    --extract $FOLDER/common_snps/common_snps_index${INDEX}_chr${CHR}_${NAME}.geno.bim \
    --make-bed --out $FOLDER/new/AFR_mesa/${NAME}_${CHR}_geno
done
