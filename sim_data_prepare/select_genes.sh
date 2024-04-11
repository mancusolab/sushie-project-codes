#!/bin/bash

source /home1/zeyunlu/init.sh
conda activate jax2

# this script is to select genes that fits for 3 ancestires and other criterion, and have more than 100 SNPs available.
hapmap=/project/nmancuso_8/data/HapMap_SNPs

output=/scratch1/zeyunlu/sushie_sim_genes/
plink=/project/nmancuso_8/zeyunlu/tools/plink2
TG=/project/nmancuso_8/data/1000GP3_multiPop_allelesAligned

tmpoutput=/scratch1/zeyunlu/genes55

NR=0
IDX=1
while [ ! -e ${output}/EAS/EAS_gene500_*_snps.bim ] && [ ! -e ${output}/EUR/EUR_gene500_*_snps.bim ] && [ ! -e ${output}/AFR/AFR_gene500_*_snps.bim ]
do
  NR=`python -c "print(int($NR+1))"`
  params=`sed "${NR}q;d" /home1/zeyunlu/github/sushie-data-codes/sim_data_prepare/shuffle_genes_list.tsv`
  set -- junk $params
  shift
  name=$1
  CHR=$2
  start=`python -c "print(int(max($3 - 500000, 1)))"`
  stop=`python -c "print( int(int($4) + 500000))"`

  for pop in EAS AFR EUR
  do
    $plink --bfile ${TG}/${pop}/1000G.${pop}.QC.allelesAligned.${CHR} \
    --chr $CHR \
    --from-bp $start \
    --to-bp $stop \
    --make-bed \
    --out ${tmpoutput}/${pop}_gene${IDX}_${name}_snps \
    --snps-only \
    --silent \
    --geno 0.01 \
    --hwe midp 1e-6 \
    --maf 0.01 \
    --allow-no-sex \
    --extract $hapmap/hm.${CHR}.snp \
    --force-intersect
  done

  if [ -e ${tmpoutput}/EAS_gene${IDX}_${name}_snps.bim ] && [ -e ${tmpoutput}/EUR_gene${IDX}_${name}_snps.bim ] && [ -e ${tmpoutput}/AFR_gene${IDX}_${name}_snps.bim ]
  then

    n_common=$(python /home1/zeyunlu/github/sushie-data-codes/sim_data_prepare/common_snps.py \
      --file1 ${tmpoutput}/EAS_gene${IDX}_${name}_snps.bim \
      --file2 ${tmpoutput}/EUR_gene${IDX}_${name}_snps.bim \
      --file3 ${tmpoutput}/AFR_gene${IDX}_${name}_snps.bim)

    threshold=500
    if [ $n_common -gt $threshold ]
    then
      mv ${tmpoutput}/EUR_gene${IDX}_${name}_snps* ${output}/EUR/
      mv ${tmpoutput}/EAS_gene${IDX}_${name}_snps* ${output}/EAS/
      mv ${tmpoutput}/AFR_gene${IDX}_${name}_snps* ${output}/AFR/
      echo "NR=${NR} IDX=${IDX} common SNPs ${n_common} ${name}"
      IDX=`python -c "print(int($IDX+1))"`
    fi
  fi
done




