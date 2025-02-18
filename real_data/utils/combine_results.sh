#!/bin/bash

out=~/data/sushie/real2

# rnaseq
for type in her corr
do
  echo rnaseq_${type}
  rm -rf ${out}/rnaseq_${type}.tsv.gz
  cd /scratch1/zeyunlu/sushie_rnaseq/sushie/${type}
  head -1 ENSG00000089818_NECAP1.normal.sushie.${type}.tsv > ${out}/rnaseq_${type}.tsv
  find . -name "*.sushie.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_${type}.tsv
  gzip  ${out}/rnaseq_${type}.tsv
done

for class in cs
do
  cd /scratch1/zeyunlu/sushie_rnaseq/sushie/${class}
  for type in normal.meta normal.mega normal.sushie indep.sushie
  do
    echo rnaseq_${type}_${class}
    rm -rf ${out}/rnaseq_${type}_${class}.tsv.gz
    head -1 ENSG00000089818_NECAP1.${type}.${class}.tsv > ${out}/rnaseq_${type}_${class}.tsv
    find . -name "*.${type}.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_${type}_${class}.tsv
    gzip ${out}/rnaseq_${type}_${class}.tsv
  done
done

for class in cs
do
  for type in susiex mesusie
  do
    cd /scratch1/zeyunlu/sushie_rnaseq/${type}/${class}
    echo rnaseq_${type}_${class}
    rm -rf ${out}/rnaseq_${type}_${class}.tsv.gz
    head -1 ${type}.ENSG00000001036_FUCA2.${class}.tsv > ${out}/rnaseq_${type}_${class}.tsv
    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_${type}_${class}.tsv
    gzip ${out}/rnaseq_${type}_${class}.tsv
  done
done

#for class in corr
#do
#  for type in xmap
#  do
#    cd /scratch1/zeyunlu/sushie_rnaseq/${type}/${class}
#    echo rnaseq_${type}_${class}
#    rm -rf ${out}/rnaseq_${type}_${class}.tsv.gz
#    head -1 ${type}.ENSG00000001036_FUCA2.${class}.tsv > ${out}/rnaseq_${type}_${class}.tsv
#    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_${type}_${class}.tsv
#    gzip ${out}/rnaseq_${type}_${class}.tsv
#  done
#done

# proteins
for type in her corr
do
  echo proteins_${type}
  rm -rf ${out}/proteins_${type}.tsv.gz
  cd /scratch1/zeyunlu/sushie_proteins/sushie/${type}
  head -1 ENSG00000285441_SL001815.normal.sushie.${type}.tsv > ${out}/proteins_${type}.tsv
  find . -name "*.sushie.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_${type}.tsv
  gzip ${out}/proteins_${type}.tsv
done

for class in cs
do
  cd /scratch1/zeyunlu/sushie_proteins/sushie/${class}
  for type in normal.meta normal.mega normal.sushie indep.sushie
  do
    echo proteins_${type}_${class}
    rm -rf ${out}/proteins_${type}_${class}.tsv.gz
    head -1 ENSG00000285441_SL001815.${type}.${class}.tsv > ${out}/proteins_${type}_${class}.tsv
    find . -name "*.${type}.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_${type}_${class}.tsv
    gzip ${out}/proteins_${type}_${class}.tsv
  done
done

for class in cs
do
  for type in susiex mesusie
  do
    cd /scratch1/zeyunlu/sushie_proteins/${type}/${class}
    echo proteins_${type}_${class}
    rm -rf ${out}/proteins_${type}_${class}.tsv.gz
    head -1 ${type}.ENSG00000115598_SL004875.${class}.tsv > ${out}/proteins_${type}_${class}.tsv
    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_${type}_${class}.tsv
    gzip ${out}/proteins_${type}_${class}.tsv
  done
done

#for class in corr
#do
#  for type in xmap
#  do
#    cd /scratch1/zeyunlu/sushie_proteins/${type}/${class}
#    echo proteins_${type}_${class}
#    rm -rf ${out}/proteins_${type}_${class}.tsv.gz
#    head -1 ${type}.ENSG00000115598_SL004875.${class}.tsv > ${out}/proteins_${type}_${class}.tsv
#    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_${type}_${class}.tsv
#    gzip ${out}/proteins_${type}_${class}.tsv
#  done
#done


# genoa
for type in her corr
do
  echo genoa_${type}
  cd /scratch1/zeyunlu/sushie_genoa/sushie/${type}
  rm -rf ${out}/genoa_${type}.tsv.gz
  head -1 ENSG00000288547.normal.sushie.${type}.tsv > ${out}/genoa_${type}.tsv
  find . -name "*.sushie.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_${type}.tsv
  gzip ${out}/genoa_${type}.tsv
done

for class in cs
do
  cd /scratch1/zeyunlu/sushie_genoa/sushie/${class}
  for type in normal.meta normal.mega normal.sushie indep.sushie
  do
    echo genoa_${type}_${class}
    rm -rf ${out}/genoa_${type}_${class}.tsv.gz
    head -1 ENSG00000288547.${type}.${class}.tsv > ${out}/genoa_${type}_${class}.tsv
    find . -name "*.${type}.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_${type}_${class}.tsv
    gzip ${out}/genoa_${type}_${class}.tsv
  done
done

for class in cs
do
  for type in susiex mesusie
  do
    cd /scratch1/zeyunlu/sushie_genoa/${type}/${class}
    echo genoa_${type}_${class}
    rm -rf ${out}/genoa_${type}_${class}.tsv.gz
    head -1 ${type}.ENSG00000266714.${class}.tsv > ${out}/genoa_${type}_${class}.tsv
    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_${type}_${class}.tsv
    gzip ${out}/genoa_${type}_${class}.tsv
  done
done

#for class in corr
#do
#  for type in xmap
#  do
#    cd /scratch1/zeyunlu/sushie_genoa/${type}/${class}
#    echo genoa_${type}_${class}
#    rm -rf ${out}/genoa_${type}_${class}.tsv.gz
#    head -1 ${type}.ENSG00000266714.${class}.tsv > ${out}/genoa_${type}_${class}.tsv
#    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_${type}_${class}.tsv
#    gzip ${out}/genoa_${type}_${class}.tsv
#  done
#done

# v5
for type in her corr
do
  echo v5_${type}
  rm -rf ${out}/v5_${type}.tsv.gz
  cd /scratch1/zeyunlu/sushie_v5/sushie/${type}
  head -1 ENSG00000284594_MIR7847.normal.sushie.${type}.tsv > ${out}/v5_${type}.tsv
  find . -name "*.sushie.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/v5_${type}.tsv
  gzip ${out}/v5_${type}.tsv
done

for class in cs
do
  cd /scratch1/zeyunlu/sushie_v5/sushie/${class}
  for type in  normal.sushie indep.sushie normal.meta normal.mega
  do
    echo v5_${type}_${class}
    rm -rf ${out}/v5_${type}_${class}.tsv.gz
    head -1 ENSG00000269893_SNHG8.${type}.${class}.tsv > ${out}/v5_${type}_${class}.tsv
    find . -name "*.${type}.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/v5_${type}_${class}.tsv
    gzip ${out}/v5_${type}_${class}.tsv
  done
done

for class in cs
do
  for type in susiex mesusie
  do
    cd /scratch1/zeyunlu/sushie_v5/${type}/${class}
    echo v5_${type}_${class}
    rm -rf ${out}/v5_${type}_${class}.tsv.gz
    head -1 ${type}.ENSG00000284594_MIR7847.${class}.tsv > ${out}/v5_${type}_${class}.tsv
    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/v5_${type}_${class}.tsv
    gzip ${out}/v5_${type}_${class}.tsv
  done
done

#for class in corr
#do
#  for type in xmap
#  do
#    cd /scratch1/zeyunlu/sushie_v5/${type}/${class}
#    echo v5_${type}_${class}
#    rm -rf ${out}/v5_${type}_${class}.tsv.gz
#    head -1 ${type}.ENSG00000269893_SNHG8.${class}.tsv > ${out}/v5_${type}_${class}.tsv
#    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/v5_${type}_${class}.tsv
#    gzip ${out}/v5_${type}_${class}.tsv
#  done
#done

# geuvadis
for type in her corr
do
  echo geuvadis_${type}
  rm -rf ${out}/geuvadis_${type}.tsv.gz
  cd /scratch1/zeyunlu/sushie_geuvadis/sushie/${type}
  head -1 ENSG00000263345.normal.sushie.${type}.tsv > ${out}/geuvadis_${type}.tsv
  find . -name "*.sushie.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/geuvadis_${type}.tsv
gzip  ${out}/geuvadis_${type}.tsv
done

for class in cs
do
  cd /scratch1/zeyunlu/sushie_geuvadis/sushie/${class}
  for type in normal.meta normal.mega normal.sushie indep.sushie
  do
    echo geuvadis_${type}_${class}
    rm -rf ${out}/geuvadis_${type}_${class}.tsv.gz
    head -1 ENSG00000263345.${type}.${class}.tsv >  ${out}/geuvadis_${type}_${class}.tsv
    find . -name "*.${type}.${class}.tsv" | xargs -n 1 tail -n +2 >>  ${out}/geuvadis_${type}_${class}.tsv
    gzip  ${out}/geuvadis_${type}_${class}.tsv
  done
done

for class in cs
do
  for type in susiex mesusie
  do
    cd /scratch1/zeyunlu/sushie_geuvadis/${type}/${class}
    echo geuvadis_${type}_${class}
    rm -rf ${out}/geuvadis_${type}_${class}.tsv.gz
    head -1 ${type}.ENSG00000262292.${class}.tsv > ${out}/geuvadis_${type}_${class}.tsv
    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/geuvadis_${type}_${class}.tsv
    gzip ${out}/geuvadis_${type}_${class}.tsv
  done
done

#for class in corr
#do
#  for type in xmap
#  do
#    cd /scratch1/zeyunlu/sushie_geuvadis/${type}/${class}
#    echo geuvadis_${type}_${class}
#    rm -rf ${out}/geuvadis_${type}_${class}.tsv.gz
#    head -1 ${type}.ENSG00000262292.${class}.tsv > ${out}/geuvadis_${type}_${class}.tsv
#    find . -name "${type}.*.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/geuvadis_${type}_${class}.tsv
#    gzip ${out}/geuvadis_${type}_${class}.tsv
#  done
#done

# interval
for type in her
do
  echo interval_${type}
  cd /scratch1/zeyunlu/sushie_interval/${type}
  rm -rf ${out}/interval_${type}.tsv.gz
  head -1 ENSG00000278505_Q8N4C9_C17orf78.8545.14.3.normal.sushie.${type}.tsv > ${out}/interval_${type}.tsv
  find . -name "*.sushie.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/interval_${type}.tsv
  gzip ${out}/interval_${type}.tsv
done

for class in cs
do
  cd /scratch1/zeyunlu/sushie_interval/${class}
  for type in normal.sushie
  do
    echo interval_${type}_${class}
    rm -rf ${out}/interval_${type}_${class}.tsv.gz
    head -1 ENSG00000278505_Q8N4C9_C17orf78.8545.14.3.${type}.${class}.tsv > ${out}/interval_${type}_${class}.tsv
    find . -name "*.${type}.${class}.tsv" | xargs -n 1 tail -n +2 >> ${out}/interval_${type}_${class}.tsv
    gzip ${out}/interval_${type}_${class}.tsv
  done
done


# enrichment
# proteins
cd /scratch1/zeyunlu/sushie_proteins/enrich
rm -rf ${out}/proteins_enrich_all.tsv.gz
head -1 enrich_proteins_ENSG00000115598_SL004875.tsv > ${out}/proteins_enrich_all.tsv
find . -name "enrich_proteins_*.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_enrich_all.tsv
gzip ${out}/proteins_enrich_all.tsv

# genoa
cd /scratch1/zeyunlu/sushie_genoa/enrich
rm -rf ${out}/genoa_enrich_all.tsv.gz
head -1 enrich_genoa_ENSG00000287151.tsv > ${out}/genoa_enrich_all.tsv
find . -name "enrich_genoa_*.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_enrich_all.tsv
gzip ${out}/genoa_enrich_all.tsv

# mesa rnaseq
cd /scratch1/zeyunlu/sushie_rnaseq/enrich
rm -rf ${out}/rnaseq_enrich_all.tsv.gz
head -1 enrich_rnaseq_ENSG00000284594_MIR7847.tsv > ${out}/rnaseq_enrich_all.tsv
find . -name "enrich_rnaseq_*.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_enrich_all.tsv
gzip ${out}/rnaseq_enrich_all.tsv


# tss
# proteins
cd /scratch1/zeyunlu/sushie_proteins/tss
rm -rf ${out}/proteins_tss.tsv.gz
head -1 tss.proteins.ENSG00000115598_SL004875.tsv > ${out}/proteins_tss.tsv
find . -name "tss.proteins.*.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_tss.tsv
gzip ${out}/proteins_tss.tsv

# rnaseq
cd /scratch1/zeyunlu/sushie_rnaseq/tss
rm -rf ${out}/rnaseq_tss.tsv.gz
head -1 tss.rnaseq.ENSG00000284594_MIR7847.tsv > ${out}/rnaseq_tss.tsv
find . -name "tss.rnaseq.*.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_tss.tsv
gzip ${out}/rnaseq_tss.tsv

# genoa
cd /scratch1/zeyunlu/sushie_genoa/tss
rm -rf ${out}/genoa_tss.tsv.gz
head -1 tss.genoa.ENSG00000287151.tsv > ${out}/genoa_tss.tsv
find . -name "tss.genoa.*.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_tss.tsv
gzip ${out}/genoa_tss.tsv


# valid
# proteins
cd /scratch1/zeyunlu/sushie_proteins/valid
rm -rf ${out}/proteins_valid.tsv.gz
head -1 valid.normal.ENSG00000115598_SL004875.tsv > ${out}/proteins_valid.tsv
find . -name "valid.*.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_valid.tsv
gzip ${out}/proteins_valid.tsv

# rnaseq
cd /scratch1/zeyunlu/sushie_rnaseq/valid
for type in normal indep meta mega susiex mesusie
do
  rm -rf ${out}/rnaseq_${type}_valid.tsv.gz
  head -1 valid.${type}.ENSG00000284194_SCO2.tsv > ${out}/rnaseq_${type}_valid.tsv
  find . -name "valid.${type}.*.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_${type}_valid.tsv
  gzip ${out}/rnaseq_${type}_valid.tsv
done

cd /scratch1/zeyunlu/sushie_genoa/valid
for type in mega
do
  rm -rf ${out}/genoa_${type}_valid.tsv.gz
  head -1 valid.${type}.ENSG00000250565.tsv > ${out}/genoa_${type}_valid.tsv
  find . -name "valid.${type}.*.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_${type}_valid.tsv
  gzip ${out}/genoa_${type}_valid.tsv
done

# mesa rnaseq
rm -rf ${out}/rnaseq_weights.tsv.gz
cd /scratch1/zeyunlu/sushie_rnaseq/sushie/weights
head -1 ENSG00000284594_MIR7847.normal.sushie.weights.tsv > ${out}/rnaseq_weights.tsv
for IDX in `seq 1 100`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ~/trash/case2/rnaseq_weights.tsv`
  set -- junk $params
  shift
  NAME=$1
  echo $IDX, $NAME
  find . -name "${NAME}.normal.sushie.weights.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_weights.tsv
done
gzip ${out}/rnaseq_weights.tsv

# mesa proteins
rm -rf ${out}/proteins_weights.tsv.gz
cd /scratch1/zeyunlu/sushie_proteins/sushie/weights
head -1 ENSG00000115598_SL004875.normal.sushie.weights.tsv > ${out}/proteins_weights.tsv
for IDX in `seq 1 100`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ~/trash/case2/proteins_weights.tsv`
  set -- junk $params
  shift
  NAME=$1
  echo $IDX, $NAME
  find . -name "${NAME}.normal.sushie.weights.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_weights.tsv
done
gzip ${out}/proteins_weights.tsv

# genoa
rm -rf ${out}/genoa_weights.tsv.gz
cd /scratch1/zeyunlu/sushie_genoa/sushie/weights
head -1 ENSG00000288547.normal.sushie.weights.tsv > ${out}/genoa_weights.tsv
for IDX in `seq 1 100`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ~/trash/case2/genoa_weights.tsv`
  set -- junk $params
  shift
  NAME=$1
  echo $IDX, $NAME
  find . -name "${NAME}.normal.sushie.weights.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_weights.tsv
done
gzip ${out}/genoa_weights.tsv

# twas
file=~/Documents/github/data/sushie_results/real/sushie_twas.tsv
head -1 ~/Downloads/twas_res/chr13/mesa.proteins.chr13.ENSG00000153487_SL009628.sushie.twas.tsv > ${file}
for idx in `seq 1 22`
do
  echo chr$idx
  cd ~/Downloads/twas_res/chr${idx}
  find . -name "*.tsv" | xargs -n 1 tail -n +2 >> ${file}
done
gzip ${file}


# fst
# rnaseq
cd /scratch1/zeyunlu/sushie_rnaseq/fst
for type in all_snp
do
  echo $type
  rm -rf ${out}/rnaseq_fst.${type}.tsv.gz
  head -1 ENSG00000284594_MIR7847_${type}.fst.summary > ${out}/rnaseq_fst.${type}.tsv
  find . -name "*_${type}.fst.summary" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_fst.${type}.tsv
  gzip ${out}/rnaseq_fst.${type}.tsv
done

cd /scratch1/zeyunlu/sushie_genoa/fst
for type in all_snp
do
  echo $type
  rm -rf ${out}/genoa_fst.${type}.tsv.gz
  head -1 ENSG00000287151_${type}.fst.summary > ${out}/genoa_fst.${type}.tsv
  find . -name "*_${type}.fst.summary" | xargs -n 1 tail -n +2 >> ${out}/genoa_fst.${type}.tsv
  gzip ${out}/genoa_fst.${type}.tsv
done

cd /scratch1/zeyunlu/sushie_proteins/fst
for type in all_snp
do
  echo $type
  rm -rf ${out}/proteins_fst.${type}.tsv.gz
  head -1 ENSG00000115598_SL004875_${type}.fst.summary > ${out}/proteins_fst.${type}.tsv
  find . -name "*_${type}.fst.summary" | xargs -n 1 tail -n +2 >> ${out}/proteins_fst.${type}.tsv
  gzip ${out}/proteins_fst.${type}.tsv
done

# loss
cd /scratch1/zeyunlu/sushie_rnaseq/loss
rm -rf ${out}/rnaseq_loss.tsv.gz
head -1 ENSG00000284594_MIR7847.loss.tsv > ${out}/rnaseq_loss.tsv
find . -name "*.loss.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_loss.tsv
gzip ${out}/rnaseq_loss.tsv

# r2
cd /scratch1/zeyunlu/sushie_proteins/r2
for type in sushie mesusie
do
  rm -rf ${out}/proteins_${type}.r2.tsv.gz
  head -1 ${type}.ENSG00000285441_SL001815.cvr2.tsv > ${out}/proteins_${type}.r2.tsv
  find . -name "${type}.*.cvr2.tsv" | xargs -n 1 tail -n +2 >> ${out}/proteins_${type}.r2.tsv
  gzip ${out}/proteins_${type}.r2.tsv
done

cd /scratch1/zeyunlu/sushie_rnaseq/r2
for type in mesusie sushie
do
  rm -rf ${out}/rnaseq_${type}.r2.tsv.gz
  head -1 ${type}.ENSG00000000419_DPM1.cvr2.tsv > ${out}/rnaseq_${type}.r2.tsv
  find . -name "${type}.*.cvr2.tsv" | xargs -n 1 tail -n +2 >> ${out}/rnaseq_${type}.r2.tsv
  gzip ${out}/rnaseq_${type}.r2.tsv
done

cd /scratch1/zeyunlu/sushie_genoa/r2
for type in sushie mesusie
do
  rm -rf ${out}/genoa_${type}.r2.tsv.gz
  head -1 ${type}.ENSG00000288547.cvr2.tsv > ${out}/genoa_${type}.r2.tsv
  find . -name "${type}.*.cvr2.tsv" | xargs -n 1 tail -n +2 >> ${out}/genoa_${type}.r2.tsv
  gzip ${out}/genoa_${type}.r2.tsv
done
