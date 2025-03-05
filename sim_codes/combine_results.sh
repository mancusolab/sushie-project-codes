#!/bin/bash
salloc --time=12:00:00 --mem=16G --partition=conti
out=~/data/sushie/sim3

cd /scratch1/zeyunlu/sushie_sim_2pop/sushie
for type in pip
do
  echo $type
  rm -rf ${out}/sushie_2pop_${type}.tsv.gz
  head -1 sushie_sim_2pop.sim4.locus9.${type}.tsv > ${out}/sushie_2pop_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sushie_2pop_${type}.tsv
  gzip ${out}/sushie_2pop_${type}.tsv
done

cd /scratch1/zeyunlu/sushie_sim_2pop/susiex
for ld in "in"
do
  for type in cs pip sens
  do
    echo $type $ld
    rm -rf ${out}/susiex_${ld}_2pop_${type}.tsv.gz
    head -1 susiex.${ld}.sim7.locus9.${type}.tsv >${out}/susiex_${ld}_2pop_${type}.tsv
    find . -name "*.${ld}.*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/susiex_${ld}_2pop_${type}.tsv
    gzip ${out}/susiex_${ld}_2pop_${type}.tsv
  done
done

cd /scratch1/zeyunlu/sushie_sim_2pop/mesusie
for ld in "in"
do
  for type in cs pip sens
  do
    echo $type $ld
    rm -rf ${out}/mesusie_${ld}_2pop_${type}.tsv.gz
    head -1 mesusie.${ld}.sim7.locus9.${type}.tsv > ${out}/mesusie_${ld}_2pop_${type}.tsv
    find . -name "*.${ld}.*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/mesusie_${ld}_2pop_${type}.tsv
    gzip ${out}/mesusie_${ld}_2pop_${type}.tsv
  done
done

cd /scratch1/zeyunlu/sushie_sim_2pop/xmap
# check which files selected as initial
for idx in `seq 1 500`
do
  for jdx in `seq 1 7`
  do
      echo sim${jdx}.locus${idx}
      ll *sim${jdx}.locus${idx}.* | wc -l
  done
done

for ld in "in" "ind"
do
  for type in cs pip rho sens
  do
    echo $type $ld
    rm -rf ${out}/xmap_${ld}_2pop_${type}.tsv.gz
    head -1 xmap.${ld}.sim6.locus2.${type}.tsv > ${out}/xmap_${ld}_2pop_${type}.tsv
    find . -name "*.${ld}.*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/xmap_${ld}_2pop_${type}.tsv
    gzip ${out}/xmap_${ld}_2pop_${type}.tsv
  done
done


cd /scratch1/zeyunlu/sushie_sim_noshared/sushie
for type in fdr rho pip
do
  echo $type
  rm -rf ${out}/sushie_noshared_${type}.tsv.gz
  head -1 sushie.sim4.locus9.${type}.tsv > ${out}/sushie_noshared_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sushie_noshared_${type}.tsv
  gzip ${out}/sushie_noshared_${type}.tsv
done

cd /scratch1/zeyunlu/sushie_sim_noshared/susiex
for type in fdr pip
do
  echo $type
  rm -rf ${out}/susiex.in_noshared_${type}.tsv.gz
  head -1 susiex.in.sim4.locus9.${type}.tsv >${out}/susiex.in_noshared_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/susiex.in_noshared_${type}.tsv
  gzip ${out}/susiex.in_noshared_${type}.tsv
done

cd /scratch1/zeyunlu/sushie_sim_noshared/xmap
# check which files selected as initial
for idx in `seq 1 500`
do
  for jdx in `seq 1 9`
  do
      echo sim${jdx}.locus${idx}
      ll *sim${jdx}.locus${idx}.* | wc -l
  done
done

cd /scratch1/zeyunlu/sushie_sim_noshared/xmap
for method in "in" "ind"
do
  for type in fdr rho pip
  do
    echo $type ${method}
    rm -rf ${out}/xmap_${method}_noshared_${type}.tsv.gz
    head -1 xmap.${method}.sim3.locus1.${type}.tsv > ${out}/xmap_${method}_noshared_${type}.tsv
    find . -name "xmap.${method}.*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/xmap_${method}_noshared_${type}.tsv
    gzip ${out}/xmap_${method}_noshared_${type}.tsv
  done
done

cd /scratch1/zeyunlu/sushie_sim_noshared/mesusie
for method in "in"
do
  for type in fdr pip
  do
    echo $type ${method}
    rm -rf ${out}/mesusie_${method}_noshared_${type}.tsv.gz
    head -1 mesusie.${method}.sim4.locus9.${type}.tsv >  ${out}/mesusie_${method}_noshared_${type}.tsv
    find . -name "mesusie.${method}.*.${type}.tsv" | xargs -n 1 tail -n +2 >>  ${out}/mesusie_${method}_noshared_${type}.tsv
    gzip ${out}/mesusie_${method}_noshared_${type}.tsv
  done
done

cd /scratch1/zeyunlu/sushie_sim_pred
for type in r2 twas
do
  echo $type
  rm -rf ${out}/sim_pred_${type}.tsv.gz
  head -1 pred.sim5.locus9.${type}.tsv > ${out}/sim_pred_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_pred_${type}.tsv
  gzip ${out}/sim_pred_${type}.tsv
done


cd /scratch1/zeyunlu/sushie_sim_3pop
for type in cs pip
do
  echo $type
  rm -rf ${out}/sim_3pop_${type}.tsv.gz
  head -1 sushie_sim_3pop.sim3.locus9.${type}.tsv > ${out}/sim_3pop_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_3pop_${type}.tsv
  gzip ${out}/sim_3pop_${type}.tsv
done

# new simulations
out=~/data/sushie/sim3
cd /scratch1/zeyunlu/rho_files/
for type in rho
do
  echo $type
  rm -rf ${out}/sushie_2pop_${type}_all.tsv.gz
  head -1 sim10.locus1.${type}.tsv > ${out}/sushie_2pop_${type}_all.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sushie_2pop_${type}_all.tsv
  gzip ${out}/sushie_2pop_${type}_all.tsv
done

out=~/data/sushie/sim3
cd /scratch1/zeyunlu/sushie_sim_fst/

for type in fst ld
do
  rm -rf ${out}/sim_causal_${type}.tsv.gz
  cat *.${type}.tsv > ${out}/sim_causal_${type}.tsv
  gzip ${out}/sim_causal_${type}.tsv
done


out=~/data/sushie/sim3
cd /scratch1/zeyunlu/tmp_gene_fst/

rm -rf ${out}/sim_gene_fst.tsv.gz
cat *.tsv > ${out}/sim_gene_fst.tsv
gzip ${out}/sim_gene_fst.tsv



out=~/data/sushie/sim3
cd /scratch1/zeyunlu/sushie_rnaseq/fst2

echo $type
rm -rf ${out}/sushie_real_fst.tsv.gz
head -1 ENSG00000284594_MIR7847_all_snp.fst.summary > ${out}/sushie_real_fst.tsv
find . -name "*.summary" | xargs -n 1 tail -n +2 >> ${out}/sushie_real_fst.tsv
gzip ${out}/sushie_real_fst.tsv



