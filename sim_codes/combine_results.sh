#!/bin/bash

out=~/data/sushie/sim

cd /scratch1/zeyunlu/sushie_sim_2pop
for type in cs pip prop rho
do
  echo $type
  rm -rf ${out}/sim_2pop_${type}.tsv.gz
  head -1 sushie_sim_2pop.sim9.locus9.${type}.tsv > ${out}/sim_2pop_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_2pop_${type}.tsv
  gzip ${out}/sim_2pop_${type}.tsv
done


cd /scratch1/zeyunlu/sushie_sim_3pop
for type in cs pip prop
do
  echo $type
  rm -rf ${out}/sim_3pop_${type}.tsv.gz
  head -1 sushie_sim_3pop.sim1.locus9.${type}.tsv > ${out}/sim_3pop_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_3pop_${type}.tsv
  gzip ${out}/sim_3pop_${type}.tsv
done


cd /scratch1/zeyunlu/sushie_sim_noshared
for type in pip prop
do
  echo $type
  rm -rf ${out}/sim_noshared_${type}.tsv.gz
  head -1 sushie_sim_noshared.sim9.locus9.${type}.tsv > ${out}/sim_noshared_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_noshared_${type}.tsv
  gzip ${out}/sim_noshared_${type}.tsv
done

cd /scratch1/zeyunlu/sushie_sim_rho
for type in rho
do
  echo $type
  rm -rf ${out}/sim_rho_rho.tsv.gz
  head -1 sushie_sim_rho.sim3.locus9.${type}.tsv > ${out}/sim_rho_rho.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_rho_rho.tsv
  gzip ${out}/sim_rho_rho.tsv
done

cd /scratch1/zeyunlu/sushie_sim_pred
for type in r2 twas
do
  echo $type
  rm -rf ${out}/sim_pred_${type}.tsv.gz
  head -1 sushie_sim_pred.sim9.locus9.${type}.tsv > ${out}/sim_pred_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_pred_${type}.tsv
  gzip ${out}/sim_pred_${type}.tsv
done

cd /scratch1/zeyunlu/sushie_sim_pred3
for type in r2 twas
do
  echo $type
  rm -rf ${out}/sim_pred3_${type}.tsv.gz
  head -1 sushie_sim_pred3.sim9.locus9.${type}.tsv > ${out}/sim_pred3_${type}.tsv
  find . -name "*.${type}.tsv" | xargs -n 1 tail -n +2 >> ${out}/sim_pred3_${type}.tsv
  gzip ${out}/sim_pred3_${type}.tsv
done
