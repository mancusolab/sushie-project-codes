#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=1Gb
#SBATCH --array=1-302
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

bigTMP=/scratch1/zeyunlu/fst_rnaseq/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 10 *int(int($NR-1)))"`
stop=$((start + 9))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/param/real_qtl.tsv`
  set -- junk $params
  shift
  NAME=$2
  CHR=$4
  P0=$8
  P1=$9
  ID=${12}

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}

  # get genotype data
  bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis

  mkdir -p $TMPDIR

  for pop in EUR AFR HIS
  do
    KEEP=${SCRATCH}/mesa_rnaseq_pt_${pop}.tsv
    awk -v val="$pop" '{print $0 "\t" val}' $KEEP > $TMPDIR/tmp.${pop}.pt
  done

  cat $TMPDIR/tmp.*.pt > $TMPDIR/tmp.all.pt
  awk '{print "0\t"$0}' $TMPDIR/tmp.all.pt > $TMPDIR/all.pt

  ${PLINK} --bfile ${bfile}  --extract /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/param/gene_level/real_qtl_${ID}.tsv \
  --fst CATPHENO --within  $TMPDIR/all.pt --out ${TMPDIR}/${ID}_all_snp

  awk -v id="$ID" 'BEGIN{OFS=FS="\t"} NR==1{print $0, "trait"} NR>1{print $0, id}' \
    ${TMPDIR}/${ID}_all_snp.fst.summary \
    > /scratch1/zeyunlu/sushie_rnaseq/fst2/${ID}_all_snp.fst.summary

  rm -rf $TMPDIR/*

done
