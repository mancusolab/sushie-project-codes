#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=1Gb
#SBATCH --array=1-32
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1274
# 1274 < 1280 = 32 * 40

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 40 *int(int($NR-1)))"`
stop=$((start + 39))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_proteins_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  ID=${15}
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  TMPDIR=/scratch1/zeyunlu/fst_proteins/${ID}

  # get genotype data
  bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis
  wfile=/scratch1/zeyunlu/sushie_proteins/sushie/weights/${ID}.normal.sushie.weights.tsv
  # csfile=/scratch1/zeyunlu/sushie_proteins/cs/${ID}.normal.sushie.cs.tsv

  if [ -f ${wfile} ]; then
    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp
    # awk 'NR > 1 {print $3}' $csfile > $TMPDIR/${ID}.cs.snp

    for pop in EUR AFR HIS
    do
      KEEP=${SCRATCH}/mesa_proteins_pt_${pop}.tsv
      awk -v val="$pop" '{print $0 "\t" val}' $KEEP > $TMPDIR/tmp.${pop}.pt
    done

    cat $TMPDIR/tmp.*.pt > $TMPDIR/tmp.all.pt
    awk '{print "0\t"$0}' $TMPDIR/tmp.all.pt > $TMPDIR/all.pt

    ${PLINK} --bfile ${bfile}  --extract $TMPDIR/${ID}.snp \
    --fst CATPHENO --within  $TMPDIR/all.pt --out ${TMPDIR}/${ID}_all_snp

    awk -v id="$ID" 'BEGIN{OFS=FS="\t"} NR==1{print $0, "trait"} NR>1{print $0, id}' \
      ${TMPDIR}/${ID}_all_snp.fst.summary \
      > /scratch1/zeyunlu/sushie_proteins/fst/${ID}_all_snp.fst.summary

    rm -rf $TMPDIR/*
#    if grep -q '^$' $TMPDIR/${ID}.cs.snp; then
#      echo "The empty file."
#    else
#      ${PLINK} --bfile ${bfile}  --extract $TMPDIR/${ID}.cs.snp \
#        --fst CATPHENO --within  $TMPDIR/all.pt --out ${TMPDIR}/${ID}_cs_snp
#      awk -v id="$ID" 'BEGIN{OFS=FS="\t"} NR==1{print $0, "trait"} NR>1{print $0, id}' \
#        ${TMPDIR}/${ID}_cs_snp.fst.summary \
#        > /scratch1/zeyunlu/sushie_proteins/fst/${ID}_cs_snp.fst.summary
#    fi
  fi
done
