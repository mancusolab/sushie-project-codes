#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=1Gb
#SBATCH --array=1-348
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 13905
# 13905 < 13920 = 348 * 40

source /home1/zeyunlu/init.sh
conda activate bcf

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

ea_pt=/project/nmancuso_8/data/GENOA/sushie/ea_373_pt.id
aa_pt=/project/nmancuso_8/data/GENOA/sushie/aa_441_pt.id

start=`python -c "print( 1 + 40 *int(int($NR-1)))"`
stop=$((start + 39))

bigTMP=/scratch1/zeyunlu/fst_genoa/tempf_${NR}

mkdir -p ${bigTMP}

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d" ${SCRATCH}/genoa_sushie_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  ID=$2
  NAME=$3
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${NAME}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}

  wfile=/scratch1/zeyunlu/sushie_genoa/sushie/weights/${ID}.normal.sushie.weights.tsv
  # csfile=/scratch1/zeyunlu/sushie_genoa/sushie/cs/${ID}.normal.sushie.cs.tsv

  if [ -f ${wfile} ]; then
    mkdir -p $TMPDIR

    # get genotype data
    ea_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/ea_chr${CHR}
    aa_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/aa_chr${CHR}

    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp
    # awk 'NR > 1 {print $3}' $csfile > $TMPDIR/${ID}.cs.snp

    ${PLINK} --bfile $ea_bfile --extract $TMPDIR/${ID}.snp --keep $ea_pt \
      --export vcf --out ${TMPDIR}/${ID}_vcf_geno_ea

    ${PLINK} --bfile $aa_bfile --extract $TMPDIR/${ID}.snp --keep $aa_pt \
      --export vcf --out ${TMPDIR}/${ID}_vcf_geno_aa

    bgzip ${TMPDIR}/${ID}_vcf_geno_ea.vcf
    bgzip ${TMPDIR}/${ID}_vcf_geno_aa.vcf
    bcftools index ${TMPDIR}/${ID}_vcf_geno_ea.vcf.gz
    bcftools index ${TMPDIR}/${ID}_vcf_geno_aa.vcf.gz
    bcftools merge ${TMPDIR}/${ID}_vcf_geno_ea.vcf.gz ${TMPDIR}/${ID}_vcf_geno_aa.vcf.gz \
      -o ${TMPDIR}/${ID}_vcf_geno_all.vcf.gz

    awk '{print $1"_"$2 "\tEUR"}' $ea_pt > $TMPDIR/tmp.EUR.pt
    awk '{print $1"_"$2 "\tAFR"}' $aa_pt > $TMPDIR/tmp.AFR.pt

    cat $TMPDIR/tmp.*.pt > $TMPDIR/tmp.all.pt
    awk '{print "0\t"$0}' $TMPDIR/tmp.all.pt > $TMPDIR/all.pt

    ${PLINK} --vcf ${TMPDIR}/${ID}_vcf_geno_all.vcf.gz \
    --fst CATPHENO --within $TMPDIR/all.pt --out ${TMPDIR}/${ID}_all_snp

    awk -v id="$ID" 'BEGIN{OFS=FS="\t"} NR==1{print $0, "trait"} NR>1{print $0, id}' \
      ${TMPDIR}/${ID}_all_snp.fst.summary \
      > /scratch1/zeyunlu/sushie_genoa/fst/${ID}_all_snp.fst.summary

    rm -rf $TMPDIR/*
#    if grep -q '^$' $TMPDIR/${ID}.cs.snp; then
#      echo "The empty file."
#    else
#      ${PLINK} --vcf ${TMPDIR}/${ID}_vcf_geno_all.vcf.gz --extract $TMPDIR/${ID}.cs.snp \
#        --fst CATPHENO --within $TMPDIR/all.pt --out ${TMPDIR}/${ID}_cs_snp
#      awk -v id="$ID" 'BEGIN{OFS=FS="\t"} NR==1{print $0, "trait"} NR>1{print $0, id}' \
#        ${TMPDIR}/${ID}_cs_snp.fst.summary \
#        > /scratch1/zeyunlu/sushie_genoa/fst/${ID}_cs_snp.fst.summary
#    fi
  fi
done
