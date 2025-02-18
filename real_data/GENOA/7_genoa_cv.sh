#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=40Gb
#SBATCH --array=2389,239
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
# 13905 = 2781* 5

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

OUT=/scratch1/zeyunlu/sushie_genoa
PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data
SCRIPTF=/project/nmancuso_8/zeyunlu/projects

ea_pt=/project/nmancuso_8/data/GENOA/sushie/ea_373_pt.id
aa_pt=/project/nmancuso_8/data/GENOA/sushie/aa_441_pt.id

ea_covar=${SCRATCH}/ea_covars_sushie.tsv
aa_covar=${SCRATCH}/aa_covars_sushie.tsv

ea_covar2=${SCRATCH}/ea_covars_sushie_2iid_col.tsv
aa_covar2=${SCRATCH}/aa_covars_sushie_2iid_col.tsv

ea_pheno=${SCRATCH}/ea_pheno_sushie.tsv.gz
aa_pheno=${SCRATCH}/aa_pheno_sushie.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_cv_genoa/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 5 *int(int($NR-1)))"`
stop=$((start + 4))

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
  
  wfile=/scratch1/zeyunlu/sushie_genoa/sushie/weights/${ID}.normal.sushie.weights.tsv
  
  if [ -f ${wfile} ]; then
    TMPDIR=${bigTMP}/${ID}
    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

    ea_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/ea_chr${CHR}
    aa_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/aa_chr${CHR}

    ${PLINK} --silent --bfile ${ea_bfile} --extract $TMPDIR/${ID}.snp --rm-dup force-first \
      --make-bed --keep ${ea_pt} --out $TMPDIR/${ID}_ea_plink

    ${PLINK} --silent --bfile ${aa_bfile} --extract $TMPDIR/${ID}.snp --rm-dup force-first \
      --make-bed --keep ${aa_pt} --out $TMPDIR/${ID}_aa_plink

    python ${SCRIPTF}/sushie-data-codes/real_data/GENOA/extract_pheno.py \
      --file $ea_pheno \
      --gene $NAME \
      --out $TMPDIR/${ID}_ea_pheno

    python ${SCRIPTF}/sushie-data-codes/real_data/GENOA/extract_pheno.py \
      --file $aa_pheno \
      --gene $NAME \
      --out $TMPDIR/${ID}_aa_pheno

    python ${SCRIPTF}/sushie-data-codes/real_data/utils/predict_old.py \
      --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed $TMPDIR/${ID}_aa_pheno.sushie.bed \
      --covar $ea_covar $aa_covar \
      --plink $TMPDIR/${ID}_ea_plink $TMPDIR/${ID}_aa_plink \
      --trait ${ID} \
      --seed $((IDX + NR*70000)) \
      --out $TMPDIR \
      --out_r2 $OUT/r2/sushie.${ID}.cvr2.tsv
  
    Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/old_pred_mesusie_2pop.R \
      $TMPDIR $ID $OUT/r2/mesusie.${ID}.cvr2.tsv

    rm -rf ${TMPDIR}
  fi
done
rm -rf ${bigTMP}
