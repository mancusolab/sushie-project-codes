#!/bin/bash
#SBATCH --time=6:30:00
#SBATCH --mem=4Gb
#SBATCH --array=1-2969
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 5938
# 5938 = 2969 * 2

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRIPTF=/project/nmancuso_8/zeyunlu/projects
DATA=/project/nmancuso_8/data/GENOA/sushie/
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

ea_pt=${DATA}/ea_373_pt.id
aa_pt=${DATA}/aa_441_pt.id

ea_covar=${SCRATCH}/ea_covars_sushie.tsv
aa_covar=${SCRATCH}/aa_covars_sushie.tsv

ea_pheno=${SCRATCH}/ea_pheno_sushie.tsv.gz
aa_pheno=${SCRATCH}/aa_pheno_sushie.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_pred_genoa/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 2 *int(int($NR-1)))"`
stop=$((start + 1))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d" ${SCRATCH}/genoa_cs_trait_list.tsv`
  set -- junk $params
  shift
  ID=$2
  NAME=$3
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${NAME}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}
  ea_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/ea_chr${CHR}
  aa_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/aa_chr${CHR}

  wfile=/scratch1/zeyunlu/sushie_genoa/sushie/weights/${ID}.normal.sushie.weights.tsv
  OUT=/project/nmancuso_8/data/sushie/aou_reg/
  if [ -f ${wfile} ]; then
    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

    $PLINK --bfile ${ea_bfile} --chr $CHR --extract $TMPDIR/${ID}.snp --keep ${ea_pt} --snps-only \
      --rm-dup force-first --make-bed --out $TMPDIR/${ID}_ea_plink --silent
    $PLINK --bfile ${aa_bfile} --chr $CHR  --extract $TMPDIR/${ID}.snp --keep ${aa_pt} --snps-only \
      --rm-dup force-first --make-bed --out $TMPDIR/${ID}_aa_plink --silent

   python ${SCRIPTF}/sushie-data-codes/real_data/GENOA/extract_pheno.py \
      --file $ea_pheno \
      --gene $NAME \
      --out $TMPDIR/${ID}_ea_pheno

    python ${SCRIPTF}/sushie-data-codes/real_data/GENOA/extract_pheno.py \
      --file $aa_pheno \
      --gene $NAME \
      --out $TMPDIR/${ID}_aa_pheno

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict_aou.py \
      --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed $TMPDIR/${ID}_aa_pheno.sushie.bed \
      --covar $ea_covar $aa_covar \
      --plink $TMPDIR/${ID}_ea_plink  $TMPDIR/${ID}_aa_plink \
      --trait ${ID} \
      --w_files /scratch1/zeyunlu/sushie_genoa/ \
      --study genoa.mrna \
      --seed $((IDX + NR*30000)) \
      --out $OUT

    rm -rf ${TMPDIR}
  fi
done
rm -rf ${bigTMP}
