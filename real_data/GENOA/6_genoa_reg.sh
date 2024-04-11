#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --mem=8Gb
#SBATCH --array=323,324,377,395,591,594,595,727,728
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 13905
# 13905 = 927 * 15

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRIPT=/project/nmancuso_8/zeyunlu/projects/sushie-data-codes
DATA=/project/nmancuso_8/data/GENOA/sushie/
SCRATCH=/scratch1/zeyunlu/sushie

ea_pt=${DATA}/ea_373_pt.id
aa_pt=${DATA}/aa_441_pt.id

ea_covar=${SCRATCH}/ea_covars_sushie.tsv.gz
aa_covar=${SCRATCH}/aa_covars_sushie.tsv.gz

ea_pheno=${SCRATCH}/ea_pheno_sushie.tsv.gz
aa_pheno=${SCRATCH}/aa_pheno_sushie.tsv.gz

start=`python -c "print( 1 + 15 *int(int($NR-1)))"`
stop=$((start + 14))

bigTMP=/scratch1/zeyunlu/temp_pred_genoa/tempf_${NR}

mkdir ${bigTMP}

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
  ea_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/ea_chr${CHR}
  aa_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/aa_chr${CHR}

  wfile=/scratch1/zeyunlu/sushie_genoa/weights/${ID}.normal.sushie.weights.tsv
  OUT=/project/nmancuso_8/data/sushie/aou_reg/genoa.mrna/chr${CHR}/genoa.mrna.chr${CHR}.${ID}.sushie.weights.tsv
  OUT2=/scratch1/zeyunlu/sushie_genoa/r2/$ID.cvr2.tsv
  if [ ! -f ${OUT2} ]; then
    if [ -f ${wfile} ]; then
      mkdir -p $TMPDIR
      awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

      $PLINK --bfile ${ea_bfile} --chr $CHR --extract $TMPDIR/${ID}.snp --keep ${ea_pt} --snps-only \
        --rm-dup force-first --make-bed --out $TMPDIR/${ID}_ea_plink --silent
      $PLINK --bfile ${aa_bfile} --chr $CHR  --extract $TMPDIR/${ID}.snp --keep ${aa_pt} --snps-only \
        --rm-dup force-first --make-bed --out $TMPDIR/${ID}_aa_plink --silent

     python ${SCRIPT}/real_data/GENOA/extract_pheno.py \
        --file $ea_pheno \
        --gene $NAME \
        --out $TMPDIR/${ID}_ea_pheno

      python ${SCRIPT}/real_data/GENOA/extract_pheno.py \
        --file $aa_pheno \
        --gene $NAME \
        --out $TMPDIR/${ID}_aa_pheno

      python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/predict.py \
        --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed $TMPDIR/${ID}_aa_pheno.sushie.bed \
        --covar $ea_covar $aa_covar \
        --plink $TMPDIR/${ID}_ea_plink  $TMPDIR/${ID}_aa_plink \
        --trait ${ID} \
        --w_files /scratch1/zeyunlu/sushie_genoa/weights/${ID} \
        --seed $IDX \
        --out $OUT \
        --out_r2 $OUT2
    fi
  fi
done
