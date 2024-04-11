#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-927
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
SCRIPT=/scratch2/zeyunlu/projects/sushie-data-codes
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

bigTMP=/scratch1/zeyunlu/temp_genoa/tempf_${NR}

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

  mkdir -p $TMPDIR

  # get genotype data
  ea_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/ea_chr${CHR}
  aa_bfile=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155/aa_chr${CHR}

  $PLINK --bfile ${ea_bfile} --chr $CHR --from-bp $P0 --to-bp $P1 --keep ${ea_pt} --snps-only \
  --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 --make-bed --out $TMPDIR/${ID}_ea_plink --silent

  $PLINK --bfile ${aa_bfile} --chr $CHR --from-bp $P0 --to-bp $P1 --keep ${aa_pt} --snps-only \
  --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 --make-bed --out $TMPDIR/${ID}_aa_plink --silent

  # if doesn't exist something went wrong; continue
  count=$(ls  $TMPDIR/${ID}_*_plink.bed | wc -l)
  if [ $count -ne 2 ]; then
    echo "Not all genotype data available "
    continue
  fi

  python ${SCRIPT}/real_data/GENOA/extract_pheno.py \
  --file $ea_pheno \
  --gene $NAME \
  --out $TMPDIR/${ID}_ea_pheno

  python ${SCRIPT}/real_data/GENOA/extract_pheno.py \
  --file $aa_pheno \
  --gene $NAME \
  --out $TMPDIR/${ID}_aa_pheno

  sushie finemap \
  --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed $TMPDIR/${ID}_aa_pheno.sushie.bed \
  --covar $ea_covar $aa_covar\
  --plink  $TMPDIR/${ID}_ea_plink  $TMPDIR/${ID}_aa_plink \
  --her --alphas --meta --mega \
  --trait $ID \
  --output $TMPDIR/${ID}.normal

  sushie finemap \
  --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed $TMPDIR/${ID}_aa_pheno.sushie.bed \
  --covar $ea_covar $aa_covar\
  --plink  $TMPDIR/${ID}_ea_plink  $TMPDIR/${ID}_aa_plink \
  --alphas --no-update --rho 0 --alphas \
  --trait $ID \
  --output $TMPDIR/${ID}.indep

  OUT=/scratch1/zeyunlu/sushie_genoa

  mv $TMPDIR/${ID}.normal.sushie.corr.tsv $OUT/corr/
  mv $TMPDIR/${ID}.normal.sushie.her.tsv $OUT/her/
  mv $TMPDIR/${ID}.*.cs.tsv $OUT/cs/
  mv $TMPDIR/${ID}.*.alphas.tsv $OUT/alphas/
  mv $TMPDIR/${ID}.*.weights.tsv $OUT/weights/

done
