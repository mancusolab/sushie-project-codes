#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-843
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 20228
# 20232 = 843 * 24

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/scratch1/zeyunlu/sushie

EUR_COVAR=${SCRATCH}/geuvadis_EUR_covar.tsv
YRI_COVAR=${SCRATCH}/geuvadis_YRI_covar.tsv

start=`python -c "print( 1 + 24 *int(int($NR-1)))"`
stop=$((start + 23))

bigTMP=/scratch1/zeyunlu/temp_geuvadis/tempf_${NR}

mkdir ${bigTMP}

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/geuvadis_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  CHR=$1
  P0=$5
  P1=$6
  ID=$8

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}

  mkdir -p $TMPDIR

  # get genotype data
  bfile=/project/nmancuso_8/data/GEUVADIS/GenotypeData/plink/raw/GEUVADIS.chr${CHR}.annotated.dbSNP_v153

  for pop in EUR YRI
  do
    echo "Getting geno and pheno data for ${pop}"

    KEEP=${SCRATCH}/geuvadis_${pop}_pt.tsv
    PHENO=${SCRATCH}/geuvadis_${pop}_expression_rpkm.tsv

    ${PLINK} --silent --bfile ${bfile} --chr $CHR --from-bp $P0 --to-bp $P1 \
      --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 \
      --make-bed --keep $KEEP --out $TMPDIR/${ID}.${pop}.geno

    python ${SCRATCH}/extract_pheno.py \
      --file $PHENO \
      --gene $ID \
      --subject $KEEP \
      --out $TMPDIR/${ID}.${pop}.pheno
  done

  count=$(ls  $TMPDIR/${ID}.*.geno.bed | wc -l)
  if [ $count -ne 2 ]; then
    echo "Not all genotype data available "
    continue
  fi

  sushie finemap \
  --pheno $TMPDIR/${ID}.EUR.pheno.sushie.pheno $TMPDIR/${ID}.YRI.pheno.sushie.pheno \
  --covar $EUR_COVAR $YRI_COVAR\
  --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.YRI.geno \
  --her --alphas --numpy --meta --mega  \
  --trait ${ID} \
  --output $TMPDIR/${ID}.normal

  sushie finemap \
  --pheno $TMPDIR/${ID}.EUR.pheno.sushie.pheno $TMPDIR/${ID}.YRI.pheno.sushie.pheno \
  --covar $EUR_COVAR $YRI_COVAR \
  --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.YRI.geno \
  --alphas --no-update --rho 0 --alphas --numpy \
  --trait ${ID}\
  --output $TMPDIR/${ID}.indep

  OUT=/scratch1/zeyunlu/sushie_geuvadis

  mv $TMPDIR/${ID}.normal.sushie.corr.tsv $OUT/corr/
  mv $TMPDIR/${ID}.normal.sushie.her.tsv $OUT/her/
  mv $TMPDIR/${ID}.*.cs.tsv $OUT/cs/
  mv $TMPDIR/${ID}.*.alphas.tsv $OUT/alphas/
  mv $TMPDIR/${ID}.*.weights.tsv $OUT/weights/

done
