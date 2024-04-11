#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-133
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 3189
# 3192 = 133 * 24

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/scratch1/zeyunlu/sushie

COVAR=${SCRATCH}/interval_covar.tsv

start=`python -c "print( 1 + 24 *int(int($NR-1)))"`
stop=$((start + 23))

bigTMP=/scratch1/zeyunlu/temp_interval/tempf_${NR}

mkdir ${bigTMP}

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/interval_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  CHR=$1
  P0=$5
  P1=$6
  NAME=${12}
  ID=${13}

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${NAME}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}

  mkdir -p $TMPDIR

  # get genotype data
  bfile=/project/nmancuso_8/data/INTERVAL/genotype/all_final_split/INTERVAL_SOMALOGIC_POSTQC_chr${CHR}_final_v1.dbsnp153

  GENO_KEEP=${SCRATCH}/interval_pt_geno.tsv
  PHENO_KEEP=${SCRATCH}/interval_pt_pheno.tsv
  PHENO=${SCRATCH}/interval_protein_levels.tsv

  ${PLINK} --silent --bfile ${bfile} --chr $CHR --from-bp $P0 --to-bp $P1 \
    --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 \
    --make-bed --keep $GENO_KEEP --out $TMPDIR/${ID}.geno

  python ${SCRATCH}/extract_pheno.py \
    --file $PHENO \
    --gene $NAME \
    --subject $PHENO_KEEP \
    --out $TMPDIR/${ID}.pheno


  count=$(ls  $TMPDIR/${ID}.geno.bed | wc -l)
  if [ $count -ne 1 ]; then
    echo "no genotype data available "
    continue
  fi

  sushie finemap \
  --pheno $TMPDIR/${ID}.pheno.sushie.pheno \
  --covar $COVAR \
  --plink $TMPDIR/${ID}.geno \
  --her --alphas --numpy  \
  --trait ${ID} \
  --output $TMPDIR/${ID}.normal

  OUT=/scratch1/zeyunlu/sushie_interval

  mv $TMPDIR/${ID}.normal.sushie.her.tsv $OUT/her/
  mv $TMPDIR/${ID}.*.cs.tsv $OUT/cs/
  mv $TMPDIR/${ID}.*.alphas.tsv $OUT/alphas/
  mv $TMPDIR/${ID}.*.weights.tsv $OUT/weights/

done
