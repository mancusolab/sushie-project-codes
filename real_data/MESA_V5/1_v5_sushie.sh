#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-948
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 21804
# 22009 = 948 * 23

source /home1/zeyunlu/init.sh
conda activate jax2

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/scratch1/zeyunlu/sushie

EUR_COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_HIS.tsv.gz

start=`python -c "print( 1 + 23 *int(int($NR-1)))"`
stop=$((start + 22))

bigTMP=/scratch1/zeyunlu/temp_v5/tempf_${NR}

mkdir ${bigTMP}

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_rnaseq_v5_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  NAME=$2
  CHR=$4
  P0=$8
  P1=$9
  ID=${12}

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  TMPDIR=${bigTMP}/${ID}

  mkdir -p $TMPDIR

  # get genotype data
  bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis

  for pop in EUR AFR HIS
  do
    echo "Getting geno and pheno data for ${pop}"

    KEEP=${SCRATCH}/mesa_rnaseq_pt_v5_${pop}.tsv
    PHENO=${SCRATCH}/mesa_rnaseq_v5_${pop}.tsv.gz

    ${PLINK} --silent --bfile ${bfile} --chr $CHR --from-bp $P0 --to-bp $P1 \
      --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 \
      --make-bed --keep $KEEP --out $TMPDIR/${ID}.${pop}.geno

    python ${SCRATCH}/extract_pheno.py \
      --file $PHENO \
      --gene $NAME \
      --subject $KEEP \
      --out $TMPDIR/${ID}.${pop}.pheno

  done

  count=$(ls  $TMPDIR/${ID}.*.geno.bed | wc -l)
  if [ $count -ne 3 ]; then
    echo "Not all genotype data available "
    continue
  fi

  sushie finemap \
  --pheno $TMPDIR/${ID}.EUR.pheno.sushie.pheno $TMPDIR/${ID}.AFR.pheno.sushie.pheno $TMPDIR/${ID}.HIS.pheno.sushie.pheno \
  --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
  --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
  --her --alphas --numpy --meta --mega \
  --trait ${ID} \
  --output $TMPDIR/${ID}.normal

  sushie finemap \
  --pheno $TMPDIR/${ID}.EUR.pheno.sushie.pheno $TMPDIR/${ID}.AFR.pheno.sushie.pheno $TMPDIR/${ID}.HIS.pheno.sushie.pheno \
  --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
  --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
  --alphas --no-update --rho 0 0 0  --alphas --numpy \
  --trait ${ID}\
  --output $TMPDIR/${ID}.indep

  OUT=/scratch1/zeyunlu/sushie_v5

  mv $TMPDIR/${ID}.normal.sushie.corr.tsv $OUT/corr/
  mv $TMPDIR/${ID}.normal.sushie.her.tsv $OUT/her/
  mv $TMPDIR/${ID}.*.cs.tsv $OUT/cs/
  mv $TMPDIR/${ID}.*.alphas.tsv $OUT/alphas/
  mv $TMPDIR/${ID}.*.weights.tsv $OUT/weights/

done
