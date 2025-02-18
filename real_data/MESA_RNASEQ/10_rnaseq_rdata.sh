#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1-4371
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 21855
# 21855 = 4371 * 5

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

OUT=/scratch1/zeyunlu/sushie_rnaseq
SCRIPTF=/project/nmancuso_8/zeyunlu/projects

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

EUR_COVAR=${SCRATCH}/mesa_rnaseq_covar_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_rnaseq_covar_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_rnaseq_covar_HIS.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_rdata_rnaseq/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 5 *int(int($NR-1)))"`
stop=$((start + 4))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_rnaseq_gene_list_noMHC.tsv`
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

    KEEP=${SCRATCH}/mesa_rnaseq_pt_${pop}.tsv
    PHENO=${SCRATCH}/mesa_rnaseq_v1_${pop}.tsv.gz

    ${PLINK} --silent --bfile ${bfile} --chr $CHR --from-bp $P0 --to-bp $P1 \
      --rm-dup force-first --maf 0.01 --geno 0.1 --hwe midp 1e-6 \
      --make-bed --keep $KEEP --out $TMPDIR/${ID}.${pop}.geno

    python ${SCRIPTF}/sushie-data-codes/real_data/utils/extract_pheno.py \
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
    --pheno $TMPDIR/${ID}.EUR.pheno $TMPDIR/${ID}.AFR.pheno $TMPDIR/${ID}.HIS.pheno \
    --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
    --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
    --her --mega \
    --trait ${ID} \
    --output $TMPDIR/${ID}.normal

  if [ ! -e "$TMPDIR/${ID}.normal.sushie.weights.tsv" ]; then
    rm -rf ${TMPDIR}
    exit 0
  fi

  # compile together the results
  python ${SCRIPTF}/sushie-data-codes/real_data/utils/compile_reg.py \
    --pheno $TMPDIR/${ID}.EUR.pheno $TMPDIR/${ID}.AFR.pheno $TMPDIR/${ID}.HIS.pheno \
    --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
    --plink $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
    --trait ${ID} \
    --w_files ${TMPDIR} \
    --seed $((IDX + NR*60000)) \
    --out $TMPDIR/${ID}.reg.weights.tsv

  awk 'NR > 1 {print $3}' $TMPDIR/${ID}.normal.sushie.weights.tsv > $TMPDIR/${ID}.snp

  for pop in EUR AFR HIS
  do
    # split people
    echo "Getting geno and pheno data for ${pop}"
    KEEP=${SCRATCH}/mesa_rnaseq_pt_${pop}.tsv
    Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/split_people.R $KEEP $((IDX + NR*50000)) $TMPDIR/${ID}.${pop}.split
    PHENO=${SCRATCH}/mesa_rnaseq_v1_${pop}.tsv.gz

    for idx in `seq 1 5`
    do
      for type in "test" "train"
      do
        ${PLINK} --silent --bfile ${bfile} --extract $TMPDIR/${ID}.snp --rm-dup force-first \
          --make-bed --keep $TMPDIR/${ID}.${pop}.split.${type}.pt.group${idx}.tsv --out $TMPDIR/${ID}.${pop}.${type}.group${idx}.geno

        python ${SCRIPTF}/sushie-data-codes/real_data/utils/extract_pheno.py \
          --file $PHENO \
          --gene $NAME \
          --subject $TMPDIR/${ID}.${pop}.split.${type}.pt.group${idx}.tsv \
          --out $TMPDIR/${ID}.${pop}.${type}.group${idx}.pheno
      done
    done
  done

  for idx in `seq 1 5`
  do
    echo "################### cross validation step ${idx} ###################"
    sushie finemap \
      --pheno $TMPDIR/${ID}.EUR.train.group${idx}.pheno $TMPDIR/${ID}.AFR.train.group${idx}.pheno $TMPDIR/${ID}.HIS.train.group${idx}.pheno \
      --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
      --plink  $TMPDIR/${ID}.EUR.train.group${idx}.geno $TMPDIR/${ID}.AFR.train.group${idx}.geno $TMPDIR/${ID}.HIS.train.group${idx}.geno \
      --mega \
      --trait ${ID} \
      --output $TMPDIR/${ID}.group${idx}.normal

    if [ ! -e "$TMPDIR/${ID}.group${idx}.normal.sushie.weights.tsv" ]; then
      continue 
    fi

    # compile together the results
    python ${SCRIPTF}/sushie-data-codes/real_data/utils/compile_weights_rdata.py \
      --pheno $TMPDIR/${ID}.EUR.train.group${idx}.pheno $TMPDIR/${ID}.AFR.train.group${idx}.pheno $TMPDIR/${ID}.HIS.train.group${idx}.pheno \
      --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
      --plink  $TMPDIR/${ID}.EUR.train.group${idx}.geno $TMPDIR/${ID}.AFR.train.group${idx}.geno $TMPDIR/${ID}.HIS.train.group${idx}.geno \
      --trait ${ID} \
      --w_files ${TMPDIR} \
      --group ${idx} \
      --seed $((IDX + NR*70000)) \
      --out $TMPDIR/imp.pheno.${ID}.group${idx}
  done

  # get all the r2 results
  python ${SCRIPTF}/sushie-data-codes/real_data/utils/predict_r2_rdata.py \
    --ans 3 \
    --wk_folder ${TMPDIR} \
    --trait ${ID} \
    --out $TMPDIR/${ID}.cvr2.tsv

  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/make_FUSION.R \
    ${TMPDIR}/${ID}.reg.weights.tsv \
    ${TMPDIR}/${ID}.normal.sushie.her.tsv \
    $TMPDIR/${ID}.cvr2.tsv \
    ${ID} \
    "398:297:261" \
    /scratch1/zeyunlu/rdata/

  rm -rf ${TMPDIR}
done

rm -rf ${bigTMP}
