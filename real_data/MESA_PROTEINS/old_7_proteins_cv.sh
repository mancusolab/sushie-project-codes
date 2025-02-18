#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1-36
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1274
# 535

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

OUT=/scratch1/zeyunlu/sushie_proteins
SCRIPTF=/project/nmancuso_8/zeyunlu/projects

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

EUR_COVAR=${SCRATCH}/mesa_proteins_covar_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_proteins_covar_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_proteins_covar_HIS.tsv.gz

# start=`python -c "print( 1 + 5 *int(int($NR-1)))"`
# stop=$((start + 4))

for IDX in `seq $NR $NR`
do
  # Extract parameters for sushie run
  # params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_proteins_cs_trait_list.tsv`
  params=`sed "${IDX}q;d"  /scratch1/zeyunlu/fish/proteins_fish.tsv`
  set -- junk $params
  shift
  ID=${15}
  CHR=$4
  P0=$8
  P1=$9

  echo SUMMARY ${NR}, ${IDX}, ${ID}, ${CHR}, ${P0}, ${P1}

  fr2=/scratch1/zeyunlu/sushie_proteins/r2/${ID}.cvr2.tsv
  if [ ! -f ${fr2} ]; then

    TMPDIR=/scratch1/zeyunlu/temp_cv_proteins/${ID}
    # get genotype data
    bfile=/project/nmancuso_8/data/sushie/plink/TOPMED.WGS.chr${CHR}.filtered.analysis
    wfile=/scratch1/zeyunlu/sushie_proteins/sushie/weights/${ID}.normal.sushie.weights.tsv

    mkdir -p $TMPDIR
    awk 'NR > 1 {print $3}' $wfile > $TMPDIR/${ID}.snp

    for pop in EUR AFR HIS
    do
      # split people
      echo "Getting geno and pheno data for ${pop}"
      KEEP=${SCRATCH}/mesa_proteins_pt_${pop}.tsv
      Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/split_people.R $KEEP $((IDX + NR*10000)) $TMPDIR/${ID}.${pop}.split
      PHENO=${SCRATCH}/mesa_proteins_v1_${pop}.tsv.gz

      for idx in `seq 1 5`
      do
        for type in "test" "train"
        do
          ${PLINK} --silent --bfile ${bfile} --extract $TMPDIR/${ID}.snp --rm-dup force-first \
            --make-bed --keep $TMPDIR/${ID}.${pop}.split.${type}.pt.group${idx}.tsv --out $TMPDIR/${ID}.${pop}.${type}.group${idx}.geno

          python ${SCRIPTF}/sushie-data-codes/real_data/utils/extract_pheno.py \
            --file $PHENO \
            --gene $ID \
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
        --mega --meta \
        --trait ${ID} \
        --output $TMPDIR/${ID}.group${idx}.normal

        if [ ! -e "$TMPDIR/${ID}.group${idx}.normal.sushie.weights.tsv" ]; then
          rm -rf ${TMPDIR}
          continue 2
        fi

        sushie finemap \
          --pheno $TMPDIR/${ID}.EUR.train.group${idx}.pheno $TMPDIR/${ID}.AFR.train.group${idx}.pheno $TMPDIR/${ID}.HIS.train.group${idx}.pheno \
          --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
          --plink  $TMPDIR/${ID}.EUR.train.group${idx}.geno $TMPDIR/${ID}.AFR.train.group${idx}.geno $TMPDIR/${ID}.HIS.train.group${idx}.geno \
          --no-update --rho 0 0 0 \
          --trait ${ID} \
          --output $TMPDIR/${ID}.group${idx}.indep

        for pop in EUR AFR HIS
        do
          KEEP=$TMPDIR/${ID}.${pop}.split.train.pt.group${idx}.tsv
          COVAR=${SCRATCH}/mesa_proteins_covar_${pop}.tsv.gz

          # perform eQTLscan using plink
          ${PLINK} --bfile $TMPDIR/${ID}.${pop}.train.group${idx}.geno --keep $KEEP --glm 'hide-covar' omit-ref \
            --covar 'iid-only' $COVAR --covar-variance-standardize --pheno 'iid-only' $TMPDIR/${ID}.${pop}.train.group${idx}.pheno \
            --out ${TMPDIR}/${ID}.${pop}.train.group${idx}
        done

        file1=${TMPDIR}/${ID}.EUR.train.group${idx}.PHENO1.glm.linear
        file2=${TMPDIR}/${ID}.AFR.train.group${idx}.PHENO1.glm.linear
        file3=${TMPDIR}/${ID}.HIS.train.group${idx}.PHENO1.glm.linear

        if [[ ! -f "$file1" || ! -f "$file2" || ! -f "$file3" ]]; then
          continue
        fi

        # run mesusie
        python ${SCRIPTF}/sushie-data-codes/real_data/utils/pred_multisusie_3pop.py \
          --ss_file ${TMPDIR}/${ID}.EUR.train.group${idx}.PHENO1.glm.linear:${TMPDIR}/${ID}.AFR.train.group${idx}.PHENO1.glm.linear:${TMPDIR}/${ID}.HIS.train.group${idx}.PHENO1.glm.linear \
          --geno_file $TMPDIR/${ID}.EUR.train.group${idx}.geno:$TMPDIR/${ID}.AFR.train.group${idx}.geno:$TMPDIR/${ID}.HIS.train.group${idx}.geno \
          --wgt_file $TMPDIR/${ID}.group${idx}.normal.sushie.weights.tsv \
          --rm_amb True \
          --trait ${ID} \
          --group ${idx} \
          --tmp_out $TMPDIR/${ID}.tmp

        file1="$TMPDIR/${ID}.tmp.group${idx}.inss.ans0.tsv"
        file2="$TMPDIR/${ID}.tmp.group${idx}.inss.ans1.tsv"
        file3="$TMPDIR/${ID}.tmp.group${idx}.inss.ans2.tsv"
        file4="$TMPDIR/${ID}.tmp.group${idx}.inld.ans0.tsv"
        file5="$TMPDIR/${ID}.tmp.group${idx}.inld.ans1.tsv"
        file6="$TMPDIR/${ID}.tmp.group${idx}.inld.ans2.tsv"
        file7="$TMPDIR/${ID}.tmp.group${idx}.inldsc.tsv"

        if [[ ! -f "$file1" || ! -f "$file2" || ! -f "$file3" || ! -f "$file4" || ! -f "$file5" || ! -f "$file6" || ! -f "$file7" ]]; then
          continue
        fi

        # prepare data for mesusie
        Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/prepare_3pop.R \
          $TMPDIR/${ID}.tmp.group${idx}.inss.ans0.tsv \
          $TMPDIR/${ID}.tmp.group${idx}.inss.ans1.tsv \
          $TMPDIR/${ID}.tmp.group${idx}.inss.ans2.tsv \
          $TMPDIR/${ID}.tmp.group${idx}.inld.ans0.tsv \
          $TMPDIR/${ID}.tmp.group${idx}.inld.ans1.tsv \
          $TMPDIR/${ID}.tmp.group${idx}.inld.ans2.tsv \
          $TMPDIR/${ID}.tmp.group${idx}.inldsc.tsv \
          402:175:277 ${ID} $TMPDIR/prepare.${ID}.group${idx}.rdata

        # run xmap and mesusie
        Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/pred_mesusie_3pop.R \
          $TMPDIR/prepare.${ID}.group${idx}.rdata \
          $TMPDIR/mesusie.${ID}.group${idx}.tmp.weights.tsv

        if [[ ! -f "$TMPDIR/mesusie.${ID}.group${idx}.tmp.weights.tsv" ]]; then
          continue
        fi

        # compile together the results
        python ${SCRIPTF}/sushie-data-codes/real_data/utils/compile_weights.py \
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
      python ${SCRIPTF}/sushie-data-codes/real_data/utils/predict_r2.py \
        --ans 3 \
        --wk_folder ${TMPDIR} \
        --trait ${ID} \
        --out $TMPDIR/${ID}.cvr2.tsv

    cp $TMPDIR/${ID}.cvr2.tsv ${OUT}/r2/

    rm -rf ${TMPDIR}

  fi
done
