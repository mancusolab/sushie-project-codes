#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1700-2181
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 21804
# 21810 = 2181 * 10

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

# SCRIPTF=/home1/zeyunlu/github
OUT=/scratch1/zeyunlu/sushie_v5

SCRIPTF=/project/nmancuso_8/zeyunlu/projects

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

EUR_COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_EUR.tsv.gz
AFR_COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_AFR.tsv.gz
HIS_COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_HIS.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_v5/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 10 *int(int($NR-1)))"`
stop=$((start + 9))

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
  --her --alphas --meta --mega \
  --trait ${ID} \
  --output $TMPDIR/${ID}.normal

  # this is usually the case that we have less than 100 SNPs
  if [ ! -e "$TMPDIR/${ID}.normal.sushie.weights.tsv" ]; then
    continue
  fi

  sushie finemap \
  --pheno $TMPDIR/${ID}.EUR.pheno $TMPDIR/${ID}.AFR.pheno $TMPDIR/${ID}.HIS.pheno \
  --covar $EUR_COVAR $AFR_COVAR $HIS_COVAR \
  --plink  $TMPDIR/${ID}.EUR.geno $TMPDIR/${ID}.AFR.geno $TMPDIR/${ID}.HIS.geno \
  --no-update --rho 0 0 0 --alphas \
  --trait ${ID} \
  --output $TMPDIR/${ID}.indep

  mv $TMPDIR/${ID}.normal.sushie.corr.tsv $OUT/sushie/corr/
  mv $TMPDIR/${ID}.normal.sushie.her.tsv $OUT/sushie//her/
  mv $TMPDIR/${ID}.*.cs.tsv $OUT/sushie/cs/
  mv $TMPDIR/${ID}.*.alphas.tsv $OUT/sushie/alphas/
  cp $TMPDIR/${ID}.*.weights.tsv $OUT/sushie//weights/

  for pop in EUR AFR HIS
  do
    KEEP=${SCRATCH}/mesa_rnaseq_pt_v5_${pop}.tsv
    COVAR=${SCRATCH}/mesa_rnaseq_covar_v5_${pop}.tsv.gz

    # perform eQTLscan using plink
    ${PLINK} --bfile $TMPDIR/${ID}.${pop}.geno --keep $KEEP --glm 'hide-covar' omit-ref \
    --covar 'iid-only' $COVAR --covar-variance-standardize --pheno 'iid-only' $TMPDIR/${ID}.${pop}.pheno \
    --out ${TMPDIR}/${ID}.${pop}
  done

  # run mesusie
  python ${SCRIPTF}/sushie-data-codes/real_data/utils/run_multisusie_3pop.py \
      --ss_file ${TMPDIR}/${ID}.EUR.PHENO1.glm.linear:${TMPDIR}/${ID}.AFR.PHENO1.glm.linear:${TMPDIR}/${ID}.HIS.PHENO1.glm.linear \
      --geno_file $TMPDIR/${ID}.EUR.geno:$TMPDIR/${ID}.AFR.geno:$TMPDIR/${ID}.HIS.geno\
      --wgt_file $TMPDIR/${ID}.normal.sushie.weights.tsv \
      --rm_amb True \
      --trait ${ID} \
      --tmp_out $TMPDIR/${ID}.tmp \
      --out $OUT/multisusie

  # run in-sample susiex
  ss1=$(awk 'NR==1 {print $1}' $TMPDIR/${ID}.tmp.md.tsv)
  ss2=$(awk 'NR==2 {print $1}' $TMPDIR/${ID}.tmp.md.tsv)
  ss3=$(awk 'NR==3 {print $1}' $TMPDIR/${ID}.tmp.md.tsv)
  bp1=$(awk 'NR==1 {print $2}' $TMPDIR/${ID}.tmp.md.tsv)
  bp2=$(awk 'NR==2 {print $2}' $TMPDIR/${ID}.tmp.md.tsv)

  # run susiex
  /project/nmancuso_8/zeyunlu/tools/SuSiEx/bin_static/SuSiEx \
  --sst_file=$TMPDIR/${ID}.tmp.inss.ans0.tsv,$TMPDIR/${ID}.tmp.inss.ans1.tsv,$TMPDIR/${ID}.tmp.inss.ans2.tsv \
  --n_gwas=${ss1},${ss2},${ss3} \
  --ref_file=$TMPDIR/${ID}.EUR.geno,$TMPDIR/${ID}.AFR.geno,$TMPDIR/${ID}.HIS.geno \
  --ld_file=$TMPDIR/${ID}.EUR.geno,$TMPDIR/${ID}.AFR.geno,$TMPDIR/${ID}.HIS.geno \
  --out_dir=$TMPDIR \
  --out_name=susiex.${ID} \
  --chr=${CHR} \
  --bp=${bp1},${bp2} \
  --chr_col=1,1,1 \
  --snp_col=2,2,2 \
  --bp_col=3,3,3 \
  --a1_col=5,5,5 \
  --a2_col=4,4,4 \
  --eff_col=7,7,7 \
  --se_col=8,8,8 \
  --pval_col=9,9,9 \
  --plink=/project/nmancuso_8/zeyunlu/tools/SuSiEx/utilities/plink \
  --maf=0.01 \
  --level=0.95 \
  --n_sig=10 \
  --pval_thresh=1 \
  --max_iter=500

  # prepare susiex out
  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/process_susiex.R \
    $TMPDIR/susiex.${ID}.cs \
    $TMPDIR/susiex.${ID}.snp \
    ${ID} ${CHR} ${OUT}/susiex

  # prepare data for xmap and mesusie
  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/prepare_3pop.R \
    $TMPDIR/${ID}.tmp.inss.ans0.tsv \
    $TMPDIR/${ID}.tmp.inss.ans1.tsv \
    $TMPDIR/${ID}.tmp.inss.ans2.tsv \
    $TMPDIR/${ID}.tmp.inld.ans0.tsv \
    $TMPDIR/${ID}.tmp.inld.ans1.tsv \
    $TMPDIR/${ID}.tmp.inld.ans2.tsv \
    $TMPDIR/${ID}.tmp.inldsc.tsv \
    374:142:261 ${ID} $TMPDIR/prepare.${ID}.rdata

  # run mesusie
  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/run_mesusie_3pop.R \
    $TMPDIR/prepare.${ID}.rdata \
    $OUT/mesusie

#  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/run_xmap_3pop.R \
#    $TMPDIR/prepare.${ID}.rdata \
#    $OUT/xmap

  rm -rf ${TMPDIR}

done
