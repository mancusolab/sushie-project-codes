#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=12Gb
#SBATCH --array=1-927
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
# 13905 = 927 * 15

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

# SCRIPTF=/home1/zeyunlu/github
OUT=/scratch1/zeyunlu/sushie_genoa

SCRIPTF=/project/nmancuso_8/zeyunlu/projects

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
SCRATCH=/project/nmancuso_8/data/sushie/meta_data

ea_pt=/project/nmancuso_8/data/GENOA/sushie/ea_373_pt.id
aa_pt=/project/nmancuso_8/data/GENOA/sushie/aa_441_pt.id

ea_covar=${SCRATCH}/ea_covars_sushie.tsv
aa_covar=${SCRATCH}/aa_covars_sushie.tsv

ea_covar2=${SCRATCH}/ea_covars_sushie_2iid_col.tsv
aa_covar2=${SCRATCH}/aa_covars_sushie_2iid_col.tsv

ea_pheno=${SCRATCH}/ea_pheno_sushie.tsv.gz
aa_pheno=${SCRATCH}/aa_pheno_sushie.tsv.gz

bigTMP=/scratch1/zeyunlu/temp_genoa/tempf_${NR}

mkdir -p ${bigTMP}

start=`python -c "print( 1 + 15 *int(int($NR-1)))"`
stop=$((start + 14))

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

  python ${SCRIPTF}/sushie-data-codes/real_data/GENOA/extract_pheno.py \
  --file $ea_pheno \
  --gene $NAME \
  --out $TMPDIR/${ID}_ea_pheno

  python ${SCRIPTF}/sushie-data-codes/real_data/GENOA/extract_pheno.py \
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

  wfile=/scratch1/zeyunlu/sushie_genoa/sushie/weights/${ID}.normal.sushie.weights.tsv

  # this is usually the case that we have less than 100 SNPs
  if [ ! -e "${wfile}" ]; then
    continue
  fi

#  sushie finemap \
#  --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed $TMPDIR/${ID}_aa_pheno.sushie.bed \
#  --covar $ea_covar $aa_covar\
#  --plink  $TMPDIR/${ID}_ea_plink  $TMPDIR/${ID}_aa_plink \
#  --no-update --rho 0 --alphas \
#  --trait $ID \
#  --output $TMPDIR/${ID}.indep
#
#  mv $TMPDIR/${ID}.normal.sushie.corr.tsv $OUT/sushie/corr/
#  mv $TMPDIR/${ID}.normal.sushie.her.tsv $OUT/sushie//her/
#  mv $TMPDIR/${ID}.*.cs.tsv $OUT/sushie/cs/
#  mv $TMPDIR/${ID}.*.alphas.tsv $OUT/sushie/alphas/
#  cp $TMPDIR/${ID}.*.weights.tsv $OUT/sushie//weights/

  awk 'BEGIN {FS=OFS="\t"} {print $1, $0}'  $TMPDIR/${ID}_ea_pheno.sushie.bed >  $TMPDIR/${ID}_ea_pheno.sushie.bed2
  awk 'BEGIN {FS=OFS="\t"} {print $1, $0}'  $TMPDIR/${ID}_aa_pheno.sushie.bed >  $TMPDIR/${ID}_aa_pheno.sushie.bed2

  ${PLINK} --bfile $TMPDIR/${ID}_ea_plink --keep ${ea_pt} --glm 'hide-covar' omit-ref \
    --covar  $ea_covar2 --covar-variance-standardize --pheno $TMPDIR/${ID}_ea_pheno.sushie.bed2 \
    --out ${TMPDIR}/${ID}.EUR

  ${PLINK} --bfile $TMPDIR/${ID}_aa_plink --keep ${aa_pt} --glm 'hide-covar' omit-ref \
    --covar $aa_covar2 --covar-variance-standardize --pheno $TMPDIR/${ID}_aa_pheno.sushie.bed2 \
    --out ${TMPDIR}/${ID}.AFR

  # run mesusie
  python ${SCRIPTF}/sushie-data-codes/real_data/utils/run_multisusie_2pop.py \
      --ss_file ${TMPDIR}/${ID}.EUR.PHENO1.glm.linear:${TMPDIR}/${ID}.AFR.PHENO1.glm.linear \
      --geno_file $TMPDIR/${ID}_ea_plink:$TMPDIR/${ID}_aa_plink \
      --wgt_file ${wfile} \
      --rm_amb True \
      --trait ${ID} \
      --tmp_out $TMPDIR/${ID}.tmp \
      --out $OUT/multisusie

  # run in-sample susiex
  ss1=$(awk 'NR==1 {print $1}' $TMPDIR/${ID}.tmp.md.tsv)
  ss2=$(awk 'NR==2 {print $1}' $TMPDIR/${ID}.tmp.md.tsv)
  bp1=$(awk 'NR==1 {print $2}' $TMPDIR/${ID}.tmp.md.tsv)
  bp2=$(awk 'NR==2 {print $2}' $TMPDIR/${ID}.tmp.md.tsv)

  # run susiex
  /project/nmancuso_8/zeyunlu/tools/SuSiEx/bin_static/SuSiEx \
  --sst_file=$TMPDIR/${ID}.tmp.inss.ans0.tsv,$TMPDIR/${ID}.tmp.inss.ans1.tsv \
  --n_gwas=${ss1},${ss2} \
  --ref_file=$TMPDIR/${ID}_ea_plink,$TMPDIR/${ID}_aa_plink \
  --ld_file=$TMPDIR/${ID}_ea_plink,$TMPDIR/${ID}_aa_plink \
  --out_dir=$TMPDIR \
  --out_name=susiex.${ID} \
  --chr=${CHR} \
  --bp=${bp1},${bp2} \
  --chr_col=1,1 \
  --snp_col=2,2 \
  --bp_col=3,3 \
  --a1_col=5,5 \
  --a2_col=4,4 \
  --eff_col=7,7 \
  --se_col=8,8 \
  --pval_col=9,9 \
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
  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/prepare_2pop.R \
    $TMPDIR/${ID}.tmp.inss.ans0.tsv \
    $TMPDIR/${ID}.tmp.inss.ans1.tsv \
    $TMPDIR/${ID}.tmp.inld.ans0.tsv \
    $TMPDIR/${ID}.tmp.inld.ans1.tsv \
    $TMPDIR/${ID}.tmp.inldsc.tsv \
    373:441 ${ID} $TMPDIR/prepare.${ID}.rdata

  # run xmap and mesusie
  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/run_mesusie_2pop.R \
    $TMPDIR/prepare.${ID}.rdata \
    $OUT/mesusie

#  Rscript ${SCRIPTF}/sushie-data-codes/real_data/utils/run_xmap_2pop.R \
#    $TMPDIR/prepare.${ID}.rdata \
#    $OUT/xmap

  rm -rf $TMPDIR

done
