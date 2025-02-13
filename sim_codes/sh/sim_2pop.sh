#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-1300
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 13000 = 520 * 25

source /home1/zeyunlu/init.sh
conda activate jax2

module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

DATAF=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/

start=`python -c "print( 1 + 10 *int(int($NR-1)))"`
stop=$((start + 9))

for IDX in `seq $start $stop`
do
  params=`sed "${IDX}q;d" /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/param/param_2pop.tsv`
  echo "NR=${IDX}: ${params}"
  set -- junk $params
  shift

  # PARAMETERS
  row=$1
  N=$2
  L1=$3
  L2=$4
  L3=$5
  h2g=$6
  rho=$7
  locus=$8
  seed=$9
  gene=${10}
  chr=${11}

  echo "Running on Locus ${locus} on Row ${row}"

  TMPDIR=/scratch1/zeyunlu/tmp_2pop/sim${row}_locus${locus}

  mkdir -p $TMPDIR

  OUT=/scratch1/zeyunlu/sushie_sim_2pop
  tmp_output=${TMPDIR}/other
  eur1=$DATAF/new/EUR_1000G/${gene}_${chr}_geno
  afr1=$DATAF/new/AFR_1000G/${gene}_${chr}_geno

  eur2=$DATAF/new/EUR_mesa/${gene}_${chr}_geno
  afr2=$DATAF/new/AFR_mesa/${gene}_${chr}_geno

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/py/run_sushie_2pop.py \
  ${eur1}:${afr1}:${eur2}:${afr2} \
  --N $N \
  --L1 $L1 \
  --L2 $L2 \
  --L3 $L3 \
  --h2g $h2g \
  --rho $rho \
  --sim $row \
  --locus $locus \
  --seed $seed \
  --tmp_output $tmp_output \
  --output $OUT/sushie/sushie_sim_2pop

  # run in-sample susiex
  ss1=$(awk 'NR==1 {print $1}' ${tmp_output}.md.sim${row}.locus${locus}.tsv)
  ss2=$(awk 'NR==2 {print $1}' ${tmp_output}.md.sim${row}.locus${locus}.tsv)
  bp1=$(awk 'NR==1 {print $2}' ${tmp_output}.md.sim${row}.locus${locus}.tsv)
  bp2=$(awk 'NR==2 {print $2}' ${tmp_output}.md.sim${row}.locus${locus}.tsv)

  /project/nmancuso_8/zeyunlu/tools/SuSiEx/bin_static/SuSiEx \
    --sst_file=${tmp_output}.inss.sim${row}.locus${locus}.ans0.tsv,${tmp_output}.inss.sim${row}.locus${locus}.ans1.tsv \
    --n_gwas=${ss1},${ss2} \
    --ref_file=${eur1},${afr1} \
    --ld_file=${eur1},${afr1} \
    --out_dir=$TMPDIR \
    --out_name=susiex.in.sim${row}.locus${locus} \
    --chr=${chr} \
    --bp=${bp1},${bp2} \
    --chr_col=1,1 \
    --snp_col=2,2 \
    --bp_col=3,3 \
    --a1_col=5,5 \
    --a2_col=4,4 \
    --eff_col=6,6 \
    --se_col=7,7 \
    --pval_col=8,8 \
    --plink=/project/nmancuso_8/zeyunlu/tools/SuSiEx/utilities/plink \
    --maf=0.01 \
    --level=0.95 \
    --n_sig=$L2 \
    --pval_thresh=1 \
    --max_iter=500

  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/process_susiex_2pop.R \
    ${TMPDIR} ${row} ${locus} ${tmp_output}.causal.sim${row}.locus${locus}.tsv ${OUT}/susiex/susiex $L2

  # prepare data for xmap and mesusie
  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/prepare_data.R \
    ${tmp_output}.inss.sim${row}.locus${locus}.ans0.tsv \
    ${tmp_output}.inss.sim${row}.locus${locus}.ans1.tsv \
    ${tmp_output}.inld.sim${row}.locus${locus}.ans0.tsv \
    ${tmp_output}.inld.sim${row}.locus${locus}.ans1.tsv \
    ${tmp_output}.xy.sim${row}.locus${locus}.ans0.tsv \
    ${tmp_output}.xy.sim${row}.locus${locus}.ans1.tsv \
    ${tmp_output}.inldsc.sim${row}.locus${locus}.tsv \
    ${tmp_output}.causal.sim${row}.locus${locus}.tsv \
    $N $L2 ${row} ${locus} ${tmp_output}.in.sim${row}.locus${locus}.rdata

   run xmap and mesusie
  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/run_mesusie.R \
    ${tmp_output}.in.sim${row}.locus${locus}.rdata \
    $OUT/mesusie/mesusie.in.sim${row}.locus${locus}
#
  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/run_xmap.R \
    ${tmp_output}.in.sim${row}.locus${locus}.rdata \
    $OUT/xmap/xmap.in.sim${row}.locus${locus}

  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/run_xmap2.R \
    ${tmp_output}.in.sim${row}.locus${locus}.rdata \
    $OUT/xmap/xmap.ind.sim${row}.locus${locus}

  rm -rf $TMPDIR
done
