#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=24Gb
#SBATCH --array=1-220
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source /home1/zeyunlu/init.sh
conda activate jax2

module load legacy/CentOS7
module load gcc/8.3.0
module load openblas/0.3.8
module load r/4.1.0

DATAF=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/

start=`python -c "print( 1 + 25 *int(int($NR-1)))"`
stop=$((start + 24))

for IDX in `seq $start $stop`
do
  params=`sed "${IDX}q;d" /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/param/param_pred.tsv`
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
  ngwas=$8
  h2ge=$9
  locus=${10}
  seed=${11}
  gene=${12}
  chr=${13}

  echo "Running on Locus ${locus} on Row ${row}"

  TMPDIR=/scratch1/zeyunlu/tmp_pred/sim${row}_locus${locus}
  mkdir -p $TMPDIR

  OUT=/scratch1/zeyunlu/sushie_sim_pred
  tmp_output=${TMPDIR}/other
  eur1=$DATAF/new/EUR_1000G/${gene}_${chr}_geno
  afr1=$DATAF/new/AFR_1000G/${gene}_${chr}_geno

  eur2=$DATAF/new/EUR_mesa/${gene}_${chr}_geno
  afr2=$DATAF/new/AFR_mesa/${gene}_${chr}_geno

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/py/run_sushie_pred1.py \
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
  --tmp_output $tmp_output

  # run xmap and mesusie
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
    $N $L2 ${row} ${locus} ${tmp_output}.in.sim${row}.locus${locus}.rdata $L3 $L1

  # run xmap and mesusie
  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/pred_mesusie.R \
    ${tmp_output}.in.sim${row}.locus${locus}.rdata \
    ${tmp_output}.mesusie.in.sim${row}.locus${locus}.tsv

  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/pred_xmap.R \
    ${tmp_output}.in.sim${row}.locus${locus}.rdata \
    ${tmp_output}.xmap.in.sim${row}.locus${locus}.tsv

  Rscript /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/R/pred_xmap2.R \
    ${tmp_output}.in.sim${row}.locus${locus}.rdata \
    ${tmp_output}.xmap.ind.sim${row}.locus${locus}.tsv

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/py/run_sushie_pred2.py \
    ${eur1}:${afr1}:${eur2}:${afr2} \
    --N $N \
    --L1 $L1 \
    --L2 $L2 \
    --L3 $L3 \
    --h2g $h2g \
    --rho $rho \
    --ngwas $ngwas \
    --h2ge $h2ge \
    --sim $row \
    --locus $locus \
    --seed $seed \
    --sushie_file ${tmp_output}.weights.sim${row}.locus${locus}.tsv \
    --mesusie_in_file ${tmp_output}.mesusie.in.sim${row}.locus${locus}.tsv \
    --xmap_in_file ${tmp_output}.xmap.in.sim${row}.locus${locus}.tsv \
    --xmap_ind_file ${tmp_output}.xmap.ind.sim${row}.locus${locus}.tsv \
    --output $OUT/pred.sim${row}.locus${locus}

  rm -rf $TMPDIR
done