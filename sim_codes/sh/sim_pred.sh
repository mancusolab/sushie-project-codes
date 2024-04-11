#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=24Gb
#SBATCH --array=1-220
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 5500 = 220 * 25

source /home1/zeyunlu/init.sh
conda activate jax2

DATAF=/scratch1/zeyunlu/sushie_sim_genes

start=`python -c "print( 1 + 25 *int(int($NR-1)))"`
stop=$((start + 24))

for IDX in `seq $start $stop`
do
  params=`sed "${IDX}q;d" /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/param/param_pred.tsv`
  echo "NR=${IDX}: ${params}"
  set -- junk $params
  shift
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

  echo "Running on Locus ${locus} on Row ${row}"

  OUT=/scratch1/zeyunlu/sushie_sim_pred3/sushie_sim_pred3

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/py/run_sushie_pred.py \
  $DATAF/EUR/EUR_gene${locus}_${gene}_snps:$DATAF/AFR/AFR_gene${locus}_${gene}_snps \
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
  --output $OUT
done