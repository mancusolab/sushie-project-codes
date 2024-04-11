#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-60
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1500 = 60 * 25

source /home1/zeyunlu/init.sh
conda activate jax2

DATAF=/scratch1/zeyunlu/sushie_sim_genes

start=`python -c "print( 1 + 25 *int(int($NR-1)))"`
stop=$((start + 24))

for IDX in `seq $start $stop`
do
  params=`sed "${IDX}q;d" /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/param/param_3pop.tsv`
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

  echo "Running on Locus ${locus} on Row ${row}"

  OUT=/scratch1/zeyunlu/sushie_sim_3pop/sushie_sim_3pop

  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/py/run_sushie_3pop.py \
  $DATAF/EUR/EUR_gene${locus}_${gene}_snps:$DATAF/AFR/AFR_gene${locus}_${gene}_snps:$DATAF/EAS/EAS_gene${locus}_${gene}_snps \
  --N $N \
  --L1 $L1 \
  --L2 $L2 \
  --L3 $L3 \
  --h2g $h2g \
  --rho $rho \
  --sim $row \
  --locus $locus \
  --seed $seed \
  --output $OUT
done