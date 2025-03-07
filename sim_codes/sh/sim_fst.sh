#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=4Gb
#SBATCH --array=1-1350
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

DATAF=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
PLINK1=/project/nmancuso_8/zeyunlu/tools/plink/plink

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

  eur1=$DATAF/new/EUR_1000G/${gene}_${chr}_geno
  afr1=$DATAF/new/AFR_1000G/${gene}_${chr}_geno

  eur2=$DATAF/new/EUR_mesa/${gene}_${chr}_geno
  afr2=$DATAF/new/AFR_mesa/${gene}_${chr}_geno
  
  TMPDIR=/scratch1/zeyunlu/tmp_fst/sim${row}_locus${locus}/
  mkdir -p $TMPDIR
  
  python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_codes/py/prepare_fst.py \
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
  --tmp_output ${TMPDIR}/snp
  
  ${PLINK1} --bfile ${eur1} --bmerge ${afr1} --make-bed --out ${TMPDIR}/merged

  awk -F'\t' 'OFS="\t" {print $1, $2, "EUR"}' ${eur1}.fam > ${TMPDIR}/all.pt
  awk -F'\t' 'OFS="\t" {print $1, $2, "AFR"}' ${afr1}.fam >> ${TMPDIR}/all.pt

  ${PLINK} --bfile ${TMPDIR}/merged  --extract ${TMPDIR}/snp.name \
  --fst CATPHENO --within  $TMPDIR/all.pt --out ${TMPDIR}/fst

  awk -F'\t' -v r="$row" -v l="$locus" 'NR==2 {print $3, r, l}' OFS='\t' ${TMPDIR}/fst.fst.summary > /scratch1/zeyunlu/sushie_sim_fst/sim${row}.locus${locus}.causal.fst.tsv
  cp ${TMPDIR}/snp.causal.tsv /scratch1/zeyunlu/sushie_sim_fst/sim${row}.locus${locus}.causal.ld.tsv
  rm -rf $TMPDIR
done
  