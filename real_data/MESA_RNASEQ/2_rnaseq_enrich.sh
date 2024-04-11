#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-4371
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

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

SCRATCH=/scratch1/zeyunlu/sushie

start=`python -c "print( 1 + 5 *int(int($NR-1)))"`
stop=$((start + 4))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_rnaseq_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  CHR=$4
  ID=${12}
  if find "/scratch1/zeyunlu/sushie_rnaseq/alphas/" -maxdepth 1 -name "${ID}.*.alphas.tsv" -print -quit | grep -q .;
  then
    TMPDIR=/scratch1/zeyunlu/trash_rnaseq/${ID}
    mkdir $TMPDIR
    cd $TMPDIR
    for JDX in `seq 77`
    do
      params=`sed "${JDX}q;d" /scratch1/zeyunlu/new_annot/mesa_anno_list.tsv`
      set -- junk $params
      shift
      PRE=$1
      SUF=$2

      echo SUMMARY ${IDX}, ${JDX}, ${ID}, ${CHR}, $PRE, $SUF

      python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/calc_enrich_all.py \
      --alpha-path /scratch1/zeyunlu/sushie_rnaseq/alphas/ \
      --anno-file /scratch1/zeyunlu/new_annot/processed/mesa/${PRE}.chr${CHR}.${SUF} \
      --gene $ID \
      --exclude /scratch1/zeyunlu/new_annot/anno_exclude.tsv \
      --enrich-out $TMPDIR/enrich.${ID}.${PRE}.${SUF}.tsv

    done
    cp /scratch1/zeyunlu/sushie/enrich_header.tsv /scratch1/zeyunlu/sushie_rnaseq/enrich/enrich_rnaseq_${ID}.tsv
    find . -name "enrich.${ID}.*.tsv" | xargs -n 1 tail -n +2 >> /scratch1/zeyunlu/sushie_rnaseq/enrich/enrich_rnaseq_${ID}.tsv
  else
    echo "No CS files."
  fi
done
