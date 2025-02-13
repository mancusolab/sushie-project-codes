#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-2200
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

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

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
  if find "/scratch1/zeyunlu/sushie_rnaseq/sushie/alphas/" -maxdepth 1 -name "${ID}.*.alphas.tsv" -print -quit | grep -q .;
  then
    TMPDIR=/scratch1/zeyunlu/trash_rnaseq/${ID}
    mkdir -p $TMPDIR
    cd $TMPDIR

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/enrich_prepare.py \
      --weights-path /scratch1/zeyunlu/sushie_rnaseq/ \
      --ans 3 \
      --gene $ID \
      --out $TMPDIR/enrich.${ID}.prepare.tsv

    if [ ! -e "$TMPDIR/enrich.${ID}.prepare.tsv" ]; then
      continue
    fi

    for JDX in `seq 77`
    do
      params=`sed "${JDX}q;d" /project/nmancuso_8/data/sushie/analysis_results/new_annot/mesa_anno_list.tsv`
      set -- junk $params
      shift
      PRE=$1
      SUF=$2

      echo SUMMARY ${IDX}, ${JDX}, ${ID}, ${CHR}, $PRE, $SUF

      python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/enrich_calc.py \
      --weights-path /scratch1/zeyunlu/sushie_rnaseq/ \
      --weights-file $TMPDIR/enrich.${ID}.prepare.tsv \
      --anno-file /project/nmancuso_8/data/sushie/analysis_results/new_annot/processed/mesa/${PRE}.chr${CHR}.${SUF} \
      --gene $ID \
      --exclude /project/nmancuso_8/data/sushie/analysis_results/new_annot/anno_exclude.tsv \
      --enrich-out $TMPDIR/enrich.${ID}.${PRE}.${SUF}.final.tsv

    done

    cp ${SCRATCH}/enrich_header.tsv /scratch1/zeyunlu/sushie_rnaseq/enrich/enrich_rnaseq_${ID}.tsv
    find . -name "enrich.${ID}.*.final.tsv" | xargs -n 1 tail -n +2 >> /scratch1/zeyunlu/sushie_rnaseq/enrich/enrich_rnaseq_${ID}.tsv
  else
    echo "No CS files."
  fi
  rm -rf $TMPDIR

done
