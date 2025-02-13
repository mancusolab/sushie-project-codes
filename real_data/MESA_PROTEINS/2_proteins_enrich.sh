#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-637
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu
#SBATCH --partition=conti

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 1274
# 1274 = 2 * 637

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/project/nmancuso_8/data/sushie/meta_data

start=`python -c "print( 1 + 2 *int(int($NR-1)))"`
stop=$((start + 1))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/mesa_proteins_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  CHR=$4
  ID=${15}
  TMPDIR=/scratch1/zeyunlu/trash_proteins/${ID}
  mkdir -p $TMPDIR

  if find "/scratch1/zeyunlu/sushie_proteins/sushie/alphas/" -maxdepth 1 -name "${ID}.*.alphas.tsv" -print -quit | grep -q .;
  then

    python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/enrich_prepare.py \
      --weights-path /scratch1/zeyunlu/sushie_proteins/ \
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
      --weights-path /scratch1/zeyunlu/sushie_proteins/ \
      --weights-file $TMPDIR/enrich.${ID}.prepare.tsv \
      --anno-file /project/nmancuso_8/data/sushie/analysis_results/new_annot/processed/mesa/${PRE}.chr${CHR}.${SUF} \
      --gene $ID \
      --exclude /project/nmancuso_8/data/sushie/analysis_results/new_annot/anno_exclude.tsv \
      --enrich-out $TMPDIR/enrich.${ID}.${PRE}.${SUF}.final.tsv

    done

    cp ${SCRATCH}/enrich_header.tsv /scratch1/zeyunlu/sushie_proteins/enrich/enrich_proteins_${ID}.tsv
    find $TMPDIR -name "enrich.${ID}.*.final.tsv" | xargs -n 1 tail -n +2 >> /scratch1/zeyunlu/sushie_proteins/enrich/enrich_proteins_${ID}.tsv
  else
    echo "No CS files."
  fi

  rm -rf ${TMPDIR}

done
