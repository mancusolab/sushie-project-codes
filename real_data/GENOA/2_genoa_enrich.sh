#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-3510
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
  SCRATCHDIR=`echo $PWD`
else
  NR=$SLURM_ARRAY_TASK_ID
fi

# 13905
# 13905 = 5 * 2781

source /home1/zeyunlu/init.sh
conda activate jax2

SCRATCH=/scratch1/zeyunlu/sushie

start=`python -c "print( 1 + 5 *int(int($NR-1)))"`
stop=$((start + 4))

for IDX in `seq $start $stop`
do
  # Extract parameters for sushie run
  params=`sed "${IDX}q;d"  ${SCRATCH}/genoa_sushie_gene_list_noMHC.tsv`
  set -- junk $params
  shift
  CHR=$4
  ID=$2
  echo STARTING ${NR}, ${IDX}, ${ID}

  if find "/scratch1/zeyunlu/sushie_genoa/cs/" -maxdepth 1 -name "${ID}.*.cs.tsv" -print -quit | grep -q .;
  then
    TMPDIR=/scratch1/zeyunlu/trash_genoa/${ID}
    mkdir $TMPDIR
    cd $TMPDIR
    for JDX in `seq 77`
    do
      params=`sed "${JDX}q;d"  /scratch1/zeyunlu/new_annot/genoa_anno_list.tsv`
      set -- junk $params
      shift
      PRE=$1
      SUF=$2
      echo SUMMARY ${NR}, ${IDX}, ${JDX}, ${ID}, ${CHR}, ${PRE}, ${SUF}

      python /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/real_data/utils/calc_enrich_all.py \
      --alpha-path /scratch1/zeyunlu/sushie_genoa/alphas/ \
      --anno-file /scratch1/zeyunlu/new_annot/processed/genoa/${PRE}.chr${CHR}.${SUF} \
      --gene $ID \
      --exclude /scratch1/zeyunlu/new_annot/anno_exclude.tsv \
      --enrich-out $TMPDIR/enrich.${ID}.${PRE}.${SUF}.tsv
    done

    cp /scratch1/zeyunlu/sushie/enrich_header.tsv /scratch1/zeyunlu/sushie_genoa/enrich/enrich_genoa_${ID}.tsv
    find . -name "enrich.${ID}.*.tsv" | xargs -n 1 tail -n +2 >> /scratch1/zeyunlu/sushie_genoa/enrich/enrich_genoa_${ID}.tsv
  else
    echo "No CS files."
  fi

done
