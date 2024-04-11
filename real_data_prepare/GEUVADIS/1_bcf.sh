#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64Gb
#SBATCH --array=1-22
#SBATCH --mail-type=all
#SBATCH --mail-user=zeyunlu@usc.edu

source /home1/zeyunlu/init.sh

CHR=$SLURM_ARRAY_TASK_ID
ORI=/project/nmancuso_8/data/GEUVADIS/GenotypeData/vcf/GEUVADIS.chr${CHR}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz
TMP1=/scratch/zeyunlu/tmp/geuvadis/vcf/GEUVADIS.chr${CHR}.vcf
TMP2=/scratch/zeyunlu/tmp/geuvadis/vcf/GEUVADIS.chr${CHR}.add.format.vcf
ANO=/project/nmancuso_8/data/dbSNP_v153/GCF_000001405.25.fixed.vcf.gz
OUT=/project/nmancuso_8/data/GEUVADIS/GenotypeData/annotated_vcf_dbsnpv153/GEUVADIS.chr${CHR}.annotated.dbSNP_v153.vcf.gz

echo "Copying"
cp $ORI ${TMP1}.gz
gunzip ${TMP1}.gz
bgzip ${TMP1}
echo "Indexing"
bcftools index ${TMP1}.gz
echo "Edit vcf"
python /home1/zeyunlu/research/sub_mefocus/geuv/edit_vcf.py ${TMP1}.gz ${TMP2}
echo "Zipping vcf"
bgzip ${TMP2}
echo "Indexing"
bcftools index ${TMP2}.gz
# tabix -p vcf ${TMP2}.gz
echo "Annotating"
bcftools annotate -c ID -a $ANO $TMP2.gz -Oz -o $OUT
