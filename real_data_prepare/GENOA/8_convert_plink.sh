#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --partition=conti

NR=$1

# Extract parameters for sushie run
params=`sed "${NR}q;d" /project/nmancuso_8/data/GENOA/sushie/ans_chr_list.tsv`
set -- junk $params
shift
POP=$1
CHR=$2

echo $POP $CHR

source /home1/zeyunlu/init.sh
conda activate bcf

ANO=/project/nmancuso_8/data/dbSNP/dbSNP_v155/split/GCF_000001405.chr${CHR}.renamed.39.gz
CHR_CONVERT=/project/nmancuso_8/data/GENOA/sushie/chr_convert.tsv

VCF=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/filter_06/${POP}_chr${CHR}_weighted_r2_06.vcf.gz
TEMPFILE=/scratch1/zeyunlu/tmpgeno/temp_${CHR}_${POP}.vcf.gz
OUT=/project/nmancuso_8/data/GENOA/processed/genotype/vcf/postimputed/annotated_dbsnp155/${POP}_chr${CHR}_weighted_r2_06_annotated_dbsnp_v155.vcf.gz

echo "index vcf"
bcftools index -f $VCF
echo "annotate chr"
bcftools annotate --rename-chrs ${CHR_CONVERT} $VCF | bgzip > $TEMPFILE
echo "index temp file"
bcftools index -f $TEMPFILE
echo "annotate ID"
bcftools annotate -c ID -a $ANO $TEMPFILE -Oz -o $OUT

path=/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2

$PLINK --vcf $OUT --make-bed --out ${path}/${POP}_chr${CHR} --silent
