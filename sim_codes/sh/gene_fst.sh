
DATAF=/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/

PLINK=/project/nmancuso_8/zeyunlu/tools/plink2
PLINK1=/project/nmancuso_8/zeyunlu/tools/plink/plink

for IDX in `seq 1 500`
do
    params=`sed "${IDX}q;d" /project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_data_prepare/sim_gene_index_list.tsv`
    echo "NR=${IDX}: ${params}"
    set -- junk $params
    shift

    row=$1
    chr=$2
    gene=$3
    
    echo $row $chr $gene
    TMPDIR=/scratch1/zeyunlu/tmp_gene_fst/

    eur1=$DATAF/new/EUR_1000G/${gene}_${chr}_geno
    afr1=$DATAF/new/AFR_1000G/${gene}_${chr}_geno

    ${PLINK1} --bfile ${eur1} --bmerge ${afr1} --make-bed --out ${TMPDIR}/${gene}_${chr}_merged

    awk -F'\t' 'OFS="\t" {print $1, $2, "EUR"}' ${eur1}.fam > ${TMPDIR}/${gene}_${chr}_all.pt
    awk -F'\t' 'OFS="\t" {print $1, $2, "AFR"}' ${afr1}.fam >> ${TMPDIR}/${gene}_${chr}_all.pt

    ${PLINK} --bfile ${TMPDIR}/${gene}_${chr}_merged \
    --fst CATPHENO --within  $TMPDIR/${gene}_${chr}_all.pt --out ${TMPDIR}/${gene}_${chr}_fst

    awk -F'\t' -v r="$row" -v l="$gene" -v c="$chr" 'NR==2 {print $3, r, c, l}' OFS='\t' ${TMPDIR}/${gene}_${chr}_fst.fst.summary > ${TMPDIR}/${gene}_${chr}_gene_fst.tsv
done
