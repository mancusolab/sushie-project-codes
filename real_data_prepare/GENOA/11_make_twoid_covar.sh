

cd /project/nmancuso_8/data/sushie/meta_data

awk 'BEGIN {FS=OFS="\t"} {print $1, $0}' ea_covars_sushie.tsv > ea_covars_sushie_2iid_col.tsv
awk 'BEGIN {FS=OFS="\t"} {print $1, $0}' aa_covars_sushie.tsv > aa_covars_sushie_2iid_col.tsv

