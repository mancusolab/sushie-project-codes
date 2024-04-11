#!/usr/bin/env python3

import pandas as pd
dbgap_mapper = pd.read_csv("/project/nmancuso_8/data/MESA/dbgap/WGS/freeze9/sample/sample-info.csv", skiprows=10,
                           header=None).rename(columns={0: "dbgap_id", 1: "geno_id"})
dbgap_mapper = dbgap_mapper[["dbgap_id", "geno_id"]]

ans_info = pd.read_csv("/project/nmancuso_8/data/MESA/processed/covariates/MESA.covar.tsv.gz", sep="\t")

ans_info = ans_info[["sidno", "race"]]
ans_info.columns = ["dbgap_id", "race"]

df_pt = dbgap_mapper.merge(ans_info, how="left", on="dbgap_id")

for ans in ["EUR", "AFR", "HIS", "EAS"]:
    tmp_pt = df_pt[df_pt.race == ans]["geno_id"]
    tmp_pt.to_csv(f"/project/nmancuso_8/data/MESA/processed/geno/plink/ld_pruned/{ans}_genoid.tsv", sep="\t",
                  header=False, index=False)

