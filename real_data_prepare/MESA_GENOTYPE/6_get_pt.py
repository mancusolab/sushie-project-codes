#!/usr/bin/env python3

import pandas as pd

all = []
for ans in ["EUR", "AFR", "HIS", "EAS"]:
    for type in ["proteins", "rnaseq"]:
        df_tmp = pd.read_csv(f"/scratch1/zeyunlu/sushie/mesa_{type}_pt_{ans}.tsv", sep="\t", header=None)
        all.append(df_tmp)

all = pd.concat(all).drop_duplicates()
all.to_csv("/scratch1/zeyunlu/sushie/mesa_all_pt.tsv", header=None, sep="\t", index=False)
