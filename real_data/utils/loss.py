#! /usr/bin/env python
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--eurfile", type=str)
parser.add_argument("--afrfile", type=str)
parser.add_argument("--hisfile", type=str)
parser.add_argument("--easfile", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

df_eur = pd.read_csv(args.eurfile, sep="\t", header=None)
df_afr = pd.read_csv(args.afrfile, sep="\t", header=None)
df_his = pd.read_csv(args.hisfile, sep="\t", header=None)
df_eas = pd.read_csv(args.easfile, sep="\t", header=None)

pop2_num = df_eur[[1]].merge(df_afr[[1]], how="inner", on=1).shape[0]

pop3_num = df_eur[[1]].merge(df_afr[[1]], how="inner", on=1).merge(df_his[[1]], how="inner", on=1).shape[0]

pop4_num = df_eur[[1]].merge(df_afr[[1]], how="inner", on=1).merge(df_his[[1]], how="inner", on=1).\
    merge(df_eas[[1]], how="inner", on=1).shape[0]

pd.DataFrame(data={"gene": [args.gene], "2pop": [pop2_num], "3pop": [pop3_num], "4pop": [pop4_num],})\
    .to_csv(args.out, sep="\t", index=False)
