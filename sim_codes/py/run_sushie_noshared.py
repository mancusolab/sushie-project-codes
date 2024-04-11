import argparse as ap
import sys
import warnings

import pandas as pd
from jax import random
from jax.config import config
from sushie.infer import infer_sushie

import helper_function as hf

config.update("jax_enable_x64", True)

# annoying warnings if on Mac ARM M1
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp


def add_annotation(dd, args):
    dd = dd.assign(
        sim=args.sim,
        locus=args.locus,
        N=args.N,
        L1=args.L1,
        L2=args.L2,
        L3=args.L3,
        h2g=args.h2g,
        rho=args.rho,
    )
    return dd


def main(args):
    argp = ap.ArgumentParser(description="")
    argp.add_argument("prefix_pop", help="Prefix to PLINK-formatted data for population")
    argp.add_argument("--N", type=str, help="sample size for phenotype")
    argp.add_argument("--L1", default=2, type=int, help="number of generative L")
    argp.add_argument("--L2", default=2, type=int, help="number of inferential L")
    argp.add_argument("--L3", default=0, type=int, help="number of additional population-specific L")
    argp.add_argument("--h2g", type=str, help="signal to noise ratio for phenotype")
    argp.add_argument("--rho", type=str, default="none", help="cov parameter")
    argp.add_argument("--seed", default=1234, type=int, help="Random seed")
    argp.add_argument("--sim", default=1, help="Simulation Index")
    argp.add_argument("--locus", default=1, help="Locus Index")
    argp.add_argument("--threshold", default=0.9, help="Credible set threshold")
    argp.add_argument("-o", "--output", default=sys.stdout)

    args = argp.parse_args(args)
    rng_key = random.PRNGKey(int(args.seed))

    plink_path = args.prefix_pop.split(":")
    N = args.N.split(":")
    h2g = args.h2g.split(":")

    if args.rho != "none":
        rho = args.rho.split(":")
        rho = [float(i) for i in rho]
    else:
        rho = 0

    N = [int(i) for i in N]
    h2g = [float(i) for i in h2g]

    output_dic = hf._gen_ld(plink_path)

    rng_key, X, y, bvec, bvec_all, s2e, g, b_covar = hf.simulation_sushie(rng_key, output_dic, N, args.L1, args.L3, h2g,
                                                                          rho)

    # perform sushie
    sushie = infer_sushie(Xs=X, ys=y, L=args.L2)
    indep = infer_sushie(Xs=X, ys=y, L=args.L2, rho=[0], no_update=True)

    # meta-sushie
    meta1 = infer_sushie(Xs=[X[0]], ys=[y[0]], L=args.L2)
    meta2 = infer_sushie(Xs=[X[1]], ys=[y[1]], L=args.L2)

    # susie
    susie = infer_sushie([jnp.concatenate((X[0], X[1]))], [jnp.concatenate((y[0], y[1]))], L=args.L2)

    # PIP
    df_pip = pd.DataFrame(data=jnp.array([[1] * args.L3, [2] * args.L3]).flatten(), columns=["ancestry"])
    df_pip["idx"] = g
    df_pip["sushie"] = sushie.pip_all[g]
    df_pip["indep"] = indep.pip_all[g]
    df_pip["meta1"] = meta1.pip_all[g]
    df_pip["meta2"] = meta2.pip_all[g]
    df_pip["meta"] = 1 - (1 - meta1.pip_all[g]) * (1 - meta2.pip_all[g])
    df_pip["susie"] = susie.pip_all[g]

    df_pip = add_annotation(df_pip, args)
    df_pip.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.pip.tsv", sep="\t", index=False)

    meta_cs = pd.concat([meta1.cs, meta2.cs], axis=0).drop_duplicates(subset="SNPIndex")

    prec1 = [
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 1].idx.values, sushie.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 1].idx.values, indep.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 1].idx.values, meta1.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 1].idx.values, meta2.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 1].idx.values, meta_cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 1].idx.values, susie.cs.SNPIndex.values.astype(int))),
    ]

    prec2 = [
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 2].idx.values, sushie.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 2].idx.values, indep.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 2].idx.values, meta1.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 2].idx.values, meta2.cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 2].idx.values, meta_cs.SNPIndex.values.astype(int))),
        jnp.sum(jnp.isin(df_pip[df_pip.ancestry == 2].idx.values, susie.cs.SNPIndex.values.astype(int))),
    ]

    df_prop = pd.DataFrame(data={
        "method": ["sushie", "indep", "meta1", "meta2", "meta", "susie"],
        "prec_ancestry1": prec1,
        "prec_ancestry2": prec2,
    })

    df_prop = add_annotation(df_prop, args)
    df_prop.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.prop.tsv", sep="\t", index=False)

    # sushie_cs = sushie.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    # indep_cs = indep.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    # meta1_cs = meta1.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    # meta2_cs = meta2.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    # meta_cs = meta_cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    # susie_cs = susie.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    #
    # df_cs = pd.DataFrame(data=jnp.arange(args.L2) + 1, columns=["CSIndex"])
    # df_cs = df_cs.merge(sushie_cs, how="left", on="CSIndex").rename(columns={"SNPIndex": "sushie"}) \
    #     .merge(indep_cs, how="left", on="CSIndex").rename(columns={"SNPIndex": "indep"}) \
    #     .merge(meta1_cs, how="left", on="CSIndex").rename(columns={"SNPIndex": "meta1"}) \
    #     .merge(meta2_cs, how="left", on="CSIndex").rename(columns={"SNPIndex": "meta2"}) \
    #     .merge(meta_cs, how="left", on="CSIndex").rename(columns={"SNPIndex": "meta"}) \
    #     .merge(susie_cs, how="left", on="CSIndex").rename(columns={"SNPIndex": "susie"}).fillna(0)
    #
    # df_cs = add_annotation(df_cs, args)
    # df_cs.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.cs.tsv", sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
