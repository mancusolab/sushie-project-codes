import argparse as ap
import sys
import warnings

import pandas as pd
from jax import random
from jax.config import config
from sushie.infer import infer_sushie
from sushie.infer_ss import infer_sushie_ss
import helper_function as hf
import numpy as np

import MultiSuSiE

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
    argp.add_argument("--tmp_output")
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

    output_dic, snps = hf._gen_ld(plink_path)
    snps.columns = ["chrom", "snp", "pos_1", "a0", "a1", "bimIDX_0", "pos_1bp", "bimIDX_1", "pos_2", "bimIDX_2", "pos_2bp", "bimIDX_3"]
    snps = snps[["chrom", "snp", "pos_1", "a0", "a1", "bimIDX_0", "bimIDX_1", "pos_2", "bimIDX_2", "bimIDX_3"]]

    rng_key, X, y, bvec, bvec_all, s2e, g, b_covar = hf.simulation_sushie(rng_key, output_dic, N, args.L1, args.L3, h2g,
                                                                          rho)
    # output X, y, summary stats, and ld for susiex, mesusie, and xmap
    ss_list = []
    for idx in range(len(X)):
        df_snps = hf.regress(X[idx], y[idx])
        df_snps2 = df_snps.copy()
        df_snps2["z"] = df_snps2.beta / df_snps2.se
        ss_list.append(df_snps2)
        snps_insample = snps[["chrom", "snp", "pos_1", "a0", "a1"]]
        snps_insample.columns = ["chrom", "snp", "pos", "a0", "a1"]
        new_snps1 = pd.concat([snps_insample, df_snps], axis=1)
        new_snps1.to_csv(f"{args.tmp_output}.inss.sim{args.sim}.locus{args.locus}.ans{idx}.tsv", index=False, sep="\t")

        snps_outsample = snps[["chrom", "snp", "pos_2", "a0", "a1"]]
        snps_outsample.columns = ["chrom", "snp", "pos", "a0", "a1"]
        new_snps2 = pd.concat([snps_outsample, df_snps], axis=1)
        new_snps2.to_csv(f"{args.tmp_output}.outss.sim{args.sim}.locus{args.locus}.ans{idx}.tsv", index=False, sep="\t")

        pd.DataFrame(output_dic["trueLD"][idx]).\
            to_csv(f"{args.tmp_output}.inld.sim{args.sim}.locus{args.locus}.ans{idx}.tsv", index=False, sep="\t",
                   header=None, float_format="%.24f")
        pd.DataFrame(jnp.concatenate([y[idx][:, jnp.newaxis], X[idx]], axis=1)). \
            to_csv(f"{args.tmp_output}.xy.sim{args.sim}.locus{args.locus}.ans{idx}.tsv", index=False, sep="\t",
                   header=None, float_format="%.24f")

    # for idx in range(len(X)):
    #     pd.DataFrame(output_dic["trueLD"][idx + 2]). \
    #         to_csv(f"{args.tmp_output}.outld.sim{args.sim}.locus{args.locus}.ans{idx}.tsv", index=False, sep="\t",
    #                header=None, float_format="%.12f")

    # compute ldscore for xmap
    ldsc1 = np.sum(output_dic["trueLD"][0] ** 2, axis=0)
    ldsc2 = np.sum(output_dic["trueLD"][1] ** 2, axis=0)
    ldsc3 = np.sum(output_dic["trueLD"][0] * output_dic["trueLD"][1], axis=0)
    pd.DataFrame(jnp.concatenate([ldsc1[:, jnp.newaxis], ldsc2[:, jnp.newaxis], ldsc3[:, jnp.newaxis]], axis=1)).\
        to_csv(f"{args.tmp_output}.inldsc.sim{args.sim}.locus{args.locus}.tsv", index=False, sep="\t", header=None)

    # ldsc1 = np.sum(output_dic["trueLD"][2] ** 2, axis=0)
    # ldsc2 = np.sum(output_dic["trueLD"][3] ** 2, axis=0)
    # ldsc3 = np.sum(output_dic["trueLD"][2] * output_dic["trueLD"][3], axis=0)
    # pd.DataFrame(jnp.concatenate([ldsc1[:, jnp.newaxis], ldsc2[:, jnp.newaxis], ldsc3[:, jnp.newaxis]], axis=1)). \
    #     to_csv(f"{args.tmp_output}.outldsc.sim{args.sim}.locus{args.locus}.tsv", index=False, sep="\t", header=None)

    # in susiex, it uses snp name for the results
    # in R, index starts with 1
    pd.DataFrame({"SNP": snps.iloc[g[0:args.L1], :].snp.values, "SNPIndex_0based": g[0:args.L1],
                  "SNPIndex_1based": g[0:args.L1]+1,
                  "CSIndex": jnp.argsort(-jnp.array(bvec**2)[:, 0])[0:args.L1] + 1}).\
        to_csv(f"{args.tmp_output}.causal.sim{args.sim}.locus{args.locus}.tsv", index=False, sep="\t")

    # tracks some meta data
    df_met = pd.DataFrame({"N": N,
                           "BP": [jnp.min(new_snps1.pos.values), jnp.max(new_snps1.pos.values)],
                           "BP2": [jnp.min(new_snps2.pos.values), jnp.max(new_snps2.pos.values)]})
    df_met.to_csv(f"{args.tmp_output}.md.sim{args.sim}.locus{args.locus}.tsv", index=False, sep="\t", header=None)

    # perform multisusie
    multisusie1_converged = True
    try:
        multisusie1 = MultiSuSiE.multisusie_rss(z_list=[np.array(ss_list[0].z.values), np.array(ss_list[1].z.values)],
                                                        R_list=output_dic["trueLD"][0:2], rho=np.array([[1, 0.1], [0.1, 1]]),
                                                        population_sizes=N, L=args.L2, max_iter=500, tol=0.0001,
                                                        min_abs_corr=0.5)
    except Exception as e:
        multisusie1_converged = False

    multisusie2_converged = True
    try:
        multisusie2 = MultiSuSiE.multisusie(X_list=[np.array(X[0]), np.array(X[1])], Y_list=[np.array(y[0]), np.array(y[1])], rho=np.array([[1, 0.1], [0.1, 1]]), L=args.L2,standardize=True,intercept=False,float_type=np.float64)
    except Exception as e:
        multisusie2_converged = False


    multisusie1_cs = []
    multisusie2_cs = []
    for idx in range(args.L2):
        if multisusie1_converged and multisusie1.sets[3][idx]:
            multisusie1_cs.append(pd.DataFrame({"CSIndex": idx + 1, "SNPIndex": multisusie1.sets[0][idx]}))
        else:
            multisusie1_cs.append(pd.DataFrame({"CSIndex": [], "SNPIndex": []}))

        if multisusie2_converged and multisusie2.sets[3][idx]:
            multisusie2_cs.append(pd.DataFrame({"CSIndex": idx+1, "SNPIndex": multisusie2.sets[0][idx]}))
        else:
            multisusie2_cs.append(pd.DataFrame({"CSIndex": [], "SNPIndex": []}))

    multisusie1_cs = pd.concat(multisusie1_cs, axis=0)
    multisusie2_cs = pd.concat(multisusie2_cs, axis=0)

    # perform sushie
    sushie = infer_sushie(Xs=X, ys=y, L=args.L2, threshold=0.95)
    indep = infer_sushie(Xs=X, ys=y, L=args.L2, rho=[0], no_update=True, threshold=0.95)
    meta1 = infer_sushie(Xs=[X[0]], ys=[y[0]], L=args.L2, threshold=0.95)
    meta2 = infer_sushie(Xs=[X[1]], ys=[y[1]], L=args.L2, threshold=0.95)
    susie = infer_sushie([jnp.concatenate((X[0], X[1]))], [jnp.concatenate((y[0], y[1]))], L=args.L2, threshold=0.95)

    # sushie summary stats
    z_ss = jnp.vstack([np.array(ss_list[0].z.values), np.array(ss_list[1].z.values)])
    lds = jnp.stack(output_dic["trueLD"][0:2], axis=0)
    sushie_ss = infer_sushie_ss(zs=z_ss, lds=lds, ns=jnp.array(N)[:, jnp.newaxis], L=args.L2, threshold=0.95)

    df_sens = pd.DataFrame({"sushie": [1 * sushie.elbo_increase],
                            "indep": [1 * indep.elbo_increase],
                            "meta": [1 * (meta1.elbo_increase and meta2.elbo_increase)],
                            "susie": [1 * susie.elbo_increase],
                            "susie_ss": [1 * sushie_ss.elbo_increase],
                            "multisusie.ss": [1 * multisusie1_converged],
                            "multisusie.ind": [1 * multisusie2_converged]})

    df_sens = add_annotation(df_sens, args)
    df_sens.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.sens.tsv", sep="\t", index=False)

    # track est_rho
    covar1 = jnp.transpose(sushie.posteriors.weighted_sum_covar)[0, 0]
    covar2 = jnp.transpose(sushie.posteriors.weighted_sum_covar)[1, 1]
    est_covar = jnp.transpose(sushie.posteriors.weighted_sum_covar)[0, 1]
    est_rho = est_covar / jnp.sqrt(covar1 * covar2)
    df_rho = pd.DataFrame(data=est_rho, columns=["est_rho"])
    df_rho["Lidx"] = df_rho.index + 1
    df_rho["method"] = "sushie"
    covar1 = jnp.transpose(sushie_ss.posteriors.weighted_sum_covar)[0, 0]
    covar2 = jnp.transpose(sushie_ss.posteriors.weighted_sum_covar)[1, 1]
    est_covar = jnp.transpose(sushie_ss.posteriors.weighted_sum_covar)[0, 1]
    est_rho = est_covar / jnp.sqrt(covar1 * covar2)
    df_rho2 = pd.DataFrame(data=est_rho, columns=["est_rho"])
    df_rho2["Lidx"] = df_rho.index + 1
    df_rho2["method"] = "sushie_ss"
    df_rho = pd.concat([df_rho, df_rho2], axis=0)
    df_rho = add_annotation(df_rho, args)
    df_rho.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.rho.tsv", sep="\t", index=False)

    meta_cs = pd.concat([meta1.cs, meta2.cs], axis=0).drop_duplicates(subset="SNPIndex")

    # PIP
    df_pip = pd.DataFrame()
    df_pip["CSIndex"] = jnp.argsort(-jnp.array(bvec**2)[:, 0]) + 1
    df_pip["SNPIndex_0based"] = g[0:args.L1]
    df_pip["SNPIndex_1based"] = g[0:args.L1] + 1
    df_pip["sushie_pip"] = sushie.pip_all[g[0:args.L1]]
    df_pip["sushie_cali"] = jnp.isin(g[0:args.L1], sushie.cs.SNPIndex.values.astype(int)).astype(int)
    df_pip["indep_pip"] = indep.pip_all[g[0:args.L1]]
    df_pip["indep_cali"] = jnp.isin(g[0:args.L1], indep.cs.SNPIndex.values.astype(int)).astype(int)
    df_pip["meta_pip"] = 1 - (1 - meta1.pip_all[g[0:args.L1]]) * (1 - meta2.pip_all[g[0:args.L1]])
    df_pip["meta_cali"] = jnp.isin(g[0:args.L1], meta_cs.SNPIndex.values.astype(int)).astype(int)
    df_pip["susie_pip"] = susie.pip_all[g[0:args.L1]]
    df_pip["susie_cali"] = jnp.isin(g[0:args.L1], susie.cs.SNPIndex.values.astype(int)).astype(int)
    df_pip["susie_ss_pip"] = sushie_ss.pip_all[g[0:args.L1]]
    df_pip["susie_ss_cali"] = jnp.isin(g[0:args.L1], sushie_ss.cs.SNPIndex.values.astype(int)).astype(int)

    if multisusie1_converged:
        df_pip["multisusie.ss_pip"] = hf.make_pip(multisusie1.alpha)[g[0:args.L1]]
        df_pip["multisusie.ss_cali"] = jnp.isin(g[0:args.L1], multisusie1_cs.SNPIndex.values.astype(int)).astype(int)
    else:
        df_pip["multisusie.ss_pip"] = jnp.nan
        df_pip["multisusie.ss_cali"] = jnp.nan

    if multisusie2_converged:
        df_pip["multisusie.ind_pip"] = hf.make_pip(multisusie2.alpha)[g[0:args.L1]]
        df_pip["multisusie.ind_cali"] = jnp.isin(g[0:args.L1], multisusie2_cs.SNPIndex.values.astype(int)).astype(int)
    else:
        df_pip["multisusie.ind_pip"] = jnp.nan
        df_pip["multisusie.ind_cali"] = jnp.nan

    df_pip = add_annotation(df_pip, args)
    df_pip.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.pip.tsv", sep="\t", index=False)

    # credible set
    df_cs = pd.DataFrame(data=jnp.arange(args.L2) + 1, columns=["CSIndex"])
    sushie_cs = sushie.cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "sushie"})
    indep_cs = indep.cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "indep"})
    meta_cs = meta_cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "meta"})
    susie_cs = susie.cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "susie"})
    sushie_ss_cs = sushie_ss.cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "susie_ss"})
    multisusie1_cs = multisusie1_cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "multisusie.ss"})
    multisusie2_cs = multisusie2_cs.groupby("CSIndex")["SNPIndex"].count().reset_index().rename(columns={"SNPIndex": "multisusie.ind"})

    df_cs = df_cs.merge(sushie_cs, how="left", on="CSIndex").merge(indep_cs, how="left", on="CSIndex").\
        merge(meta_cs, how="left", on="CSIndex") \
        .merge(susie_cs, how="left", on="CSIndex")\
        .merge(sushie_ss_cs, how="left", on="CSIndex")\
        .merge(multisusie1_cs, how="left", on="CSIndex")\
        .merge(multisusie2_cs, how="left", on="CSIndex")

    df_cs = add_annotation(df_cs, args)
    df_cs.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.cs.tsv", sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
