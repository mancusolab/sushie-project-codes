import argparse as ap
import sys
import warnings

import numpy as np
import pandas as pd
from jax.config import config
from sushie.infer import infer_sushie
from sushie.infer_ss import infer_sushie_ss
from jax import random
import helper_function as hf

config.update("jax_enable_x64", True)

# annoying warnings if on Mac ARM M1
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp


def add_annotation(dd, args):
    dd = dd.assign(sim=args.sim,
                   locus=args.locus,
                   N=args.N,
                   L1=args.L1,
                   L2=args.L2,
                   L3=args.L3,
                   h2g=args.h2g,
                   rho=args.rho,
                   h2ge=args.h2ge,
                   ngwas=args.ngwas,
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

    old_N = [int(i) for i in N]
    N = [int(i) + 200 for i in N]
    h2g = [float(i) for i in h2g]

    output_dic, snps = hf._gen_ld(plink_path)
    snps.columns = ["chrom", "snp", "pos_1", "a0", "a1", "bimIDX_0", "pos_1bp", "bimIDX_1", "pos_2", "bimIDX_2",
                    "pos_2bp", "bimIDX_3"]
    snps = snps[["chrom", "snp", "pos_1", "a0", "a1", "bimIDX_0", "bimIDX_1", "pos_2", "bimIDX_2", "bimIDX_3"]]

    rng_key, total_X, total_y, _, bvec_all, _, g, _ = hf.simulation_sushie(rng_key, output_dic, N, args.L1, args.L3,
                                                                           h2g,
                                                                           rho)
    n_pop = len(total_X)
    X = [total_X[idx][0:old_N[idx]] for idx in range(n_pop)]
    y = [total_y[idx][0:old_N[idx]] for idx in range(n_pop)]

    # output X, y, summary stats, and ld for susiex
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
                   header=None, float_format="%.12f")
        pd.DataFrame(jnp.concatenate([y[idx][:, jnp.newaxis], X[idx]], axis=1)). \
            to_csv(f"{args.tmp_output}.xy.sim{args.sim}.locus{args.locus}.ans{idx}.tsv", index=False, sep="\t",
                   header=None, float_format="%.12f")

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
    pd.DataFrame({"SNP": snps.iloc[g, :].snp.values, "SNPIndex_0based": g, "SNPIndex_1based": g+1}).\
        to_csv(f"{args.tmp_output}.causal.sim{args.sim}.locus{args.locus}.tsv", index=False, sep="\t")

    # re-train the weights
    sushie = infer_sushie(X, y, L=args.L2)

    indep = infer_sushie(X, y, L=args.L2, rho=[0], no_update=True)

    # meta-sushie
    meta1 = infer_sushie([X[0]], [y[0]], L=args.L2)
    meta2 = infer_sushie([X[1]], [y[1]], L=args.L2)

    # susie
    rb_X = jnp.concatenate((X[0], X[1]))
    rb_y = jnp.concatenate((y[0], y[1]))

    susie = infer_sushie([rb_X], [rb_y], L=args.L2)

    sushie_weight = jnp.sum(sushie.posteriors.post_mean, axis=0)
    indep_weight = jnp.sum(indep.posteriors.post_mean, axis=0)
    meta1_weight = jnp.sum(meta1.posteriors.post_mean, axis=0)
    meta2_weight = jnp.sum(meta2.posteriors.post_mean, axis=0)
    susie_weight = jnp.sum(susie.posteriors.post_mean, axis=0)

    _, h2g, _, _, _ = hf.estimate_her(rb_X, rb_y)

    enet1, _, _ = hf.fit_enet(args.seed, rb_X, rb_y, h2g)

    lasso1, _, _ = hf.fit_lasso(args.seed, rb_X, rb_y, h2g)

    ridge1, _, _ = hf.fit_ridge(args.seed, rb_X, rb_y, h2g)

    z_ss = jnp.vstack([np.array(ss_list[0].z.values), np.array(ss_list[1].z.values)])
    lds = jnp.stack(output_dic["trueLD"][0:2], axis=0)
    sushie_ss = infer_sushie_ss(zs=z_ss, lds=lds, ns=jnp.array(N)[:, jnp.newaxis], L=args.L2, threshold=0.95)
    sushie_ss_weight = jnp.sum(sushie_ss.posteriors.post_mean, axis=0)

    pop1_weight = np.concatenate([sushie_weight[:, 0][:, np.newaxis],
                                  indep_weight[:, 0][:, np.newaxis],
                                  meta1_weight,
                                  susie_weight,
                                  sushie_ss_weight[:, 0][:, np.newaxis],
                                  enet1[:, np.newaxis],
                                  lasso1[:, np.newaxis],
                                  ridge1[:, np.newaxis],
                                  ], axis=1)

    pop2_weight = np.concatenate([sushie_weight[:, 1][:, np.newaxis],
                                  indep_weight[:, 1][:, np.newaxis],
                                  meta2_weight,
                                  susie_weight,
                                  sushie_ss_weight[:, 1][:, np.newaxis],
                                  enet1[:, np.newaxis],
                                  lasso1[:, np.newaxis],
                                  ridge1[:, np.newaxis],
                                  ], axis=1)

    pd.DataFrame(np.concatenate([pop1_weight, pop2_weight], axis=1)).\
        to_csv(f"{args.tmp_output}.weights.sim{args.sim}.locus{args.locus}.tsv", index=False, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
