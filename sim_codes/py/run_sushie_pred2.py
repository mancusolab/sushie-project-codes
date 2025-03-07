import argparse as ap
import sys
import warnings

import numpy as np
import pandas as pd
from jax.config import config
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
    argp.add_argument("--ngwas", type=int, help="cov parameter")
    argp.add_argument("--h2ge", type=float, help="cov parameter")
    argp.add_argument("--seed", default=1234, type=int, help="Random seed")
    argp.add_argument("--sim", default=1, help="Simulation Index")
    argp.add_argument("--locus", default=1, help="Locus Index")
    argp.add_argument("--threshold", default=0.9, help="Credible set threshold")
    argp.add_argument("--sushie_file")
    argp.add_argument("--mesusie_in_file")
    argp.add_argument("--xmap_in_file")
    argp.add_argument("--xmap_ind_file")
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

    rng_key, total_X, total_y, _, bvec_all, _, _, _ = hf.simulation_sushie(rng_key, output_dic, N, args.L1, args.L3,
                                                                           h2g,
                                                                           rho)

    n_pop = len(total_X)

    valid_X = [total_X[idx][old_N[idx]:(old_N[idx] + 200)] for idx in range(n_pop)]
    valid_y = [total_y[idx][old_N[idx]:(old_N[idx] + 200)] for idx in range(n_pop)]

    pop1_weight = []
    pop2_weight = []
    method_name = []

    # Read the TSV file into a DataFrame
    df_sushie = pd.read_csv(args.sushie_file, sep="\t")
    pop1_weight.append(jnp.array(df_sushie)[:, 0:8])
    pop2_weight.append(jnp.array(df_sushie)[:, 8:16])
    method_name.append(["sushie", "indep", "meta", "susie", "sushie_ss", "enet", "lasso", "ridge"])

    df_mesusie_in = pd.read_csv(args.mesusie_in_file, sep="\t", header=None)
    pop1_weight.append(jnp.array(df_mesusie_in)[:, 0][:, jnp.newaxis])
    pop2_weight.append(jnp.array(df_mesusie_in)[:, 1][:, jnp.newaxis])
    method_name.append(["mesusie.in"])

    df_xmap_in = pd.read_csv(args.xmap_in_file, sep="\t", header=None)
    pop1_weight.append(jnp.array(df_xmap_in)[:, 0][:, jnp.newaxis])
    pop2_weight.append(jnp.array(df_xmap_in)[:, 1][:, jnp.newaxis])
    method_name.append(["xmap.in"])

    df_xmap_ind = pd.read_csv(args.xmap_ind_file, sep="\t", header=None)
    pop1_weight.append(jnp.array(df_xmap_ind)[:, 0][:, jnp.newaxis])
    pop2_weight.append(jnp.array(df_xmap_ind)[:, 1][:, jnp.newaxis])
    method_name.append(["xmap.ind"])

    method_list = [item for sublist in method_name for item in sublist]
    pop1_weight = np.concatenate(pop1_weight, axis=1)
    pop2_weight = np.concatenate(pop2_weight, axis=1)
    est_y1 = jnp.einsum("ij,jk->ik", valid_X[0], pop1_weight)
    est_y2 = jnp.einsum("ij,jk->ik", valid_X[1], pop2_weight)

    rsq_1, _ = hf.ols(valid_y[0][:, np.newaxis], est_y1)
    rsq_2, _ = hf.ols(valid_y[1][:, np.newaxis], est_y2)

    df_r2 = pd.DataFrame(
        data=np.concatenate(
            [rsq_1[:, np.newaxis], rsq_2[:, np.newaxis]],
            axis=1),
        index=method_list,
        columns=["ancestry1_weight1", "ancestry2_weight2"]).reset_index(
        names="method")

    df_r2 = add_annotation(df_r2, args)

    df_r2.to_csv(f"{args.output}.r2.tsv", sep="\t", index=False)

    # for pop1 twas
    rng_key, gwas1_key, gwas2_key = random.split(rng_key, 3)
    gwas1, _ = hf.sim_gwas(output_dic["L"][0], args.ngwas, bvec_all[0], args.h2ge, gwas1_key)
    gwas2, _ = hf.sim_gwas(output_dic["L"][1], args.ngwas, bvec_all[1], args.h2ge, gwas2_key)

    z_1 = hf.compute_twas(gwas1, pop1_weight, output_dic["LD"][0])
    z_2 = hf.compute_twas(gwas2, pop2_weight, output_dic["LD"][1])

    df_z = pd.DataFrame(data=np.append(z_1[np.newaxis, :], z_2[np.newaxis, :], axis=0),
                        columns=method_list,
                        index=[1, 2]).reset_index(names="ancestry")

    df_z = add_annotation(df_z, args)

    df_z.to_csv(f"{args.output}.twas.tsv", sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
