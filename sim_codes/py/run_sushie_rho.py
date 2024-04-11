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
    ans = []
    for rhos in [0, 0.1, 0.5, 0.8]:
        sushie = infer_sushie(Xs=X, ys=y, L=args.L2, rho=[rhos])

        covar1 = jnp.transpose(sushie.posteriors.weighted_sum_covar)[0, 0]
        covar2 = jnp.transpose(sushie.posteriors.weighted_sum_covar)[1, 1]
        est_covar = jnp.transpose(sushie.posteriors.weighted_sum_covar)[0, 1]
        est_rho = est_covar / jnp.sqrt(covar1 * covar2)
        ans.append(est_rho)

    df_rho = pd.DataFrame(data={"est_rho0": ans[0],
                                "est_rho1": ans[1],
                                "est_rho5": ans[2],
                                "est_rho8": ans[3]})
    df_rho["Lidx"] = df_rho.index + 1
    df_rho = add_annotation(df_rho, args)
    df_rho.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.rho.tsv", sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
