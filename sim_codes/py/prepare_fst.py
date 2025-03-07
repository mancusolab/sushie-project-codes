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

    rng_key, _, _, _, _, _, g, _ = hf.simulation_sushie(rng_key, output_dic, N, args.L1, args.L3, h2g, rho)
    snps[["snp"]].iloc[g].to_csv(f"{args.tmp_output}.name", sep="\t", index=False, header=None)

    l2_diff = np.abs(np.linalg.norm(output_dic["trueLD"][0][:,g][g,:]) - np.linalg.norm(output_dic["trueLD"][1][:,g][g,:]))
    diff_l2 = np.linalg.norm(output_dic["trueLD"][0][:,g][g,:] - output_dic["trueLD"][1][:,g][g,:])
    
    df=pd.DataFrame({"sim": args.sim, "locus": args.locus, "l2_diff": l2_diff, "diff_l2": diff_l2}, index=[0])
    df.to_csv(f"{args.tmp_output}.causal.tsv", sep="\t", index=False, header=False)
    
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
