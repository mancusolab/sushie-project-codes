import argparse as ap
import sys
import warnings

import pandas as pd
from jax.config import config
from sushie.infer import _infer_test


config.update("jax_enable_x64", True)

# annoying warnings if on Mac ARM M1
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp


def main(args):
    argp = ap.ArgumentParser(description="")
    argp.add_argument("--pheno", type=str, nargs="+")
    argp.add_argument("--geno", type=str, nargs="+")
    argp.add_argument("--trait", type=str)
    argp.add_argument("--output", type=str)

    args = argp.parse_args(args)

    pheno = []
    geno = []
    for idx in range(2):
        tmp_pheno = jnp.array(pd.read_csv(args.pheno[idx], sep="\t", header=None))
        tmp_geno = jnp.array(pd.read_csv(args.geno[idx], sep="\t", header=None))
        pheno.append(tmp_pheno)
        geno.append(tmp_geno)

    import time
    start_time = time.time()
    _infer_test(Xs=geno, ys=pheno, L=10, no_scale=True)
    end_time = time.time()

    time1 = (end_time - start_time)

    df_res = pd.DataFrame({"method": ["sushie_ind"], "time": [time1], "trait": [args.trait]})
    df_res.to_csv(args.output, index=False, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
