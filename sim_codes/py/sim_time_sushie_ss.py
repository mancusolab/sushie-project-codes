import argparse as ap
import sys
import warnings

import pandas as pd
from jax.config import config
from sushie.infer_ss import _infer_ss_test


config.update("jax_enable_x64", True)

# annoying warnings if on Mac ARM M1
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp


def main(args):
    argp = ap.ArgumentParser(description="")
    argp.add_argument("--ss", type=str, nargs="+")
    argp.add_argument("--ld", type=str, nargs="+")
    argp.add_argument("--trait", type=str)
    argp.add_argument("--output", type=str)

    args = argp.parse_args(args)

    ss = []
    ld = []
    for idx in range(2):
        tmp_ss = pd.read_csv(args.ss[idx], sep="\t")
        tmp_ld = jnp.array(pd.read_csv(args.ld[idx], sep="\t", header=None))
        ss.append(tmp_ss.T_STAT.values)
        ld.append(tmp_ld)

    z_ss = jnp.vstack(ss)
    lds = jnp.stack(ld, axis=0)


    import time
    start_time = time.time()
    _infer_ss_test(zs=z_ss, lds=lds, ns=jnp.array([398, 297])[:, jnp.newaxis], L=10)
    end_time = time.time()

    time1 = (end_time - start_time)

    df_res = pd.DataFrame({"method": ["sushie_ss"], "time": [time1], "trait": [args.trait]})
    df_res.to_csv(args.output, index=False, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
