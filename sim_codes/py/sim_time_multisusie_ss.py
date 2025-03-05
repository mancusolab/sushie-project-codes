import argparse as ap
import sys
import warnings

import pandas as pd
import numpy as np
import MultiSuSiE

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
        tmp_ld = np.array(pd.read_csv(args.ld[idx], sep="\t", header=None))
        ss.append(tmp_ss.T_STAT.values)
        ld.append(tmp_ld)

    import time
    start_time = time.time()
    MultiSuSiE.outer_wrapper(z_list=[np.array(ss[0]), np.array(ss[1])],
                              R_list=ld, rho=np.array([[1, 0.1], [0.1, 1]]),
                              population_sizes = [398, 297], L = 10, max_iter = 500, tol = 0.0001,
                              min_abs_corr = 0.5)
    end_time = time.time()

    time1 = (end_time - start_time)

    df_res = pd.DataFrame({"method": ["multisusie_ss"], "time": [time1], "trait": [args.trait]})
    df_res.to_csv(args.output, index=False, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
