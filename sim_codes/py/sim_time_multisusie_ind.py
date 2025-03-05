import argparse as ap
import sys
import warnings

import pandas as pd
import numpy as np
import MultiSuSiE


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
        tmp_pheno = np.array(pd.read_csv(args.pheno[idx], sep="\t", header=None))
        tmp_geno = np.array(pd.read_csv(args.geno[idx], sep="\t", header=None))
        pheno.append(np.squeeze(tmp_pheno))
        geno.append(tmp_geno)

    import time
    start_time = time.time()
    MultiSuSiE.ind_wrapper(X_list=[np.array(geno[0]), np.array(geno[1])],Y_list=[np.array(pheno[0]), np.array(pheno[1])],rho=np.array([[1, 0.1], [0.1, 1]]),L=10,standardize=True,intercept=False,float_type=np.float64)
    end_time = time.time()

    time1 = (end_time - start_time)

    df_res = pd.DataFrame({"method": ["multisusie_ind"], "time": [time1], "trait": [args.trait]})
    df_res.to_csv(args.output, index=False, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
