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


def simulation_sushie(rng_key, output_dic, N, L, L3, h2g, rho):
    mu = output_dic["mu"]
    L_cho = output_dic["L"]

    p, _ = L_cho[0].shape
    n_pop = len(N)

    X = []
    b_var = []

    for idx in range(n_pop):
        rng_key, x_key = random.split(rng_key, 2)
        tmp_X = L_cho[idx].dot((random.normal(x_key, shape=(N[idx], p)) + mu[idx]).T).T
        tmp_X -= jnp.mean(tmp_X, axis=0)
        tmp_X /= jnp.std(tmp_X, axis=0)
        X.append(tmp_X)
        b_var.append(h2g[idx] / (L + L3))

    b_covar = jnp.diag(jnp.array(b_var))
    ct = 0
    for row in range(1, n_pop):
        for col in range(n_pop):
            if col < row:
                _cov = jnp.sqrt(b_var[row] * b_var[col])
                b_covar = b_covar.at[row, col].set(rho[ct] * _cov)
                b_covar = b_covar.at[col, row].set(rho[ct] * _cov)
                ct += 1

    ct = 0
    # Make sure we have signficant region
    while ct < 10000:
        ct = ct + 1
        rng_key, b_key = random.split(rng_key, 2)

        bvec = random.multivariate_normal(b_key, jnp.zeros((n_pop,)), b_covar, shape=(L,))

        b_indep = jnp.zeros((n_pop, L3))
        for idx in range(n_pop):
            rng_key, b_indep_key = random.split(rng_key, 2)
            b_indep = b_indep.at[idx, :].set(jnp.sqrt(b_var[idx]) * random.normal(b_indep_key, shape=(L3,)))

        rng_key, gamma_key = random.split(rng_key, 2)
        gamma = random.choice(gamma_key, p, shape=(L + n_pop * L3,), replace=False)

        bvec_all = []
        for idx in range(n_pop):
            start_pos = L + L3 * idx
            end_pos = start_pos + L3 * (idx + 1)
            tmp_bvec = jnp.zeros(p).at[gamma[0:L]].set(jnp.transpose(bvec)[idx]).at[gamma[start_pos:end_pos]].set(
                b_indep[idx])
            bvec_all.append(tmp_bvec)

        g = []
        s2e = []
        y = []
        lrt = []
        for idx in range(n_pop):
            tmp_g = X[idx] @ bvec_all[idx]
            tmp_s2g = jnp.var(tmp_g)
            tmp_s2e = ((1 / h2g[idx]) - 1) * tmp_s2g
            rng_key, y_key = random.split(rng_key, 2)
            tmp_y = tmp_g + jnp.sqrt(tmp_s2e) * random.normal(y_key, shape=(N[idx],))
            _, _, _, tmp_lrt, _ = hf.estimate_her(X[idx][0:int(N[0] / 3), :], tmp_y[0:int(N[0] / 3)])

            g.append(tmp_g)
            s2e.append(tmp_s2e)
            y.append(tmp_y)
            lrt.append(tmp_lrt)

        # 2.7055
        keep_greater = True
        for idx in range(n_pop):
            keep_greater = False if lrt[idx] < 2.7055 else keep_greater

        if keep_greater:
            break

    return rng_key, X, y, bvec, bvec_all, s2e, gamma, b_covar


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
    argp.add_argument("--N", type=int, help="sample size for phenotype")
    argp.add_argument("--L1", default=2, type=int, help="number of generative L")
    argp.add_argument("--L2", default=2, type=int, help="number of inferential L")
    argp.add_argument("--L3", default=0, type=int, help="number of additional population-specific L")
    argp.add_argument("--h2g", type=float, help="signal to noise ratio for phenotype")
    argp.add_argument("--rho", type=float, default=0.8, help="cov parameter")
    argp.add_argument("--seed", default=1234, type=int, help="Random seed")
    argp.add_argument("--sim", default=1, help="Simulation Index")
    argp.add_argument("--locus", default=1, help="Locus Index")
    argp.add_argument("--threshold", default=0.9, help="Credible set threshold")
    argp.add_argument("-o", "--output", default=sys.stdout)

    args = argp.parse_args(args)
    rng_key = random.PRNGKey(int(args.seed))

    plink_path = args.prefix_pop.split(":")

    output_dic, _ = hf._gen_ld(plink_path)

    rng_key, X, y, bvec, bvec_all, s2e, g, b_covar = simulation_sushie(rng_key, output_dic, [args.N] * 3, args.L1,
                                                                       args.L3, [args.h2g] * 3, [args.rho] * 3)

    rng_key, p1, p2 = random.split(rng_key, 3)

    # perform sushie for 1 pop
    pop_idx = random.choice(p1, 3, shape=(1,), replace=False)
    sushie1 = infer_sushie([X[jnp.take(pop_idx, 0)]], [y[jnp.take(pop_idx, 0)]], L=args.L2, threshold=0.95)

    # perform sushie for 2 pop,

    pop_idx2 = random.choice(p2, 3, shape=(2,), replace=False)
    X_2pop = [X[jnp.take(pop_idx2, 0)][0:int(args.N / 2), :], X[jnp.take(pop_idx2, 1)][0:int(args.N / 2), :]]
    y_2pop = [y[jnp.take(pop_idx2, 0)][0:int(args.N / 2)], y[jnp.take(pop_idx2, 1)][0:int(args.N / 2)]]
    sushie2 = infer_sushie(X_2pop, y_2pop, L=args.L2, threshold=0.95)

    # perform sushie for 3 pop,
    X_3pop = [X[0][0:int(args.N / 3), :], X[1][0:int(args.N / 3), :], X[2][0:int(args.N / 3), :]]
    y_3pop = [y[0][0:int(args.N / 3)], y[1][0:int(args.N / 3)], y[2][0:int(args.N / 3)]]
    sushie3 = infer_sushie(X_3pop, y_3pop, L=args.L2, threshold=0.95)

    # PIP
    df_pip = pd.DataFrame()
    df_pip["CSIndex"] = jnp.argsort(-jnp.array(bvec ** 2)[:, 0]) + 1
    df_pip["SNPIndex_0based"] = g[0:args.L1]
    df_pip["SNPIndex_1based"] = g[0:args.L1] + 1
    df_pip["sushie1_pip"] = sushie1.pip_all[g[0:args.L1]]
    df_pip["sushie1_cali"] = jnp.isin(g[0:args.L1], sushie1.cs.SNPIndex.values.astype(int)).astype(int)
    df_pip["sushie2_pip"] = sushie2.pip_all[g[0:args.L1]]
    df_pip["sushie2_cali"] = jnp.isin(g[0:args.L1], sushie2.cs.SNPIndex.values.astype(int)).astype(int)
    df_pip["sushie3_pip"] = sushie3.pip_all[g[0:args.L1]]
    df_pip["sushie3_cali"] = jnp.isin(g[0:args.L1], sushie3.cs.SNPIndex.values.astype(int)).astype(int)

    df_pip = add_annotation(df_pip, args)
    df_pip.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.pip.tsv", sep="\t", index=False)

    sushie_cs1 = sushie1.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    sushie_cs2 = sushie2.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()
    sushie_cs3 = sushie3.cs.groupby("CSIndex")["SNPIndex"].count().reset_index()

    df_cs = pd.DataFrame(data=jnp.arange(args.L2) + 1, columns=["CSIndex"])
    df_cs = df_cs.merge(sushie_cs1, how="left", on="CSIndex").rename(columns={"SNPIndex": "sushie1"}) \
        .merge(sushie_cs2, how="left", on="CSIndex").rename(columns={"SNPIndex": "sushie2"}) \
        .merge(sushie_cs3, how="left", on="CSIndex").rename(columns={"SNPIndex": "sushie3"})

    df_cs = add_annotation(df_cs, args)
    df_cs.to_csv(f"{args.output}.sim{args.sim}.locus{args.locus}.cs.tsv", sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
