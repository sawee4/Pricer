#!/usr/bin/env python3
"""Plot the convergence studies written by the example programs.

Usage (from the repository root, after building):

    cd docs/data
    ../../build/examples/convergence_study
    ../../build/examples/variance_reduction_study
    ../../build/examples/binomial_convergence
    ../../build/examples/mc_greeks_study
    ../../build/examples/exotics_study
    cd ../..
    python scripts/plot_results.py            # reads docs/data, writes docs/figures

Requires matplotlib only (the CSV files are read with the standard library).
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt  # noqa: E402

# Categorical colours in a fixed order (colour-blind safe on a light background),
# plus a neutral grey for reference lines (theory, closed forms).
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
REFERENCE = "#52514e"
TEXT = "#0b0b0b"
GRID = "#e4e3df"

plt.rcParams.update({
    "figure.figsize": (6.4, 4.0),
    "figure.dpi": 150,
    "savefig.bbox": "tight",
    "font.size": 10,
    "axes.edgecolor": REFERENCE,
    "axes.labelcolor": TEXT,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": REFERENCE,
    "ytick.color": REFERENCE,
    "lines.linewidth": 2.0,
    "lines.markersize": 5,
    "legend.frameon": False,
})


def read_csv(path):
    """Return the CSV as a dict of column name -> list of floats."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{path} is empty")
    return {key: [float(row[key]) for row in rows] for key in rows[0]}


def save(fig, out_dir, name):
    path = out_dir / name
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def plot_mc_convergence(data_dir, out_dir):
    d = read_csv(data_dir / "mc_convergence.csv")
    n, rmse, se = d["n_paths"], d["rmse"], d["mean_std_error"]
    fig, ax = plt.subplots()
    # Reference line c / sqrt(N) anchored on the first reported standard error.
    c = se[0] * n[0] ** 0.5
    ax.loglog(n, [c / x ** 0.5 for x in n], "--", color=REFERENCE, linewidth=1.2,
              label=r"$\propto 1/\sqrt{N}$")
    ax.loglog(n, se, "o-", color=BLUE, label="mean reported standard error")
    ax.loglog(n, rmse, "s-", color=ORANGE, label="RMSE vs Black-Scholes")
    ax.set_xlabel("number of paths $N$")
    ax.set_ylabel("error")
    ax.set_title("Monte Carlo convergence (Hull call, 50 seeds per $N$)")
    ax.legend()
    save(fig, out_dir, "mc_convergence.png")


def plot_variance_reduction(data_dir, out_dir):
    d = read_csv(data_dir / "variance_reduction.csv")
    fig, ax = plt.subplots()
    ax.semilogy(d["strike"], d["factor_control"], "o-", color=BLUE,
                label=r"control variate $e^{-rT}S_T$")
    ax.semilogy(d["strike"], d["factor_antithetic"], "s-", color=ORANGE,
                label="antithetic variates")
    ax.axhline(1.0, color=REFERENCE, linestyle="--", linewidth=1.2, label="plain Monte Carlo")
    ax.axvline(100.0, color=GRID, linewidth=1.0)
    ax.set_xlabel("strike $K$ (spot 100)")
    ax.set_ylabel(r"variance reduction factor $\mathrm{Var}_{plain}/\mathrm{Var}$")
    ax.set_title("Variance reduction across moneyness ($N$ = 100 000)")
    ax.legend()
    save(fig, out_dir, "variance_reduction.png")


def plot_binomial(data_dir, out_dir):
    d = read_csv(data_dir / "binomial_convergence.csv")
    fig, ax = plt.subplots()
    ax.plot(d["n_steps"], d["error"], color=BLUE, linewidth=1.0, label="CRR tree, $n$ steps")
    ax.plot(d["n_steps"], d["error_avg_n_n1"], color=ORANGE, linewidth=1.5,
            label="average of $n$ and $n+1$ steps")
    ax.axhline(0.0, color=REFERENCE, linewidth=1.0)
    ax.set_xlabel("number of steps $n$")
    ax.set_ylabel("price error vs Black-Scholes")
    ax.set_ylim(-0.06, 0.06)
    ax.set_title("Binomial tree convergence (European call, K = 105)")
    ax.legend()
    save(fig, out_dir, "binomial_convergence.png")


def plot_mc_greeks(data_dir, out_dir):
    d = read_csv(data_dir / "mc_greeks_bump.csv")
    fig, ax = plt.subplots()
    ax.loglog(d["relative_bump"], d["std_delta_independent"], "s-", color=ORANGE,
              label="independent draws (blows up like $1/h$)")
    ax.loglog(d["relative_bump"], d["std_delta_crn"], "o-", color=BLUE,
              label="common random numbers")
    ax.invert_xaxis()
    ax.set_xlabel("relative bump $h/S$")
    ax.set_ylabel("std of the delta estimate (50 seeds)")
    ax.set_title("Bump-and-revalue delta: noise vs bump size")
    ax.legend()
    save(fig, out_dir, "mc_greeks_bump.png")


def plot_barrier(data_dir, out_dir):
    d = read_csv(data_dir / "barrier_monitoring.csv")
    m, price, se = d["monitoring_dates"], d["mc_price"], d["mc_std_error"]
    fig, ax = plt.subplots()
    ax.axhline(d["continuous_price"][0], color=REFERENCE, linestyle="--", linewidth=1.2,
               label="continuous monitoring (closed form)")
    ax.semilogx(m, d["bgk_price"], "-", color=ORANGE, marker="D",
                label="closed form at BGK-shifted barrier")
    ax.errorbar(m, price, yerr=[1.96 * s for s in se], fmt="o", color=BLUE, capsize=4,
                label="Monte Carlo, discrete monitoring (95% CI)")
    ax.set_xticks(m)
    ax.set_xticklabels([str(int(x)) for x in m])
    ax.minorticks_off()
    ax.set_xlabel("monitoring dates $M$")
    ax.set_ylabel("price")
    ax.set_title("Down-and-out call (K = 100, H = 90)")
    ax.legend()
    save(fig, out_dir, "barrier_monitoring.png")


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=root / "docs" / "data")
    parser.add_argument("--out", type=Path, default=root / "docs" / "figures")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for plot in (plot_mc_convergence, plot_variance_reduction, plot_binomial, plot_mc_greeks,
                 plot_barrier):
        plot(args.data, args.out)


if __name__ == "__main__":
    main()
