"""
overall folding of MBP is slow (10 fold compared to smFRET apparent folding rate)
however, the foldons are faster, and they are near smFRET rates

the foldons and bulk protein folding is a stretched exponential. This indicates
more heterogeneity in the folding pathway. Possibly HDX-MS shows more
heterogeneity as it can resolve more folding microstates
ie a folding intermediates may have identical smFRET efficiencies while
these do have differences in secondary structure and therefore HDX-MS
uptake.

These stretched exponentials suggest heterogeneous folding pathways which
occur from a folding pathway with many kinetic traps (ie prolines, misfolds, etc)

thus the foldon folding is heterogeneous but the bulk folding is
more characterized by a single rate-determining step after which all the remaining folding cascades

See commentary:
https://www.pnas.org/doi/10.1073/pnas.2009596117

"""

# %%
from dataclasses import dataclass
from smitfit.model import Model
from uncertainties import ufloat

import numpy as np
import pandas as pd
import ultraplot as uplt
import scipy
from pyhdx.datasets import DataVault
from pyhdx.models import HDXMeasurement

from sympy import exp

from ava.toolbox.plot.constants import Y_FMT_KWARGS
from hal.config import cfg
from hal.repro import reproduce
from ava.toolbox.constants import EXCLUDE_N_TERM
from smitfit.result import Result as FitResult
from smitfit.symbol import Symbols
from smitfit import Function
from smitfit.lmfit import Minimize
from hal.io import Output, save_csv, save_fig, save_yaml
from cmap import Colormap

# %%
OUTPUT_DIR = reproduce(globals())

# data loading
DS_NAME = "20221207_1304_MBP_folding_Karamanou"
vault = DataVault(cache_dir=cfg.root / "data")
DATASET = vault.load_dataset(DS_NAME)

# plotting
XLIM = (5e-2, 3e5)
T_EVAL = np.logspace(np.log10(XLIM[0]), np.log10(XLIM[1]), 100)

# intervals are inclusive, exclusive
# from `figures_folding_per_residue.py`
FOLDONS = [(8, 19), (30, 38), (54, 68), (107, 114), (266, 279)]
FOLDON_COLORS = [
    "#ff0000",
    "#ffc000",
    "#ed7d31",
    "#ff8ad8",
    "#a9d18e",
]
FOLDON_LABELS = ["β", "γ", "α", "ε", "δ"]

DATASET.states
STATES = ["MBP_wt_p101_2022", "MBP_wt_p101_SecB", "MBP_Y283D", "MBP_Y283D_SecB"]
OUTPUT_FILES = [
    "datapoints.csv",
    "fitted_curves.csv",
    "fit_results.yaml",
    "subplots_graph.png",
    "combined_graph.png",
    "combined_graph_linear.png",
]
OVERWRITE = True


@dataclass
class Result:
    foldon: tuple[int, int]
    fit_result: FitResult
    fit_data: dict


# %%


def mean_relaxation_time(result: FitResult) -> float:
    tau = 1 / result.parameters["k"]
    beta = result.parameters["beta"]
    tau_mean = (tau / beta) * scipy.special.gamma(1 / beta)
    return tau_mean


def annotate(ax: uplt.Axes, result: FitResult):
    k = result.parameters["k"]
    beta = result.parameters["beta"]
    tau = mean_relaxation_time(result)

    annotation = f"k = {k:.2e} s⁻¹\nβ = {beta:.2f}\n\nτ = {tau:.2f} s"
    ax.text(
        0.05,
        0.95,
        annotation,
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),  # type: ignore
    )


def make_graph(
    df_datapoints: pd.DataFrame, df_fits: pd.DataFrame, results: list[Result]
):
    fig, ax = uplt.subplots()
    for i, result in enumerate(results):
        ax.scatter(
            df_datapoints.index,
            df_datapoints[str(result.foldon)],
            label=FOLDON_LABELS[i],
            color=colors[i],
        )
        ax.plot(
            df_fits.index,
            df_fits[str(result.foldon)].to_numpy(),
            color=colors[i],
        )

    ax.scatter(df_datapoints.index, df_datapoints["bulk"], color="#5d5d5d")
    ax.plot(df_fits.index, df_fits["bulk"].to_numpy(), color="#5d5d5d")

    return fig, ax


def make_table(combined_results):
    items = []
    params = ["b", "k", "beta"]
    for foldon_label in sorted(FOLDON_LABELS):
        idx = FOLDON_LABELS.index(foldon_label)
        foldon = FOLDONS[idx]
        item = {"name": foldon_label}
        key = str(foldon)
        for param in params:
            result = combined_results[key]
            param_dict = result["fit_parameters"] | result["fixed_parameters"]
            uf = ufloat(param_dict[param], result["errors"].get(param, np.nan))
            item[param] = "{:.3g}".format(uf)
        items.append(item)

    key = "bulk"
    item = {"name": "bulk"}
    for param in params:
        result = combined_results[key]
        param_dict = result["fit_parameters"] | result["fixed_parameters"]
        uf = ufloat(param_dict[param], result["errors"].get(param, np.nan))
        item[param] = "{:.3g}".format(uf)
    items.append(item)

    df = pd.DataFrame(items).set_index("name")

    return df


# %%
s = Symbols("a b beta k t y")
rhs = s.a + (s.b - s.a) * (1 - exp(-((s.k * s.t) ** s.beta)))  # type: ignore
f = Function(rhs)
model = Model({s.y: rhs})
guess = dict(a=0.1, b=0.9, k=0.001, beta=1)

# a: 0.18302442194111535
# b: 0.8728527253417224

parameters = f.define_parameters(guess)
parameters["a"].fix()
parameters["b"].set_bounds(lower_bound=0, upper_bound=1)
parameters["k"].set_bounds(lower_bound=1e-5, upper_bound=1e1)

cmap = Colormap(f"colorbrewer:dark2_{len(FOLDONS)}")
colors = [cmap(i) for i in range(len(FOLDONS))]

# %%
p101_mean_b = None
for state in DATASET.states:
    print(state)
    if state == "MBP_wt_p101":  # this is old data
        continue

    # fit_params = parameters.copy()
    # if "wt" in state:
    #     fit_params["b"].set_guess(0.9)
    # elif "Y283D" in state:
    #     fit_params["b"].set_guess(1.0)

    save_dir = OUTPUT_DIR / state
    output = Output(save_dir, files=OUTPUT_FILES, overwrite=OVERWRITE)

    data = {
        state: HDXMeasurement.from_dataset(DATASET, state) for state in DATASET.states
    }
    hdxm = data[state]

    df = hdxm.rfu_residues
    drop_range = range(EXCLUDE_N_TERM[0], EXCLUDE_N_TERM[1] + 1)
    df = df.drop(drop_range)

    # Create a boolean mask, initially all True
    mask = np.ones(len(df), dtype=bool)
    for start, end in FOLDONS:
        # Update mask to False for indices in the interval
        mask[start:end] = False

    # take only residues which are not in the foldons
    bulk_df = df[mask].dropna()
    bulk_mean = bulk_df.mean(axis=0)

    b_vals = []
    xdata, ydata = dict(t=np.array(bulk_mean.index)), dict(y=np.array(bulk_mean.values))

    foldon_parameters = parameters.copy()
    if state == "MBP_wt_p101_SecB":
        assert p101_mean_b is not None
        foldon_parameters["b"].set_guess(p101_mean_b).fix()

    bulk_result = Minimize(model, foldon_parameters, xdata, ydata).fit()
    b_vals.append(bulk_result.parameters["b"])

    # foldon_parameters["a"].fix().set_guess(bulk_result.parameters["a"])
    # for p101 secb, we set the value of parameter b to the mean of p101 folding

    results: list[Result] = []
    for foldon in FOLDONS:
        section = df.loc[foldon[0] : foldon[1] - 1]  # pandas indexing is weird
        values = section.dropna().mean(axis=0)

        xdata, ydata = dict(t=np.array(values.index)), dict(y=np.array(values.values))
        result = Minimize(model, foldon_parameters, xdata, ydata).fit()
        b_vals.append(result.parameters["b"])
        results.append(
            Result(foldon=foldon, fit_result=result, fit_data={**xdata, **ydata})
        )

    if state == "MBP_wt_p101_2022":
        p101_mean_b = np.mean(b_vals)

    # %%
    # create dataframe of fit data to plot
    df_datapoints = pd.DataFrame(
        {str(result.foldon): result.fit_data["y"] for result in results},
        index=results[0].fit_data["t"],
    )
    df_datapoints["bulk"] = bulk_mean.values
    save_csv(df_datapoints, output["datapoints.csv"])

    # %%
    # create dataframe of fitted curves
    df_fits = pd.DataFrame(index=pd.Index(T_EVAL, name="t"))
    for result in results:
        y_eval = f(**result.fit_result.parameters, t=T_EVAL)
        df_fits[str(result.foldon)] = y_eval

    df_fits["bulk"] = f(**bulk_result.parameters, t=T_EVAL)
    save_csv(df_fits, output["fitted_curves.csv"])

    # create subplot figure with all individual curves
    fig, axes = uplt.subplots(ncols=3, nrows=2)
    for ax, result in zip(axes, results):
        ax.scatter(df_datapoints.index, df_datapoints[str(result.foldon)])
        ax.plot(df_fits.index, df_fits[str(result.foldon)], color="k")

        annotate(ax, result.fit_result)
        ax.format(title=str(result.foldon))

    axes[-1].scatter(df_datapoints.index, df_datapoints["bulk"])
    axes[-1].plot(df_fits.index, df_fits["bulk"], color="k")

    annotate(axes[-1], bulk_result)
    axes.format(
        xlabel="Folding time (s)", ylabel="DOF", xscale="log", xlim=XLIM, ylim=(0, 1)
    )
    save_fig(fig, output["subplots_graph.png"])

    fig, ax = make_graph(df_datapoints, df_fits, results)
    ax.format(
        xscale="log",
        xlim=XLIM,
        ylim=(0, 1),
        xlabel="Folding time (s)",
        ylabel="DOF",
        **Y_FMT_KWARGS,
    )
    ax.legend(loc="r", ncols=1)
    save_fig(fig, output["combined_graph.png"])

    fig, ax = make_graph(df_datapoints, df_fits, results)
    ax.format(
        xlim=(0, 900),
        ylim=(0, 1),
        xlabel="Folding time (s)",
        ylabel="DOF",
        **Y_FMT_KWARGS,
    )
    ax.legend(loc="r", ncols=1)
    save_fig(fig, output["combined_graph_linear.png"])

    combined_results = {
        "bulk": bulk_result.to_dict(),
        **{str(result.foldon): result.fit_result.to_dict() for result in results},
    }
    save_yaml(combined_results, output["fit_results.yaml"])

    df = make_table(combined_results)
    s = df.to_html(save_dir / "result_table.html")

# %%
first_vals = [0.00375, 0.00434, 0.00383]
first_errors = [0.00146, 0.00116, 0.00076]

ufs = [ufloat(v, e) for v, e in zip(first_vals, first_errors)]
np.mean(ufs)
# 0.003973333333333333+/-0.0006712177987310328

# %%
second_vals = [0.00238, 0.00220]
second_errors = [0.00040, 0.00045]

ufs = [ufloat(v, e) for v, e in zip(second_vals, second_errors)]
np.mean(ufs)
# 0.0022900000000000004+/-
# 0.0003010398644698074
# %%
