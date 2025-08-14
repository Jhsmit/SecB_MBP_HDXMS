# %%
from pathlib import Path

import colorcet as cc
import pandas as pd
import ultraplot as uplt
from pyhdx.datasets import DataVault
from pyhdx.models import HDXMeasurement
from pyhdx.plot import CMAP_NORM_DEFAULTS
import numpy as np

from hal.config import cfg
from hal.fmt import load_pplt_config
from hal.repro import reproduce
from ava.toolbox.plot.constants import SAVE_KWARGS
from ava.toolbox.plot.utils import colormesh
from ava.prod.figures_folding_per_residue.plot_data import PlotData


# %%
DS_NAME = "20221207_1304_MBP_folding_Karamanou"
# DG_DIR = "bottom_up_deltaGs"
DG_CUTOFF = 17500  # dg J/mol below which are considered flexible
FOLDON_CUTOFF = 30.000002  # time in minutes below which foldons are defined
DOF_CMAP = uplt.Colormap(cc.cm.gouldian, reverse=True)
DOF_CMAP.set_bad("#dedede")
DOF_NORM = uplt.Norm("linear", 0.0, 1.0)
FOLDON_COLOR = "r"
DG_CMAP, DG_NORM = CMAP_NORM_DEFAULTS["dG"]

REFASPECT = 12  # aspect ratio bars


# %%
OUTPUT_DIR = reproduce(globals())
T_MAX_HOURS = 1
SETTINGS = {
    "MBP_wt_p101": {"bottom_up_state": "MBP_p101", "tmax": 60.0 * T_MAX_HOURS},
    "MBP_wt_p101_2022": {"bottom_up_state": "MBP_p101", "tmax": 60.0 * T_MAX_HOURS},
    "MBP_wt_p101_SecB": {"bottom_up_state": "MBP_p101", "tmax": 60.0 * T_MAX_HOURS},
    "MBP_Y283D": {"bottom_up_state": "MBP_Y283D", "tmax": 60.0 * T_MAX_HOURS},
    "MBP_Y283D_SecB": {"bottom_up_state": "MBP_Y283D", "tmax": 60.0 * T_MAX_HOURS},
}

load_pplt_config("paper")

custom_labels = [
    ("5 s", 5),
    ("10 s", 10),
    ("20 s", 20),
    # ("40 s", 40),
    ("1 min", 60),
    ("2.5 min", 150),
    ("5 min", 300),
    ("10 min", 600),
    # ("15 min", 900),
    # ("20 min", 1200),
    ("30 min", 1800),
    ("1 h", 3600),
    ("3 h", 10800),
    ("20 h", 72000),
]
tick_label, tick_pos = zip(*[(label, pos / 60) for label, pos in custom_labels])


# %%
def convert_save(df: pd.DataFrame, path: Path) -> None:
    df_copy = df.copy()
    df_copy.columns = df_copy.columns / 60.0

    df_copy.to_csv(path)


# %%
vault = DataVault(cache_dir=cfg.root / "data")
ds = vault.load_dataset(DS_NAME)

mydata = {}
for state in ds.states:
    data = PlotData(
        hdxm=HDXMeasurement.from_dataset(ds, state),
        settings=SETTINGS[state],  # type: ignore
    )
    mydata[state] = data


# %%

for state in ds.states:
    data: PlotData = mydata[state]

    # save dof data to file
    convert_save(data.hdxm.rfu_residues, OUTPUT_DIR / f"{state}_dof.csv")
    convert_save(data.hdxm.rfu_residues_sd, OUTPUT_DIR / f"{state}_dof_std.csv")

    # make the figure
    fig, ax = uplt.subplots(refwidth="100mm", refheight="30mm", sharex=False)
    mesh = ax.pcolormesh(
        data.dof_df.transpose(), edgefix=False, cmap=DOF_CMAP, N=256, vmin=0, vmax=1
    )
    ax.plot(data.get_tfold(0.50), color=FOLDON_COLOR, lw=1)

    px_foldon = ax.panel_axes("t", width="3mm", space="1mm", share=False)
    px_foldon.pcolormesh(
        data.r_edges, [0, 1], data.foldons_by_integration, cmap="greys"
    )
    px_foldon.format(
        yticks=[],
        xlim=data.xlim,
        xticks=[],
    )

    if "secb" not in state.lower():
        px_dg = ax.panel_axes("t", width="3mm", space="0mm", share=False)
        px_dg.pcolor(
            data.r_edges,
            [0, 1],
            data.dg_df["dG"].to_numpy().reshape(1, -1),
            cmap=DG_CMAP,
            vmin=DG_NORM.vmin,
            vmax=DG_NORM.vmax,
        )
        px_dg.format(
            yticks=[],
            xlim=data.xlim,
            xticks=[],
        )

    # cbar = ax.colorbar(DOF_CMAP, locator=pplt.Locator("fixed", pplt.arange(0, 1, 0.2)), width="2mm")
    # cbar.set_label("Degree of Folding")

    ax_d = ax.dualy(lambda x: x, yformatter=uplt.Formatter("log"))
    ax_d.format(yticklabels=[])
    ax.format(
        ylim=data.ylim,
        xlim=data.xlim,
        ylabel="Folding time",
        xlabel="Residue Number",
        yscale="log",
        yticks=tick_pos,
        yticklabels=tick_label,
        ytickminor=False,
        grid=False,
    )

    fig.savefig(OUTPUT_DIR / f"{state}_folding.png", **SAVE_KWARGS)
    fig.savefig(OUTPUT_DIR / f"{state}_folding.pdf", **SAVE_KWARGS)
    uplt.close(fig)

# %%


def get_column(columns, value):
    idx = np.argmin([abs(col - value) for col in columns])
    return columns[idx]


BAR_DIR = OUTPUT_DIR / "bars"
BAR_DIR.mkdir(exist_ok=True)

states_and_exposures = [
    ("MBP_wt_p101_2022", 10.0 / 60),
    ("MBP_wt_p101_2022", 2.5),
    ("MBP_wt_p101_SecB", 2.5),
    ("MBP_Y283D", 10.0),
    ("MBP_Y283D", 30.0),
    ("MBP_Y283D_SecB", 30.0),
]

for state, exposure in states_and_exposures:
    df = mydata[state].dof_df
    col = get_column(df.columns, exposure)

    values = df[col]
    values = values.reindex(pd.RangeIndex(0, values.index.max() + 1))

    fig, ax = uplt.subplots(refaspect=REFASPECT, axheight="10mm")
    mesh = colormesh(ax, values, DOF_CMAP, DOF_NORM)

    ax.format(
        yticks=[],
        xlabel="Residue Number",
        title=f"{state} at {exposure:.3f} min",
        grid=False,
    )
    ax.colorbar(
        mesh, loc="b", width="2mm", ticks=uplt.arange(0.0, 1.0, 0.25), label="DOF"
    )

    fig.savefig(BAR_DIR / f"BAR_{state}_{exposure:.3f}.png", **SAVE_KWARGS)
    fig.savefig(BAR_DIR / f"BAR_{state}_{exposure:.3f}.pdf", **SAVE_KWARGS)

# %%
