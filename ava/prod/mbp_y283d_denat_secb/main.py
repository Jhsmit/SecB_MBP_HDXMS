# %%

import colorcet as cc
import numpy as np
import pandas as pd
import ultraplot as uplt
from cmap import Colormap
from pyhdx.datasets import DataVault
from pyhdx.models import HDXMeasurement
from pyhdx.plot import CMAP_NORM_DEFAULTS, linear_bars

from hal.config import cfg
from hal.fmt import load_pplt_config
from hal.repro import reproduce
from ava.toolbox.constants import EXCLUDE_N_TERM
from ava.toolbox.plot.constants import SAVE_KWARGS
from ava.toolbox.plot.utils import colormesh

# %%
uplt.__version__
# %%
DS_NAME = "1729844783_Y283D_SecB_denat_karamanou"
FOLDON_CUTOFF = 30.000002  # time in minutes below which foldons are defined
DOF_CMAP = uplt.Colormap(cc.cm.gouldian, reverse=True)
NO_COVERAGE = "#dedede"
REFASPECT = 12

DOF_CMAP.set_bad(NO_COVERAGE)
DOF_NORM = uplt.Norm("linear", 0, 1)

DDOF_NORM = uplt.Norm("linear", vmin=-0.25, vmax=0.25)
DDOF_CMAP = Colormap("tol:sunset").to_matplotlib()
DDOF_CMAP.set_bad(NO_COVERAGE)

DG_CMAP, DG_NORM = CMAP_NORM_DEFAULTS["dG"]


# %%
OUTPUT_DIR = reproduce(globals())
load_pplt_config("paper")

vault = DataVault(cache_dir=cfg.root / "data")
ds = vault.load_dataset(DS_NAME)

# %%
mydata = {state: HDXMeasurement.from_dataset(ds, state) for state in ds.states}

# these data have only one column
# and likely their exposure labels are not correct
series_dict = {}
for state in ds.states:
    df = mydata[state].rfu_residues
    df.loc[EXCLUDE_N_TERM[0] : EXCLUDE_N_TERM[1], :] = np.nan
    s = df.iloc[:, 0]
    series_dict[state] = s


# %%
data_dict = {"SecB denaturase Y283D": series_dict}
fig, ax, cbar = linear_bars(data_dict, uplt.Norm("linear", 0, 1), DOF_CMAP)
cbar.set_label("Degree of Folding")

fig.savefig(OUTPUT_DIR / "Y283D_SecB_folding_norm_1.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_folding_norm_1.pdf", **SAVE_KWARGS)

# %%
data_dict = {"SecB denaturase Y283D": series_dict}
fig, ax, cbar = linear_bars(data_dict, uplt.Norm("linear", 0, 0.5), DOF_CMAP)
cbar.set_label("Degree of Folding")

fig.savefig(OUTPUT_DIR / "Y283D_SecB_folding_norm_05.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_folding_norm_05.pdf", **SAVE_KWARGS)

# %%

fig, ax = uplt.subplots(refaspect=3.5, axheight="75mm")
for state in ds.states:
    s = series_dict[state]
    ax.plot(s.index, s, label=state)
ax.format(xlabel="Residue Number", ylabel="DOF")
ax.legend()

fig.savefig(OUTPUT_DIR / "Y283D_SecB_DOF_lines.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_DOF_lines.pdf", **SAVE_KWARGS)

# %%
fig, ax = uplt.subplots(refaspect=3.5, axheight="75mm")
for state in ds.states:
    if state == "folding_secb":
        continue
    s = series_dict[state]
    ax.plot(s.index, s, label=state)
ax.format(xlabel="Residue Number", ylabel="DOF")
ax.legend()

fig.savefig(OUTPUT_DIR / "Y283D_SecB_DOF_lines_denaturase.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_DOF_lines_denaturase.pdf", **SAVE_KWARGS)


# %%

ref = "folding"
others = set(ds.states) - set([ref])

fig, ax = uplt.subplots(refaspect=3.5, axheight="50mm")
s_ref = series_dict[ref]
for state in others:
    s = series_dict[state]
    d = s - s_ref
    ax.plot(d.index, d, label=state)
ax.axhline(0, color="#8c8c8c", ls="--")
ax.format(
    xlabel="Residue Number",
    ylabel="Denaturase (ΔDOF)",
    title="SecB Denaturase activity",
)
ax.legend()

fig.savefig(OUTPUT_DIR / "Y283D_SecB_delta_DOF_lines.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_delta_DOF_lines.pdf", **SAVE_KWARGS)

# %%

ref = "folding"
others = set(ds.states) - set([ref])

fig, ax = uplt.subplots(refaspect=3.5, axheight="50mm")
s_ref = series_dict[ref]
for state in others:
    if state == "folding_secb":
        continue
    s = series_dict[state]
    d = s - s_ref
    ax.plot(d.index, d, label=state)
ax.axhline(0, color="#8c8c8c", ls="--")
ax.format(
    xlabel="Residue Number",
    ylabel="Denaturase (ΔDOF)",
    title="SecB Denaturase activity",
)
ax.legend()

fig.savefig(OUTPUT_DIR / "Y283D_SecB_delta_DOF_lines_denaturase.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_delta_DOF_lines_denaturase.pdf", **SAVE_KWARGS)

# %%


fig, ax = uplt.subplots(refaspect=REFASPECT, axheight="10mm")

test_state = "folding_then_secb"
s = series_dict[test_state]
d = s - s_ref
values = d.reindex(pd.RangeIndex(0, d.index.max() + 1))

mesh = colormesh(ax, values, DDOF_CMAP, DDOF_NORM)
ax.format(yticks=[], xlabel="Residue Number")
ax.colorbar(mesh, loc="b", width="2mm", ticks=uplt.arange(-0.2, 0.2, 0.1), label="ΔDOF")

fig.savefig(OUTPUT_DIR / "Y283D_SecB_delta_DOF_colorbar_denaturase.png", **SAVE_KWARGS)
fig.savefig(OUTPUT_DIR / "Y283D_SecB_delta_DOF_colorbar_denaturase.pdf", **SAVE_KWARGS)

# %%

for state in ds.states:
    fig, ax = uplt.subplots(refaspect=REFASPECT, axheight="10mm")
    s = series_dict[state]
    values = s.reindex(pd.RangeIndex(0, s.index.max() + 1))

    mesh = colormesh(ax, values, DOF_CMAP, DOF_NORM)
    ax.format(yticks=[], xlabel="Residue Number")
    ax.colorbar(mesh, loc="b", width="2mm", ticks=uplt.arange(0, 1, 0.25), label="DOF")

    fig.savefig(OUTPUT_DIR / f"Y283D_SecB_DOF_{state}.png", **SAVE_KWARGS)
    fig.savefig(OUTPUT_DIR / f"Y283D_SecB_DOF_{state}.pdf", **SAVE_KWARGS)


export = {}
for state in ds.states:
    export[f"{state}_dof"] = mydata[state].rfu_residues.iloc[:, 0]
    export[f"{state}_dof_sd"] = mydata[state].rfu_residues_sd.iloc[:, 0]


# %%
df = pd.DataFrame(export)
df.to_csv(OUTPUT_DIR / "Y283D_SecB_folding.csv")
# %%
