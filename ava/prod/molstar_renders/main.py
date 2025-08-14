"""
convert to ipynb, then convert to html with:

jupyter nbconvert ava\stage\molstar_renders\main.ipynb --to html --execute --no-input

"""

# %%

from ipymolstar import PDBeMolstar
from hal.config import cfg
from hal.repro import reproduce
from pyhdx.fileIO import csv_to_dataframe
from pyhdx.plot import CMAP_NORM_DEFAULTS
from ultraplot import to_hex
import pandas as pd
from ipywidgets import VBox, Label, Layout

# %%

packages = ["pyhdx", "hdxms_datasets"]
OUTPUT_PATH = reproduce(globals(), packages=packages)
DG_CMAP, DG_NORM = CMAP_NORM_DEFAULTS["dG"]
DG_PTH = cfg.prod_output_path("bottom_up_deltaGs")


# %%


def make_color_data(df: pd.DataFrame):
    data = []
    for i in range(len(df)):
        row = df.reset_index().iloc[i, :].to_dict()
        r_number = row["r_number"]

        hex_color = to_hex(DG_CMAP(DG_NORM(row["dG"])), keep_alpha=False)
        elem = {"residue_number": r_number, "color": hex_color}
        data.append(elem)

    color_data = {"data": data, "nonSelectedColor": None}

    return color_data


def make_tooltips(r_number, values, qty="ΔG", unit="kJ/mol"):
    tooltips = {
        "data": [
            {"residue_number": int(resi), "tooltip": f"{qty}: {value:.2f} {unit}"}
            for resi, value in zip(r_number, values)
        ]
    }
    return tooltips


# %%
protein_names = ["MBP_p101", "MBP_wt", "MBP_Y283D"]
frames = {}
for protein_name in protein_names:
    df = csv_to_dataframe(DG_PTH / f"gibbs_fit_{protein_name}" / "fit_result.csv")[
        protein_name
    ].reset_index()
    frames[protein_name] = df

columns = {}
for protein_name, df in frames.items():
    color_data = make_color_data(df)
    tooltips = make_tooltips(df["r_number"], df["dG"] * 1e-3)

    view = PDBeMolstar(
        molecule_id="1anf",
        hide_water=True,
        hide_carbs=True,
        color_data=color_data,
        tooltips=tooltips,
        layout=Layout(width="100%"),
    )
    title = Label(value=protein_name)
    columns[protein_name] = VBox([title, view], layout=Layout(width="32.5%"))
