# %%
from dataclasses import dataclass
from functools import cached_property
from typing import TypedDict

import numpy as np
import pandas as pd
from pyhdx.fileIO import csv_to_dataframe
from pyhdx.models import HDXMeasurement

from hal.config import cfg
from ava.toolbox.constants import EXCLUDE_N_TERM

# %%

DG_CUTOFF = 17500  # dg J/mol below which are considered flexible


DS_NAME = "20221207_1304_MBP_folding_Karamanou"
FOLDON_CUTOFF = 30.000002  # time in minutes below which foldons are defined


# DOF_CMAP = uplt.Colormap(cc.cm.gouldian, reverse=True)
# DOF_CMAP.set_bad("#dedede")
# FOLDON_COLOR = "r"
# DG_CMAP, DG_NORM = CMAP_NORM_DEFAULTS["dG"]

# %%


class Settings(TypedDict):
    bottom_up_state: str
    tmax: float


@dataclass(frozen=True)
class PlotData:
    hdxm: HDXMeasurement
    settings: Settings

    @cached_property
    def _dof_df(self) -> pd.DataFrame:
        df = self.hdxm.rfu_residues
        df.columns = df.columns / 60.0
        df.loc[EXCLUDE_N_TERM[0] : EXCLUDE_N_TERM[1], :] = np.nan
        return df

    @cached_property
    def dof_df(self) -> pd.DataFrame:
        """Degree of folding dataframe"""
        return self._dof_df.reindex(self.index)

    @cached_property
    def tfold_df(self) -> pd.DataFrame:
        fracs = [0.05, 0.25, 0.5, 0.75, 0.95]
        df_tfold = pd.DataFrame(
            {
                frac: [
                    np.interp(frac, row, self._dof_df.columns)
                    for row in self._dof_df.to_numpy()
                ]
                for frac in fracs
            },
            index=self._dof_df.index,
        )
        return df_tfold.reindex(self.index)

    @cached_property
    def _dg_df(self) -> pd.DataFrame:
        dg_state = self.settings["bottom_up_state"]
        dg_pth = (
            cfg.prod_output_path("bottom_up_deltaGs")
            / f"gibbs_fit_{dg_state}"
            / "fit_result.csv"
        )
        dg_df = csv_to_dataframe(dg_pth)[dg_state]

        return dg_df

    @cached_property
    def dg_df(self) -> pd.DataFrame:
        return self._dg_df.reindex(self.index)

    @cached_property
    def index(self) -> pd.RangeIndex:
        rmin = min(self._dof_df.index.min(), self._dg_df.index.min())
        rmax = max(self._dof_df.index.max(), self._dg_df.index.max())
        return pd.RangeIndex(rmin, rmax + 1, step=1, name="r_number")

    @cached_property
    def xlim(self) -> tuple[float, float]:
        return 0, self.index.max()

    @cached_property
    def ylim(self) -> tuple[float, float]:
        return self.dof_df.columns.min(), self.settings["tmax"]

    @cached_property
    def r_edges(self) -> np.ndarray:
        r_array = self.index.to_numpy()
        r_edges = np.append(r_array - 0.5, r_array[-1] + 0.5)
        return r_edges

    @cached_property
    def flexible_parts(self) -> np.ndarray:
        flexible_parts = self.dg_df["dG"].to_numpy() < DG_CUTOFF
        return flexible_parts

    @cached_property
    def foldons_by_direct_thd(self) -> np.ndarray:
        mean = self.dof_df.mean(axis=0)
        std = 0.1  # self.dof_df.std(axis=0)
        std = self.dof_df.std(axis=0)
        foldons = self.dof_df > mean + std

        foldon_img = (foldons.sum(axis=1) > len(mean) / 2).to_numpy().astype(float)

        foldon_img[self.flexible_parts] = 0.5

        return foldon_img.reshape(1, -1)

    @cached_property
    def foldons_by_integration(self) -> np.ndarray:
        # we define foldons to be present only in the first 30 minutes of folding
        idx = list(self.dof_df.columns).index(FOLDON_CUTOFF)
        dof_early = self.dof_df.iloc[:, : idx + 1]

        # integrate over time to get total degree of folding in the 30 minute regime
        dof_time = np.trapz(dof_early.to_numpy(), x=list(dof_early.columns), axis=1)
        # take out flexible parts
        dof_time[self.flexible_parts] = np.nan
        thd = np.nanmean(dof_time) + np.nanstd(dof_time)

        # Set foldons to 1 if the integrated degree of folding is 1 std above the mean
        foldon_thd = (dof_time > thd).astype(float)

        # Set flexible parts to 0.5
        # foldon_thd[data.flexible_parts] = 0.5

        return foldon_thd.reshape(1, -1)

    def get_tfold(self, frac) -> pd.Series:
        s = self.tfold_df[frac]
        # s[self.flexible_parts] = np.nan

        return s
