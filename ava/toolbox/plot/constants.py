from typing import Any

import ultraplot as uplt

SCATTER_KWARGS: dict[str, Any] = {"s": 20}
FIT_CONTOUR_KWARGS: dict[str, Any] = dict(color="magenta", lw=0.4)
CONTOURF_KWARGS_LINES: dict[str, Any] = dict(levels=15, lw=0.4, ec="w", cmap="viridis")
CONTOURF_KWARGS: dict[str, Any] = dict(levels=15, cmap="viridis")

HIST_KWARGS: dict[str, Any] = dict(color="#4d4d4d", bins="fd")
SAVE_KWARGS: dict[str, Any] = dict(dpi=300, transparent=True)
Y_FMT_KWARGS: dict[str, Any] = dict(
    yformatter=uplt.Formatter("{x:.1f}"), yticks=uplt.arange(0, 1, 0.2)
)
X_FMT_KWARGS: dict[str, Any] = dict(
    xformatter=uplt.Formatter("{x:.1f}"), xticks=uplt.arange(0, 1, 0.2)
)

# figure 2 width in mm per second of folding time
F2_WIDTH_PER_SECOND = 60 / 1800.0

STATE_COLORS = {
    "I1": "#FF8A24",
    "Native": "#006600",
    "N": "#006600",
    "Np": "#006633",
    "I2": "#0AC2B5",
    "RC": "#ADADAD",
}

SIGMA = 10

BW_TIME = 100
BW_E = 0.05
BW_TE = (BW_TIME, BW_E)
