# %%
import ultraplot as uplt
import numpy as np
from hal.repro import reproduce

# %%

OUTPUT = reproduce(globals())

smfret_rates = {"WT": 8.2e-3, "Y283D": 1.9e-4}
hdx_rates = {"WT": 9e-4, "Y283D": 1.7e-4}
hdx_beta = {"WT": 0.78, "Y283D": 0.71}


# %%

time = np.logspace(0, 5, 1000)
time


# %%

fig, ax = uplt.subplots()


k = smfret_rates["WT"]
ax.plot(time, 1 - np.exp(-k * time), label="smFRET WT", color="blue")

k = smfret_rates["Y283D"]
ax.plot(time, 1 - np.exp(-k * time), label="smFRET Y283D", color="orange")

k, beta = hdx_rates["WT"], hdx_beta["WT"]
ax.plot(
    time,
    1 - np.exp(-((k * time) ** beta)),
    label="HDX WT",
    color="blue",
    linestyle="--",
)

k, beta = hdx_rates["Y283D"], hdx_beta["Y283D"]
ax.plot(
    time,
    1 - np.exp(-((k * time) ** beta)),
    label="HDX Y283D",
    color="orange",
    linestyle="--",
)

ax.format(
    xscale="log",
    xlabel="Time (s)",
    ylabel="Fraction folded",
    title="Comparison of smFRET and HDX rates",
    xlim=(1, 1e5),
    ylim=(0, 1),
    xformatter="log",
)

ax.legend(ncols=1, loc="r")
fig.savefig(OUTPUT / "smfret_vs_hdx_rates.png", dpi=300, transparent=True)

# %%
