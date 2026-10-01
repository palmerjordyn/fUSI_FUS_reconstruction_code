#!/usr/bin/env python3
"""Bar plot paired whole-period CBV in the transferred anatomical LGN."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import config
import numpy as np
from scipy.stats import ttest_rel


OUT = config.OUTPUT_DIR / config.RECONSTRUCTION_TYPE / "ttl_aligned_10hz"

z = np.load(OUT / "theta_analysis_results.npz")
x = z["anatomical_lgn_cbv_trials"]

contra = np.nanmean(x[:, :100], axis=1)
ipsi = np.nanmean(x[:, 100:], axis=1)

test = ttest_rel(ipsi, contra)

means = np.array([ipsi.mean(), contra.mean()])
sems = np.array([ipsi.std(ddof=1), contra.std(ddof=1)]) / np.sqrt(len(ipsi))


fig, ax = plt.subplots(
    figsize=(6.4, 5.7),
    constrained_layout=True
)

ax.bar(
    [1, 2],
    means,
    width=.62,
    color=["#f32323", "#2f61d5"],
    edgecolor="black",
    linewidth=1.5
)

ax.errorbar(
    [1, 2],
    means,
    yerr=sems,
    fmt="none",
    ecolor="black",
    lw=1.5,
    capsize=6
)

ax.axhline(0, color="#999", lw=.8)

low = float(np.min(means - sems))
high = float(np.max(means + sems))
span = max(high - low, 1)

y = high + .20 * span
tick = y - .06 * span

ax.plot(
    [1, 1, 2, 2],
    [tick, y, y, tick],
    color="black",
    lw=1.5
)

stars = (
    "***" if test.pvalue < .001
    else "**" if test.pvalue < .01
    else "*" if test.pvalue < .05
    else "n.s."
)

ax.text(
    1.5,
    y + .025 * span,
    stars,
    ha="center",
    va="bottom",
    fontsize=15,
    fontweight="bold"
)

ax.set(
    xlim=(.35, 2.65),
    xticks=[1, 2],
    xticklabels=["Ipsilateral", "Contralateral"],
    ylabel="Run-adjusted anatomical LGN CBV change (%)",
    title="Entire 10-second stimulation period"
)

ax.set_ylim(
    min(low - .18 * span, -2),
    y + .16 * span
)

ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=12, width=1.3)

ax.title.set_fontsize(15)
ax.title.set_fontweight("bold")
ax.yaxis.label.set_fontsize(12)

fig.savefig(
    OUT / "anatomical_lgn_whole_10s_bar.png",
    dpi=300
)

fig.savefig(
    OUT / "anatomical_lgn_whole_10s_bar.pdf"
)

print({
    "n": len(ipsi),
    "ipsilateral_mean": float(means[0]),
    "contralateral_mean": float(means[1]),
    "t": float(test.statistic),
    "df": len(ipsi) - 1,
    "p": float(test.pvalue)
})