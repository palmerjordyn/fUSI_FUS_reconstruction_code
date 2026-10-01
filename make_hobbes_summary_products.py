#!/usr/bin/env python3

"""Publication-style overview and ROI traces for the Sep-10 theta-burst data."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import config

from reanalyze_hobbes_20260831_frame_indexed import mri_slice


ROOT = config.OUTPUT_DIR / config.RECONSTRUCTION_TYPE
OUT = ROOT / "ttl_aligned_10hz"


def mean_sem(x):
    return (
        np.nanmean(x, axis=0),
        np.nanstd(x, axis=0, ddof=1)
        / np.sqrt(np.sum(np.isfinite(x), axis=0)),
    )


def style(ax):
    ax.tick_params(
        labelsize=11,
        width=1.2,
    )

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontsize(12)


def plot_roi_traces(
    time,
    raw,
    cbv,
    ids,
    prefix,
    label,
):
    # Remove multiplicative run gain from raw power so runs can be pooled without
    # changing its meaning into percent CBV.
    raw_norm = raw.copy().astype(float)

    for run in np.unique(ids[:, 0]):
        q = ids[:, 0] == run
        raw_norm[q] /= np.nanmean(raw_norm[q])

    traces = (
        (
            raw_norm,
            "Run-normalized power Doppler (a.u.)",
            "power_doppler",
        ),
        (
            cbv,
            f"Run-adjusted {label} CBV change (%)",
            "cbv",
        ),
    )

    for values, ylabel, suffix in traces:
        mean, sem = mean_sem(values)

        fig, ax = plt.subplots(
            figsize=(10.5, 4.6),
            constrained_layout=True,
        )

        left = time < 10
        right = ~left

        ax.axvspan(
            0,
            10,
            color="#2f61d5",
            alpha=0.055,
        )

        ax.axvspan(
            10,
            20,
            color="#f32323",
            alpha=0.055,
        )

        ax.fill_between(
            time[left],
            mean[left] - sem[left],
            mean[left] + sem[left],
            color="#2f61d5",
            alpha=0.2,
        )

        ax.plot(
            time[left],
            mean[left],
            color="#2f61d5",
            lw=2.3,
            label="Contralateral",
        )

        ax.fill_between(
            time[right],
            mean[right] - sem[right],
            mean[right] + sem[right],
            color="#f32323",
            alpha=0.2,
        )

        ax.plot(
            np.r_[time[left][-1], time[right]],
            np.r_[mean[left][-1], mean[right]],
            color="#f32323",
            lw=2.3,
            label="Ipsilateral",
        )

        for start in (0, 10):
            for pulse in np.arange(
                start,
                start + 8.21,
                0.2,
            ):
                ax.axvline(
                    pulse,
                    color="#e5b422",
                    lw=0.55,
                    alpha=0.16,
                    zorder=0,
                )

        ax.axvline(
            10,
            color="#555",
            lw=1,
        )

        if suffix == "cbv":
            ax.axhline(
                0,
                color="#888",
                lw=0.8,
            )

        ax.set(
            xlim=(0, 20),
            xlabel="Time (s)",
            ylabel=ylabel,
        )

        ax.legend(
            frameon=False,
            loc="best",
        )

        style(ax)

        fig.savefig(
            OUT / f"{prefix}_{suffix}_trace.png",
            dpi=300,
        )

        fig.savefig(
            OUT / f"{prefix}_{suffix}_trace.pdf"
        )

        plt.close(fig)


def main():
    z = np.load(
        OUT / "theta_analysis_results.npz"
    )

    b = np.load(
        ROOT / "reconstruction_10hz/temporal_median_bmode.npz"
    )["bmode_db"]

    e = np.load(
        OUT / "exploratory_biomarker_results.npz"
    )

    time = z["time_s"]
    ids = z["trial_run_and_number"]

    plot_roi_traces(
        time,
        z["anatomical_lgn_raw_power_trials"],
        z["anatomical_lgn_cbv_trials"],
        ids,
        "anatomical_lgn",
        "anatomical LGN",
    )

    plot_roi_traces(
        time,
        e["raw_power_trials"],
        e["cbv_trials"],
        ids,
        "exploratory_biomarker",
        "biomarker",
    )

    mri, xo, zo = mri_slice()

    xm, zm = (
        z["mri_x_mm"],
        z["mri_depth_mm"],
    )

    anatomical = z[
        "anatomical_lgn_mask"
    ].astype(bool)

    biomarker = e[
        "roi_mask"
    ].astype(bool)

    t = z["observed_t"]

    view_x, view_z = (
        (-65, 75),
        (-8, 82),
    )

    def anatomy(ax, title):
        ax.imshow(
            mri,
            cmap="gray",
            origin="upper",
            extent=[
                xo[0],
                xo[-1],
                zo[-1],
                zo[0],
            ],
            aspect="equal",
        )

        ax.set(
            xlim=view_x,
            ylim=(
                view_z[1],
                view_z[0],
            ),
            xlabel="Distance from midline (mm)",
            ylabel="MRI depth (mm)",
            title=title,
        )

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(19, 6.7),
        constrained_layout=True,
    )

    anatomy(
        axes[0],
        "A. Monkey MRI",
    )

    axes[0].contour(
        xm,
        zm,
        anatomical,
        [0.5],
        colors="#00ffff",
        linewidths=2,
    )

    anatomy(
        axes[1],
        "B. Registered Functional Ultrasound",
    )

    bn = np.clip(
        (b + 60) / 60,
        0,
        1,
    )

    axes[1].pcolormesh(
        xm,
        zm,
        bn,
        cmap="viridis",
        vmin=0,
        vmax=1,
        alpha=np.where(
            bn > 0.03,
            0.68,
            0,
        ),
        shading="nearest",
    )

    axes[1].contour(
        xm,
        zm,
        anatomical,
        [0.5],
        colors="#00ffff",
        linewidths=2,
    )

    anatomy(
        axes[2],
        "C. Changes with Stimulation",
    )

    axes[2].pcolormesh(
        xm,
        zm,
        bn,
        cmap="gray",
        vmin=0,
        vmax=1,
        alpha=np.where(
            bn > 0.03,
            0.20,
            0,
        ),
        shading="nearest",
    )

    axes[2].pcolormesh(
        xm,
        zm,
        t,
        cmap="coolwarm",
        vmin=-6,
        vmax=6,
        alpha=0.78,
        shading="nearest",
    )

    axes[2].contour(
        xm,
        zm,
        anatomical,
        [0.5],
        colors="#00ffff",
        linewidths=2,
    )

    axes[2].contour(
        xm,
        zm,
        biomarker,
        [0.5],
        colors="#39ff14",
        linewidths=2.3,
    )

    axes[2].text(
        15.5,
        30,
        "Exploratory\nsonication biomarker",
        color="white",
        weight="bold",
        fontsize=11,
        va="center",
    )

    for ax in axes:
        ax.title.set_fontsize(15)

    fig.savefig(
        OUT / "three_panel_mri_registered_fusi_contrast.png",
        dpi=250,
    )

    fig.savefig(
        OUT / "three_panel_mri_registered_fusi_contrast.pdf"
    )

    plt.close(fig)


if __name__ == "__main__":
    main()