#!/usr/bin/env python3

"""Create exploratory Sep-10 LGN biomarker ROI products and 10-Hz movie."""

from pathlib import Path
import json
import os

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["animation.ffmpeg_path"] = "/opt/homebrew/bin/ffmpeg"

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import t as td, ttest_rel

import sys
import config

sys.path.insert(0, str(Path(__file__).parent))

from reanalyze_hobbes_20260831_frame_indexed import mri_slice


#ROOT = Path("outputs/hobbes20260831_theta/None")
ROOT = config.OUTPUT_DIR / config.RECONSTRUCTION_TYPE
AN = ROOT / "ttl_aligned_10hz"
REC = ROOT / "reconstruction_10hz"


def tmap(x):
    n = np.sum(np.isfinite(x), 0)
    m = np.nanmean(x, 0)
    sd = np.nanstd(x, 0, ddof=1)

    return np.divide(
        m,
        sd / np.sqrt(n),
        out=np.zeros_like(m),
        where=(n >= 3) & (sd > 0),
    )


def main():
    z = np.load(
        AN / "theta_analysis_results.npz"
    )

    t = z["observed_t"]
    anat = z["anatomical_lgn_mask"].astype(bool)
    ntr = len(z["trial_run_and_number"])

    # ------------------------------------------------------------------
    # FDR selection within the anatomical LGN
    # ------------------------------------------------------------------

    p = 2 * td.sf(
        np.abs(t),
        ntr - 1,
    )

    idx = np.flatnonzero(anat)

    order = np.argsort(
        p.flat[idx]
    )

    ps = p.flat[idx][order]

    passed = np.flatnonzero(
        ps
        <=
        0.05
        * np.arange(1, len(ps) + 1)
        / len(ps)
    )

    roi = np.zeros_like(anat)

    if len(passed):
        roi.flat[
            idx[order[:passed[-1] + 1]]
        ] = True

    # Keep only positive t-values.
    roi &= t > 0

    if not np.any(roi):
        print(
            "WARNING: No positive LGN pixels survived within-ROI BH FDR. "
            "Continuing without a significant LGN ROI."
        )

        roi = np.zeros_like(
            anat,
            dtype=bool,
        )

    # ------------------------------------------------------------------
    # Load trial data
    # ------------------------------------------------------------------

    ids = z["trial_run_and_number"]
    seq = z["power_trials"].astype(np.float32)
    grid = z["time_s"]

    # These are measured Intan stimulation rising edges, aligned to the
    # first trigger of each paired block. Use their pooled median in the
    # display rather than depicting the 8.2-s train as continuous sonication.
    pulse_times = np.nanmedian(
        z["pulse_times_s"],
        axis=0,
    )

    full = np.load(
        REC / "power_doppler_10hz.npz"
    )["power_doppler"]

    run_labels = [1]
    n_runs = len(run_labels)

    # ------------------------------------------------------------------
    # Anatomical-LGN biomarker calculation
    #
    # This entire analysis is only performed if a positive FDR-selected
    # LGN ROI exists.
    # ------------------------------------------------------------------

    if np.any(roi):

        print(
            f"Using {int(roi.sum())} positive FDR-significant LGN pixels "
            "for the anatomical-LGN biomarker."
        )

        base = np.array(
            [
                np.mean(
                    np.median(
                        full[r, :80],
                        axis=0,
                    )[roi]
                )
                for r in range(n_runs)
            ]
        )

        raw = np.nanmean(
            seq[..., roi],
            axis=-1,
        )

        cbv = np.empty_like(raw)

        for rr, run in enumerate(run_labels):
            cbv[ids[:, 0] == run] = (
                100
                * (
                    raw[ids[:, 0] == run]
                    / base[rr]
                    - 1
                )
            )

        offsets = []
        adj = cbv.copy()

        for run in run_labels:
            off = float(
                np.nanmean(
                    cbv[ids[:, 0] == run]
                )
            )

            offsets.append(off)
            adj[ids[:, 0] == run] -= off

        mean = np.nanmean(
            adj,
            0,
        )

        nn = np.sum(
            np.isfinite(adj),
            0,
        )

        sem = (
            np.nanstd(
                adj,
                0,
                ddof=1,
            )
            / np.sqrt(nn)
        )

        c10 = np.nanmean(
            adj[:, :100],
            1,
        )

        i10 = np.nanmean(
            adj[:, 100:],
            1,
        )

        test = ttest_rel(
            i10,
            c10,
        )

    else:

        print(
            "No significant LGN ROI available; "
            "skipping anatomical-LGN biomarker calculations."
        )

        adj = np.full(
            (len(ids), len(grid)),
            np.nan,
            dtype=np.float32,
        )

        raw = np.full(
            (len(ids), len(grid)),
            np.nan,
            dtype=np.float32,
        )

        mean = np.full(
            len(grid),
            np.nan,
            dtype=np.float32,
        )

        sem = np.full(
            len(grid),
            np.nan,
            dtype=np.float32,
        )

        offsets = []

    # ------------------------------------------------------------------
    # Fixed ROI QC figure
    #
    # This is still created even if roi is empty. The anatomical LGN
    # contour remains visible, while the FDR-selected ROI contour simply
    # won't appear.
    # ------------------------------------------------------------------

    mri, xo, zo = mri_slice()

    xm, zm = (
        z["mri_x_mm"],
        z["mri_depth_mm"],
    )

    fig, ax = plt.subplots(
        figsize=(7, 9),
        constrained_layout=True,
    )

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

    ax.pcolormesh(
        xm,
        zm,
        t,
        cmap="coolwarm",
        vmin=-6,
        vmax=6,
        alpha=0.76,
        shading="nearest",
    )

    ax.contour(
        xm,
        zm,
        anat.astype(float),
        [0.5],
        colors="#00ffff",
        linewidths=1.3,
    )

    ax.contour(
        xm,
        zm,
        roi.astype(float),
        [0.5],
        colors="#39ff14",
        linewidths=2,
    )

    ax.set(
        xlim=(-5, 25),
        ylim=(50, -5),
        xlabel="Distance from midline (mm)",
        ylabel="MRI depth (mm)",
        title="Exploratory Sonication Biomarker within anatomical LGN",
    )

    fig.savefig(
        AN / "exploratory_biomarker_roi_qc.png",
        dpi=240,
    )

    plt.close(fig)

    # ------------------------------------------------------------------
    # Trace and bar summary
    #
    # Only create this figure if a significant LGN ROI exists.
    # Otherwise i10, c10, test, and vals do not exist and should not
    # be referenced.
    # ------------------------------------------------------------------

    if np.any(roi):

        fig, axes = plt.subplots(
            1,
            2,
            figsize=(13.2, 5.6),
            gridspec_kw={"width_ratios": [1.55, 1]},
            constrained_layout=True,
        )

        ax = axes[0]
        a = grid < 10
        b = ~a

        ax.axvspan(
            0,
            8.4,
            color="#f2c94c",
            alpha=0.12,
        )

        ax.axvspan(
            10,
            18.4,
            color="#f2c94c",
            alpha=0.12,
        )

        ax.fill_between(
            grid[a],
            mean[a] - sem[a],
            mean[a] + sem[a],
            color="blue",
            alpha=0.2,
        )

        ax.plot(
            grid[a],
            mean[a],
            color="blue",
            lw=2.5,
            label="Contralateral",
        )

        ax.fill_between(
            grid[b],
            mean[b] - sem[b],
            mean[b] + sem[b],
            color="red",
            alpha=0.2,
        )

        ax.plot(
            np.r_[grid[a][-1], grid[b]],
            np.r_[mean[a][-1], mean[b]],
            color="red",
            lw=2.5,
            label="Ipsilateral",
        )

        ax.axvline(
            10,
            color="#666",
            lw=1,
        )

        ax.axhline(
            0,
            color="#888",
            lw=0.8,
        )

        ax.set(
            xlim=(0, 20),
            xlabel="Time (s)",
            ylabel="Run-adjusted biomarker CBV change (%)",
            title="Theta-burst fUSI Neurovascular Response: Subject H",
        )

        ax.legend()

        vals = [
            i10.mean(),
            c10.mean(),
        ]

        errs = [
            i10.std(ddof=1) / np.sqrt(len(i10)),
            c10.std(ddof=1) / np.sqrt(len(c10)),
        ]

        ax = axes[1]

        ax.bar(
            [1, 2],
            vals,
            0.6,
            color=[
                (0.85, 0.2, 0.2),
                (0.2, 0.4, 0.85),
            ],
            edgecolor="black",
            linewidth=1.5,
        )

        ax.errorbar(
            [1, 2],
            vals,
            yerr=errs,
            fmt="none",
            ecolor="black",
            lw=1.5,
            capsize=6,
        )

        ax.axhline(
            0,
            color="#999",
            lw=0.8,
        )

        top = max(
            np.asarray(vals) + errs
        )

        yr = max(
            np.ptp(
                np.r_[
                    np.asarray(vals) + errs,
                    np.asarray(vals) - errs,
                ]
            ),
            1,
        )

        yb = top + 0.12 * yr
        yt = top + 0.08 * yr

        ax.plot(
            [1, 2],
            [yb, yb],
            "k-",
            lw=1.5,
        )

        ax.plot(
            [1, 1],
            [yt, yb],
            "k-",
            lw=1.5,
        )

        ax.plot(
            [2, 2],
            [yt, yb],
            "k-",
            lw=1.5,
        )

        stars = (
            "***"
            if test.pvalue < 0.001
            else "**"
            if test.pvalue < 0.01
            else "*"
            if test.pvalue < 0.05
            else "n.s."
        )

        ax.text(
            1.5,
            yb + 0.02 * yr,
            stars,
            ha="center",
            fontsize=16,
            fontweight="bold",
        )

        ax.set(
            xlim=(0.3, 2.7),
            xticks=[1, 2],
            xticklabels=[
                "Ipsilateral",
                "Contralateral",
            ],
            ylabel="Run-adjusted biomarker CBV change (%)",
            title="Mean CBV: entire 10-second period",
        )

        ax.set_ylim(
            min(
                -5,
                min(
                    np.asarray(vals) - errs
                ) - 0.1 * yr,
            ),
            yb + 0.12 * yr,
        )

        for q in axes:
            q.title.set_fontsize(16)
            q.title.set_fontweight("bold")

            q.xaxis.label.set_fontsize(14)
            q.xaxis.label.set_fontweight("bold")

            q.yaxis.label.set_fontsize(14)
            q.yaxis.label.set_fontweight("bold")

            q.tick_params(
                labelsize=12,
                width=1.5,
            )

            for lab in (
                q.get_xticklabels()
                + q.get_yticklabels()
            ):
                lab.set_fontweight("bold")

        fig.savefig(
            AN / "exploratory_biomarker_response_and_bar.png",
            dpi=300,
        )

        fig.savefig(
            AN / "exploratory_biomarker_response_and_bar.pdf"
        )

        plt.close(fig)

    else:

        print(
            "Skipping exploratory biomarker response/bar figure "
            "because no significant LGN ROI is available."
        )

    # ------------------------------------------------------------------
    # Framewise whole-image t maps
    #
    # IMPORTANT: This does NOT depend on roi.
    #
    # Therefore tm can still be calculated even when there are no
    # positive FDR-significant LGN pixels.
    # ------------------------------------------------------------------

    logs = z["log_change_trials"].astype(np.float32)

    geom = np.load(
        "outputs/hobbes20260831/lgn_roi/lgn_roi_mask.npz"
    )

    xn, zn = (
        geom["native_x_mm"][0],
        geom["native_depth_mm"][:, 0],
    )

    smm = 0.5 / 2.354820045

    sigma = (
        smm / abs(np.median(np.diff(zn))),
        smm / abs(np.median(np.diff(xn))),
    )

    tm = []

    for k in range(len(grid)):

        arr = np.stack(
            [
                gaussian_filter(
                    q,
                    sigma=sigma,
                    mode="nearest",
                )
                for q in logs[:, k]
                if np.isfinite(q[0, 0])
            ]
        )

        tm.append(
            tmap(arr)
        )

    tm = np.stack(tm)

    print(
        "Calculated whole-image time-resolved t maps."
    )

    print(
        "tm shape:",
        tm.shape,
    )

    # ------------------------------------------------------------------
    # Create whole-image t-map movie
    # ------------------------------------------------------------------

    fig = plt.figure(
        figsize=(9, 12),
        constrained_layout=True,
    )

    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[11, 2.5],
    )

    ax = fig.add_subplot(gs[0])
    tr = fig.add_subplot(gs[1])

    movie = (
        AN
        / "theta_burst_10hz_whole_image_tmap_movie.mp4"
    )

    writer = FFMpegWriter(
        fps=10,
        codec="libx264",
        bitrate=5500,
        extra_args=[
            "-pix_fmt",
            "yuv420p",
        ],
    )

    with writer.saving(
        fig,
        str(movie),
        dpi=150,
    ):

        for k, x in enumerate(grid):

            ax.clear()
            tr.clear()

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

            im = ax.pcolormesh(
                xm,
                zm,
                tm[k],
                cmap="coolwarm",
                vmin=-5,
                vmax=5,
                alpha=0.78,
                shading="nearest",
            )

            ax.set(
                xlim=(-5, 25),
                ylim=(50, -5),
                xlabel="Distance from midline (mm)",
                ylabel="MRI depth (mm)",
            )

            tr.axvspan(
                0,
                10,
                color="blue",
                alpha=0.08,
            )

            tr.axvspan(
                10,
                20,
                color="red",
                alpha=0.08,
            )

            for pulse in pulse_times:
                tr.axvline(
                    pulse,
                    color="#e5b422",
                    lw=0.65,
                    alpha=0.28,
                    zorder=0,
                )

            tr.text(
                5,
                0.74,
                "Contralateral",
                transform=tr.get_xaxis_transform(),
                ha="center",
                va="center",
                color="blue",
                weight="bold",
                fontsize=15,
            )

            tr.text(
                15,
                0.74,
                "Ipsilateral",
                transform=tr.get_xaxis_transform(),
                ha="center",
                va="center",
                color="red",
                weight="bold",
                fontsize=15,
            )

            tr.axvline(
                x,
                color="black",
                lw=3,
            )

            tr.set(
                xlim=(0, 20),
                ylim=(0, 1),
                yticks=[],
                xlabel="Time (s)",
            )

            tr.set_title(
                "Measured sonication triggers",
                weight="bold",
            )

            tr.spines[
                ["left", "right", "top"]
            ].set_visible(False)

            if k == 0:
                fig.colorbar(
                    im,
                    ax=ax,
                    label="t statistic",
                    fraction=0.046,
                    pad=0.03,
                )

            writer.grab_frame()

    plt.close(fig)

    print(
        "Saved movie:",
        movie,
    )

    # ------------------------------------------------------------------
    # Save exploratory results
    #
    # tm is saved regardless of whether an LGN ROI survived FDR.
    # ------------------------------------------------------------------

    np.savez_compressed(
        AN / "exploratory_biomarker_results.npz",
        roi_mask=roi,
        time_s=grid,
        cbv_trials=adj,
        raw_power_trials=raw,
        mean=mean,
        sem=sem,
        t_maps=tm,
        trial_run_and_number=ids,
    )

    # ------------------------------------------------------------------
    # Summary
    #
    # The summary has two possible forms depending on whether an
    # FDR-selected LGN ROI exists.
    # ------------------------------------------------------------------

    if np.any(roi):

        summary = {
            "status": "exploratory and circular",
            "definition": (
                "positive pixels surviving BH FDR q<.05 within transferred "
                "anatomical LGN for whole-10-s ipsilateral-minus-contralateral contrast"
            ),
            "pixels": int(roi.sum()),
            "retained_pairs": len(ids),
            "ipsilateral_mean_percent": float(vals[0]),
            "contralateral_mean_percent": float(vals[1]),
            "paired_t": float(test.statistic),
            "df": len(ids) - 1,
            "paired_p": float(test.pvalue),
            "run_offsets_percent": offsets,
            "movie": str(movie),
            "movie_overlay": (
                "whole-image t map only; biomarker ROI, label, and trace omitted"
            ),
        }

    else:

        summary = {
            "status": "no significant LGN ROI",
            "definition": (
                "positive pixels surviving BH FDR q<.05 within transferred "
                "anatomical LGN for whole-10-s ipsilateral-minus-contralateral contrast"
            ),
            "pixels": 0,
            "retained_pairs": len(ids),
            "ipsilateral_mean_percent": None,
            "contralateral_mean_percent": None,
            "paired_t": None,
            "df": len(ids) - 1,
            "paired_p": None,
            "run_offsets_percent": [],
            "movie": str(movie),
            "movie_overlay": (
                "whole-image time-resolved t map; "
                "no significant positive LGN biomarker ROI"
            ),
        }

    (
        AN / "exploratory_biomarker_summary.json"
    ).write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    print(
        json.dumps(
            summary,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()