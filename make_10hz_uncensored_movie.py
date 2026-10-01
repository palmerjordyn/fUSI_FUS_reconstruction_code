#!/usr/bin/env python3

"""All-trial 10-Hz t-map movie; intentionally performs no artifact censoring."""

from pathlib import Path
import os

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["animation.ffmpeg_path"] = os.environ.get(
    "HOBBES_FFMPEG",
    "ffmpeg",
)

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
from scipy.ndimage import gaussian_filter

import sys
import config 

sys.path.insert(0, str(Path(__file__).parent))

from reanalyze_hobbes_20260831_frame_indexed import mri_slice


ROOT = config.OUTPUT_DIR / config.RECONSTRUCTION_TYPE
print(ROOT)

OUT = ROOT / "ttl_aligned_10hz"

GRID = np.arange(0.05, 20, 0.1)
CENTERS = (np.arange(8) + 0.5) / 10
RUNS = np.array([1, 3, 4])

CENSOR_Z = float(
    os.environ.get(
        "HOBBES_CENSOR_Z",
        "nan",
    )
)


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


def nearest(samples, times):
    out = np.full(
        (len(GRID),) + samples.shape[1:],
        np.nan,
        np.float32,
    )

    for k, g in enumerate(GRID):
        j = np.argmin(abs(times - g))

        if abs(times[j] - g) <= 0.06:
            out[k] = samples[j]

    return out


def main():
    p = np.load(
        ROOT / "reconstruction_10hz/power_doppler_10hz.npz"
    )

    power = p["power_doppler"].astype(np.float32)

    timing = np.load(
        config.OUTPUT_DIR / "design/selected_run_timing.npz"
    )

    fu = timing["fusi_ttl_s"]
    blocks = timing["block_pulse_ttl_s"]

    geom = np.load(
        "outputs/hobbes20260831/lgn_roi/lgn_roi_mask.npz"
    )

    xm, zm = (
        geom["mri_x_mm"],
        geom["mri_depth_mm"],
    )

    roi = np.load(
        "outputs/hobbes20260831/ttl_aligned_corrected/"
        "ttl_aligned_analysis_results.npz"
    )["lgn_mask_shifted"].astype(bool)

    base = np.stack(
        [
            np.median(
                power[r, :80],
                axis=0,
            )
            for r in range(3)
        ]
    )

    seq = []
    logs = []
    ids = []
    pulse_times = []
    removed = {}

    for r, run in enumerate(RUNS):
        child = (
            fu[r, :, None] + CENTERS
        ).ravel()

        if np.isfinite(CENSOR_Z):
            parent_z = (
                np.abs(
                    p["whole_field_robust_z"][r]
                )
                .reshape(320, 8)
                .max(1)
            )

            bad = sorted(
                {
                    int((f - 11) // 20 + 1)
                    for f in np.flatnonzero(
                        parent_z >= CENSOR_Z
                    ) + 1
                    if 11 <= f <= 310
                }
            )
        else:
            bad = []

        removed[str(int(run))] = bad

        for trial in range(1, 16):
            if trial in bad:
                continue

            origin = blocks[
                r,
                2 * (trial - 1),
                0,
            ]

            rel = child - origin

            valid = (
                (rel >= -0.1)
                & (rel <= 20.1)
            )

            s = nearest(
                power[r, valid],
                rel[valid],
            )

            seq.append(s)

            logs.append(
                np.log10(
                    np.maximum(s, 1e-30)
                )
                - np.log10(
                    np.maximum(
                        base[r],
                        1e-30,
                    )
                )
            )

            ids.append(
                (
                    run,
                    trial,
                )
            )

            pulse_times.append(
                blocks[
                    r,
                    2 * (trial - 1) : 2 * trial,
                ].reshape(-1)
                - origin
            )

    seq = np.stack(seq)
    logs = np.stack(logs)
    ids = np.asarray(ids)

    pulses = np.nanmedian(
        np.stack(pulse_times),
        axis=0,
    )

    raw = np.nanmean(
        seq[..., roi],
        axis=-1,
    )

    cbv = np.empty_like(raw)

    for rr, run in enumerate(RUNS):
        hit = ids[:, 0] == run

        cbv[hit] = (
            100
            * (
                raw[hit]
                / np.mean(base[rr][roi])
                - 1
            )
        )

        cbv[hit] -= np.nanmean(cbv[hit])

    mean = np.nanmean(
        cbv,
        0,
    )

    n = np.sum(
        np.isfinite(cbv),
        0,
    )

    sem = (
        np.nanstd(
            cbv,
            0,
            ddof=1,
        )
        / np.sqrt(n)
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

    maps = []

    for k in range(len(GRID)):
        a = np.stack(
            [
                gaussian_filter(
                    v,
                    sigma=sigma,
                    mode="nearest",
                )
                for v in logs[:, k]
                if np.isfinite(v[0, 0])
            ]
        )

        maps.append(
            tmap(a)
        )

    maps = np.stack(maps)

    mri, xo, zo = mri_slice()

    tag = (
        f"z{CENSOR_Z:g}_cutoff"
        if np.isfinite(CENSOR_Z)
        else "ALL_45_TRIALS_uncensored"
    )

    movie = OUT / (
        f"theta_burst_10hz_anatomical_lgn_{tag}.mp4"
    )

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
        for k, t in enumerate(GRID):
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
                maps[k],
                cmap="coolwarm",
                vmin=-5,
                vmax=5,
                alpha=0.78,
                shading="nearest",
            )

            ax.contour(
                xm,
                zm,
                roi,
                [0.5],
                colors="#00ffff",
                linewidths=2.2,
            )

            iy, ix = (
                np.argwhere(roi)
                .mean(0)
                .round()
                .astype(int)
            )

            ax.text(
                float(xm[iy, ix]) + 2,
                float(zm[iy, ix]),
                "Anatomical\nLGN ROI",
                color="white",
                weight="bold",
                fontsize=17,
                va="center",
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

            for pulse in pulses:
                tr.axvline(
                    pulse,
                    color="#e5b422",
                    lw=0.65,
                    alpha=0.28,
                    zorder=0,
                )

            tr.text(
                5,
                0.96,
                "Contralateral",
                transform=tr.get_xaxis_transform(),
                ha="center",
                va="top",
                color="blue",
                weight="bold",
            )

            tr.text(
                15,
                0.96,
                "Ipsilateral",
                transform=tr.get_xaxis_transform(),
                ha="center",
                va="top",
                color="red",
                weight="bold",
            )

            tr.plot(
                GRID,
                mean,
                color="black",
                lw=2,
            )

            tr.fill_between(
                GRID,
                mean - sem,
                mean + sem,
                color="black",
                alpha=0.18,
            )

            tr.axvline(
                t,
                color="#d62728",
                lw=2.5,
            )

            tr.axhline(
                0,
                color="#777",
                lw=0.7,
            )

            title = (
                f"Whole-field artifact cutoff |z| ≥ {CENSOR_Z:g}; "
                f"n={len(seq)} paired trials"
                if np.isfinite(CENSOR_Z)
                else f"No artifact censoring; n={len(seq)} paired trials"
            )

            tr.set(
                xlim=(0, 20),
                xlabel="Time (s)",
                ylabel="Anatomical LGN CBV change (%)",
                title=title,
            )

            tr.grid(alpha=0.2)

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

    result = OUT / (
        f"theta_10hz_{tag}_results.npz"
    )

    np.savez_compressed(
        result,
        time_s=GRID,
        trial_run_and_number=ids,
        anatomical_lgn_cbv_trials=cbv,
        mean=mean,
        sem=sem,
        t_maps=maps,
        censor_z=CENSOR_Z,
        removed_pairs=np.asarray(
            [
                removed[str(r)]
                for r in RUNS
            ],
            dtype=object,
        ),
    )

    print(
        {
            "movie": str(movie),
            "retained_pairs": len(seq),
            "removed_pairs": removed,
        }
    )


if __name__ == "__main__":
    main()
