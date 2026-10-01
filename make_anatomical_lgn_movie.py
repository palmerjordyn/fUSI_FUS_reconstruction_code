#!/usr/bin/env python3

"""10-Hz whole-image response movie with the transferred anatomical LGN trace."""

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

import sys
import config

sys.path.insert(0, str(Path(__file__).parent))

from reanalyze_hobbes_20260831_frame_indexed import mri_slice


#OUT = Path("outputs/hobbes20260831_theta/Butter2/ttl_aligned_10hz")
OUT = config.OUTPUT_DIR / config.RECONSTRUCTION_TYPE / "ttl_aligned_10hz"

def main():
    z = np.load(
        OUT / "theta_analysis_results.npz"
    )

    e = np.load(
        OUT / "exploratory_biomarker_results.npz"
    )

    time = z["time_s"]
    roi = z["anatomical_lgn_mask"].astype(bool)

    cbv = z["anatomical_lgn_cbv_trials"]

    mean = np.nanmean(
        cbv,
        axis=0,
    )

    n = np.sum(
        np.isfinite(cbv),
        axis=0,
    )

    sem = (
        np.nanstd(
            cbv,
            axis=0,
            ddof=1,
        )
        / np.sqrt(n)
    )

    tm = e["t_maps"]
    #tm = z["observed_t"]

    xm, zm = (
        z["mri_x_mm"],
        z["mri_depth_mm"],
    )

    pulses = np.nanmedian(
        z["pulse_times_s"],
        axis=0,
    )

    mri, xo, zo = mri_slice()

    movie = (
        OUT
        / "theta_burst_10hz_anatomical_lgn_response_movie.mp4"
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
        for k, t in enumerate(time):
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

            ax.contour(
                xm,
                zm,
                roi.astype(float),
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
                    lw=2,
                    alpha=0.9,
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
                time,
                mean,
                color="black",
                lw=2,
            )

            tr.fill_between(
                time,
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

            tr.set(
                xlim=(0, 20),
                xlabel="Time (s)",
                ylabel="Anatomical LGN CBV change (%)",
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

    print(movie)


if __name__ == "__main__":
    main()