#!/usr/bin/env python3
"""Authoritative BF-frame-index analysis: 10 rest + 15*(10 A + 10 B) + 10 rest."""

from pathlib import Path
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")

if os.environ.get("HOBBES_FFMPEG"):
    matplotlib.rcParams["animation.ffmpeg_path"] = os.environ["HOBBES_FFMPEG"]

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.stats import t as student_t


ROOT = Path("outputs/hobbes20260831")
OUT = ROOT / "frame_indexed_heuristic"
OUT.mkdir(parents=True, exist_ok=True)

SRC = ROOT / "lgn_roi" / "five_hz" / "power_doppler_5hz_and_censoring.npz"
FLASH = ROOT / "fusi_reconstruction" / "flash_censoring.npz"
MRI_FILE = Path("/Users/jpalmer/Desktop/FUSI/hobbes20230207.mat")

GRID = np.arange(.1, 10, .2)
PAIR_TIME = np.r_[GRID, 10 + GRID]
PULSES = np.array([0, 2.2, 4.4, 6.6, 10, 12.2, 14.4, 16.6])


def interp_block(images):
    source = np.ravel(
        np.arange(10)[:, None] + np.array([.1, .3, .5, .7])[None, :]
    )
    hi = np.searchsorted(source, GRID)
    hi = np.clip(hi, 1, len(source) - 1)
    lo = hi - 1

    a = ((GRID - source[lo]) / (source[hi] - source[lo])).astype(np.float32)

    return images[lo] * (1 - a[:, None, None]) + images[hi] * a[:, None, None]


def tstat(x):
    m = x.mean(0)
    sd = x.std(0, ddof=1)

    return m / np.maximum(sd / np.sqrt(len(x)), 1e-30)


def mri_slice():
    sys = loadmat(
        MRI_FILE,
        squeeze_me=True,
        struct_as_record=False
    )["sys"]

    px = min(
        abs(np.median(np.diff(sys.ux))),
        abs(np.median(np.diff(sys.uz)))
    ) * 1e3

    xo = np.arange(-80, 80 + px / 2, px)
    zo = np.arange(-21.5, 100 + px / 2, px)

    X, Z = np.meshgrid(xo, zo)

    th = np.deg2rad(50)
    depth = .012
    zf = float(sys.zDist) + .007 + depth * np.cos(th)
    yf = .04005 - depth * np.sin(th)

    itp = RegularGridInterpolator(
        (sys.ux, sys.uy, sys.uz),
        sys.aImg.astype(float),
        bounds_error=False,
        fill_value=np.nan
    )

    m = itp(
        np.column_stack([
            (X * 1e-3).ravel(),
            (yf - Z.ravel() * 1e-3 * np.sin(th)),
            (zf + Z.ravel() * 1e-3 * np.cos(th))
        ])
    ).reshape(X.shape)

    ok = np.isfinite(m)
    lo, hi = np.nanpercentile(m[ok], [1, 99.8])

    m = np.clip((m - lo) / (hi - lo), 0, 1)
    m[~ok] = 0

    return m, xo, zo


def main():
    z = np.load(SRC)
    power = z["power_doppler"].astype(np.float32)
    flash = np.load(FLASH)
    rz = np.abs(z["whole_field_robust_z"]).reshape(3, 320, 4).max(2)

    user = {r: [] for r in (1, 2, 3)}

    for r, f in zip(
        flash["user_flash_frame_run"],
        flash["user_flash_frames_1_based"]
    ):
        user[int(r)].append(int(f))

    auto = {
        r + 1: (np.flatnonzero(rz[r] >= 5) + 1).tolist()
        for r in range(3)
    }

    rejected = {}
    keep_parent = np.ones((3, 320), bool)
    rows = []

    for run in range(1, 4):
        trials = sorted({
            (f - 11) // 20 + 1
            for f in auto[run]
            if 11 <= f <= 310
        })

        rejected[run] = trials

        for tr in trials:
            keep_parent[
                run - 1,
                10 + (tr - 1) * 20:10 + tr * 20
            ] = False

        for f in auto[run]:
            if f <= 10 or f >= 311:
                keep_parent[run - 1, f - 1] = False

        for f in sorted(set(user[run] + auto[run])):
            tr = (f - 11) // 20 + 1 if 11 <= f <= 310 else None

            rows.append({
                "run": run,
                "parent_bf_frame_1_based": f,
                "user_flagged": f in user[run],
                "automatic_max_abs_robust_z": float(rz[run - 1, f - 1]),
                "paired_trial_removed": tr or "rest only"
            })

    left = []
    right = []
    ids = []

    for run in range(3):
        for tr in range(1, 16):
            start = 10 + (tr - 1) * 20

            if not keep_parent[run, start:start + 20].all():
                continue

            left.append(
                interp_block(
                    power[run, start * 4:(start + 10) * 4]
                )
            )

            right.append(
                interp_block(
                    power[run, (start + 10) * 4:(start + 20) * 4]
                )
            )

            ids.append((run + 1, tr))

    left = np.stack(left)
    right = np.stack(right)
    ids = np.asarray(ids)
    seq = np.concatenate([left, right], 1)

    base_runs = np.stack([
        np.median(power[r, :40], axis=0)
        for r in range(3)
    ])

    base = base_runs[ids[:, 0] - 1, None]

    logseq = (
        np.log10(np.maximum(seq, 1e-30))
        - np.log10(np.maximum(base, 1e-30))
    )

    roi = np.load(ROOT / "lgn_roi" / "lgn_roi_mask.npz")
    lgn0 = roi["lgn_mask"].astype(bool)
    xn = roi["native_x_mm"][0]
    zn = roi["native_depth_mm"][:, 0]
    xm = roi["mri_x_mm"]
    zm = roi["mri_depth_mm"]

    shift_rows = int(
        round(3 / abs(np.median(np.diff(zn))))
    )

    lgn = np.zeros_like(lgn0)
    lgn[:-shift_rows] = lgn0[shift_rows:]

    sigma_mm = .5 / 2.354820045
    sigma = (
        sigma_mm / abs(np.median(np.diff(zn))),
        sigma_mm / abs(np.median(np.diff(xn)))
    )

    smooth = np.empty_like(logseq)

    for i in range(len(logseq)):
        for k in range(100):
            smooth[i, k] = gaussian_filter(
                logseq[i, k],
                sigma=sigma,
                mode="nearest"
            )

    tmaps = np.stack([
        tstat(smooth[:, k])
        for k in range(100)
    ])

    base_lgn = np.mean(
        base_runs[:, lgn],
        axis=1
    )[ids[:, 0] - 1, None]

    lgn_trials = 100 * (
        np.mean(seq[..., lgn], axis=-1) / base_lgn - 1
    )

    lm = lgn_trials.mean(0)
    lse = lgn_trials.std(0, ddof=1) / np.sqrt(len(lgn_trials))

    # Full continuous ROI trace.
    fig, ax = plt.subplots(
        figsize=(13, 5),
        constrained_layout=True
    )

    ax.plot(PAIR_TIME, lm, color="black", lw=2)

    ax.fill_between(
        PAIR_TIME,
        lm - lse,
        lm + lse,
        color="black",
        alpha=.18
    )

    for p0 in PULSES:
        ax.axvspan(
            p0,
            p0 + .1,
            color="#f2c94c",
            alpha=.7
        )

    ax.axhline(0, color="#777", lw=.8)

    ax.set(
        xlim=(0, 20),
        xlabel="Time through complete 20-frame pair (s)",
        ylabel="LGN CBV change from initial 10-frame baseline (%)",
        title=f"Frame-indexed continuous LGN response (n={len(ids)} retained pairs; ROI shifted 3 mm superiorly)"
    )

    ax.grid(alpha=.2)

    fig.savefig(
        OUT / "lgn_full_20s_response_frame_indexed.png",
        dpi=220
    )

    plt.close(fig)

    # Smoothed instantaneous maps, retaining the earlier 0.5--0.9 s pulse window.
    phase = np.arange(.1, 2.2, .2)
    lp = []
    rp = []

    for onset in (0, 2.2, 4.4, 6.6):
        ix = np.array([
            np.argmin(abs(GRID - (onset + p)))
            for p in phase
        ])

        lp.append(logseq[:, :50][:, ix])
        rp.append(logseq[:, 50:][:, ix])

    lp = np.stack(lp, 1).mean(1)
    rp = np.stack(rp, 1).mean(1)

    wi = (phase >= .5) & (phase <= .9)

    contrasts = {
        "Second − first": (rp[:, wi] - lp[:, wi]).mean(1),
        "Second − baseline": rp[:, wi].mean(1),
        "First − baseline": lp[:, wi].mean(1)
    }

    smooth_con = {
        name: np.stack([
            gaussian_filter(a, sigma=sigma, mode="nearest")
            for a in arr
        ])
        for name, arr in contrasts.items()
    }

    mri, xo, zo = mri_slice()

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16, 10),
        constrained_layout=True
    )

    for ax, (name, arr) in zip(axes, smooth_con.items()):
        tt = tstat(arr)

        ax.imshow(
            mri,
            cmap="gray",
            origin="upper",
            extent=[xo[0], xo[-1], zo[-1], zo[0]],
            aspect="equal"
        )

        im = ax.pcolormesh(
            xm,
            zm,
            tt,
            cmap="coolwarm",
            vmin=-5,
            vmax=5,
            alpha=.78,
            shading="nearest"
        )

        pos = tt > 4

        if np.any(pos):
            ax.contour(
                xm,
                zm,
                pos.astype(float),
                levels=[.5],
                colors="yellow",
                linewidths=1.2
            )

        ax.set(
            xlim=(-5, 25),
            ylim=(50, -5),
            title=name + "\n0.5–0.9 s after each pulse",
            xlabel="Distance from midline (mm)",
            ylabel="MRI depth (mm)"
        )

        plt.colorbar(
            im,
            ax=ax,
            label="t statistic",
            fraction=.046
        )

    fig.suptitle(
        "Frame-indexed instantaneous contrasts — 0.5 mm FWHM (t > 4 outlined)",
        fontsize=17
    )

    fig.savefig(
        OUT / "instantaneous_smoothed_t_maps_frame_indexed.png",
        dpi=220
    )

    plt.close(fig)

    # Neutral 20-second movie with shifted LGN contour and synchronized trace.
    fig = plt.figure(
        figsize=(9, 12),
        constrained_layout=True
    )

    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[11, 2.5]
    )

    ax = fig.add_subplot(gs[0])
    trace = fig.add_subplot(gs[1])

    out = OUT / "full_20s_movie_frame_indexed_with_lgn_roi.mp4"

    writer = FFMpegWriter(
        fps=5,
        codec="libx264",
        bitrate=5500,
        extra_args=["-pix_fmt", "yuv420p"]
    )

    with writer.saving(fig, str(out), dpi=150):
        for k, t in enumerate(PAIR_TIME):
            ax.clear()
            trace.clear()

            ax.imshow(
                mri,
                cmap="gray",
                origin="upper",
                extent=[xo[0], xo[-1], zo[-1], zo[0]],
                aspect="equal"
            )

            im = ax.pcolormesh(
                xm,
                zm,
                tmaps[k],
                cmap="coolwarm",
                vmin=-5,
                vmax=5,
                alpha=.78,
                shading="nearest"
            )

            pos = tmaps[k] > 4

            if np.any(pos):
                ax.contour(
                    xm,
                    zm,
                    pos.astype(float),
                    levels=[.5],
                    colors="yellow",
                    linewidths=1.2
                )

            ax.contour(
                xm,
                zm,
                lgn.astype(float),
                levels=[.5],
                colors="#00ffff",
                linewidths=2
            )

            iy, ix = np.argwhere(lgn).mean(0)
            iy = int(round(iy))
            ix = int(round(ix))

            ax.text(
                float(xm[iy, ix]) + .8,
                float(zm[iy, ix]),
                "LGN ROI",
                color="#00ffff",
                weight="bold",
                fontsize=10
            )

            active = np.any(
                (t >= PULSES) &
                (t <= PULSES + .11)
            )

            ax.set(
                xlim=(-5, 25),
                ylim=(50, -5),
                xlabel="Distance from midline (mm)",
                ylabel="MRI depth (mm)",
                title=f"Full 20 s sequence | t={t:.1f} s | Sonication: {'ON' if active else 'off'}"
            )

            ax.text(
                .02,
                .98,
                "Frame-indexed; first 10 BF frames baseline; 0.5 mm FWHM; t > 4 outlined",
                transform=ax.transAxes,
                ha="left",
                va="top",
                color="white",
                fontsize=8.5,
                bbox=dict(
                    facecolor="black",
                    alpha=.6,
                    pad=4,
                    edgecolor="none"
                )
            )

            trace.plot(
                PAIR_TIME,
                lm,
                color="black",
                lw=2
            )

            trace.fill_between(
                PAIR_TIME,
                lm - lse,
                lm + lse,
                color="black",
                alpha=.18
            )

            for p0 in PULSES:
                trace.axvspan(
                    p0,
                    p0 + .1,
                    color="#f2c94c",
                    alpha=.75
                )

            trace.axvline(
                t,
                color="#d62728",
                lw=2.5
            )

            trace.axhline(
                0,
                color="#777",
                lw=.7
            )

            trace.set(
                xlim=(0, 20),
                xlabel="Time through complete sequence (s)",
                ylabel="LGN CBV change (%)",
                title=f"Pooled LGN response (mean ± SEM; n={len(ids)} paired trials)"
            )

            trace.grid(alpha=.2)

            if k == 0:
                fig.colorbar(
                    im,
                    ax=ax,
                    label="t statistic",
                    fraction=.046,
                    pad=.03
                )

            writer.grab_frame()

    plt.close(fig)

    with (OUT / "outlier_flashes_frame_indexed.csv").open(
        "w",
        newline=""
    ) as f:
        w = csv.DictWriter(
            f,
            fieldnames=rows[0]
        )
        w.writeheader()
        w.writerows(rows)

    summary = {
        "labeling": "BF frames 1-10 rest; 11-310 fifteen repetitions of first 10 then second 10; 311-320 rest",
        "raw_bf_frames_per_run": 320,
        "reconstructed_images_per_parent_bf": 4,
        "baseline_parent_frames_1_based": [1, 10],
        "retained_pairs": len(ids),
        "retained_by_run": {
            str(r): int(np.sum(ids[:, 0] == r))
            for r in (1, 2, 3)
        },
        "user_flagged_frames": user,
        "automatic_frames_abs_robust_z_ge_5": auto,
        "removed_paired_trials": {
            str(k): v
            for k, v in rejected.items()
        },
        "run3_final_rest_contamination": "Frames 311-320 excluded from any rest summaries; initial baseline only used",
        "outputs": {
            "movie": out.name,
            "lgn_trace": "lgn_full_20s_response_frame_indexed.png",
            "instantaneous_maps": "instantaneous_smoothed_t_maps_frame_indexed.png"
        }
    }

    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    np.savez_compressed(
        OUT / "frame_indexed_analysis_results.npz",
        grid_10s=GRID,
        time_20s=PAIR_TIME,
        left_power=left,
        right_power=right,
        trial_run_and_number=ids,
        keep_parent_bf=keep_parent,
        lgn_mask_shifted=lgn,
        lgn_cbv_trials_percent=lgn_trials,
        lgn_cbv_mean_percent=lm,
        lgn_cbv_sem_percent=lse,
        t_maps_20s=tmaps
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()