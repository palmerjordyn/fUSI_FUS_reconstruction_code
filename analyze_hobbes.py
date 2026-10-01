#!/usr/bin/env python3

"""TTL-aligned 10-Hz analysis of Sep-10 theta-burst runs 1, 3, and 4."""

from pathlib import Path
import csv
import json

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, label
from scipy.stats import t as tdist, ttest_rel

import sys
import config

sys.path.insert(0, str(Path(__file__).parent))

from reanalyze_hobbes_20260831_frame_indexed import mri_slice

ROOT = config.OUTPUT_DIR / config.RECONSTRUCTION_TYPE

OUT = ROOT / "ttl_aligned_10hz"
OUT.mkdir(parents=True, exist_ok=True)

PD = ROOT / "reconstruction_10hz/power_doppler_10hz.npz"
TIMING = config.OUTPUT_DIR / "design/selected_run_timing.npz"

GRID = np.arange(0.05, 20, 0.1)
CENTERS = (np.arange(8) + 0.5) / 10
#RUN_LABELS = np.array([1, 2, 3])
RUN_LABELS = config.RHS_IDENTIFIER
n_runs = len(RUN_LABELS)
n_PARENT = max(config.BF)

CLUSTER_P = 0.01
N_PERM = 10000
ALPHA = 0.05


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


def tmap(x):
    n = np.sum(np.isfinite(x), axis=0)
    m = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0, ddof=1)

    return np.divide(
        m,
        sd / np.sqrt(n),
        out=np.zeros_like(m),
        where=(n >= 3) & (sd > 0),
    )


def cluster_test(d):
    obs = tmap(d)

    thr = float(
        tdist.ppf(
            1 - CLUSTER_P / 2,
            len(d) - 1,
        )
    )

    structure = np.ones((3, 3), int)

    def get(tt):
        ans = []

        for sign in (1, -1):
            lab, n = label(
                sign * tt > thr,
                structure,
            )

            for j in range(1, n + 1):
                q = lab == j
                ans.append(
                    (
                        sign,
                        q,
                        float(np.sum(np.abs(tt[q]))),
                    )
                )

        return ans

    observed = get(obs)

    rng = np.random.default_rng(20260910)
    null = np.zeros(N_PERM)

    for p in range(N_PERM):
        signs = rng.choice(
            (-1, 1),
            size=(len(d), 1, 1),
        )

        cc = get(
            d * signs if False else tmap(d * signs)
        )

        null[p] = max(
            (x[2] for x in cc),
            default=0,
        )

    sig = np.zeros(obs.shape, np.int16)
    rows = []
    rank = 0

    for sign, q, mass in sorted(
        observed,
        key=lambda x: x[2],
        reverse=True,
    ):
        pv = (
            1 + np.sum(null >= mass)
        ) / (N_PERM + 1)

        if pv < ALPHA:
            rank += 1
            sig[q] = rank

            iy, ix = np.unravel_index(
                np.argmax(
                    np.where(
                        q,
                        np.abs(obs),
                        -np.inf,
                    )
                ),
                obs.shape,
            )

            rows.append(
                (
                    rank,
                    sign,
                    pv,
                    mass,
                    int(q.sum()),
                    float(obs[iy, ix]),
                    iy,
                    ix,
                )
            )

    return obs, thr, null, sig, rows


def main():
    p = np.load(PD)

    power = p["power_doppler"].astype(np.float32)

    rz = (
        np.abs(p["whole_field_robust_z"])
        .reshape(n_runs, n_PARENT, 8)
        .max(2)
    )

    timing = np.load(TIMING)

    fu = timing["fusi_ttl_s"]
    blocks = timing["block_pulse_ttl_s"]

    reject = {}
    keep = np.ones((n_runs, n_PARENT), bool)
    artifact_frames = {}

    for r in range(n_runs):
        frames = np.flatnonzero(rz[r] >= 5) + 1

        artifact_frames[
            RUN_LABELS[r]
        ] = frames.tolist()

        trials = sorted(
            {
                int((f - 11) // 20 + 1)
                for f in frames
                if 11 <= f <= 310
            }
        )

        reject[r] = trials

        for tr in trials:
            keep[
                r,
                10 + (tr - 1) * 20 : 10 + tr * 20,
            ] = False

    base = np.stack(
        [
            np.median(
                power[r, :80],
                axis=0,
            )
            for r in range(n_runs)
        ]
    )

    seqs = []
    logs = []
    ids = []
    pulse_rel = []

    for r in range(n_runs):
        child = (
            fu[r, :, None] + CENTERS
        ).ravel()

        parent = np.repeat(
            np.arange(n_PARENT),
            8,
        )

        for tr in range(1, 16):
            if tr in reject[r]:
                continue
                
            origin = blocks[
                r,
                2 * (tr - 1),
                0,
            ]

            rel = child - origin

            valid = (
                (rel >= -0.1)
                & (rel <= 20.1)
                & keep[r, parent]
            )

            gr = nearest(
                power[r, valid],
                rel[valid],
            )

            seqs.append(gr)

            logs.append(
                np.log10(
                    np.maximum(gr, 1e-30)
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
                    r,
                    tr,
                )
            )

            pulse_rel.append(
                blocks[
                    r,
                    2 * (tr - 1) : 2 * tr,
                ].reshape(-1)
                - origin
            )

    seqs = np.stack(seqs)
    logs = np.stack(logs)
    ids = np.asarray(ids)
    pulse_rel = np.stack(pulse_rel)

    roi0 = np.load(
        "outputs/hobbes20260831/ttl_aligned_corrected/"
        "ttl_aligned_analysis_results.npz"
    )

    anatomical = roi0[
        "lgn_mask_shifted"
    ].astype(bool)

    geom = np.load(
        "outputs/hobbes20260831/lgn_roi/lgn_roi_mask.npz"
    )

    xm, zm = (
        geom["mri_x_mm"],
        geom["mri_depth_mm"],
    )

    xn, zn = (
        geom["native_x_mm"][0],
        geom["native_depth_mm"][:, 0],
    )

    sigma_mm = 0.5 / 2.354820045

    sigma = (
        sigma_mm / abs(np.median(np.diff(zn))),
        sigma_mm / abs(np.median(np.diff(xn))),
    )

    smooth = np.full_like(logs, np.nan)

    for i in range(len(logs)):
        for k in range(len(GRID)):
            if np.isfinite(logs[i, k, 0, 0]):
                smooth[i, k] = gaussian_filter(
                    logs[i, k],
                    sigma=sigma,
                    mode="nearest",
                )

    contra = np.nanmean(
        smooth[:, :100],
        axis=1,
    )

    ipsi = np.nanmean(
        smooth[:, 100:],
        axis=1,
    )

    d = ipsi - contra

    obs, thr, null, siglabels, clusters = cluster_test(d)

    rows = []

    for (
        rank,
        sign,
        pv,
        mass,
        npix,
        peak,
        iy,
        ix,
    ) in clusters:
        q = siglabels == rank

        rows.append(
            {
                "cluster": rank,
                "direction": (
                    "ipsilateral > contralateral"
                    if sign > 0
                    else "contralateral > ipsilateral"
                ),
                "fwer_p": pv,
                "cluster_mass": mass,
                "pixels": npix,
                "peak_t": peak,
                "peak_x_mm": float(xm[iy, ix]),
                "peak_depth_mm": float(zm[iy, ix]),
                "x_min_mm": float(xm[q].min()),
                "x_max_mm": float(xm[q].max()),
                "depth_min_mm": float(zm[q].min()),
                "depth_max_mm": float(zm[q].max()),
                "overlap_anatomical_lgn_pixels": int(
                    np.sum(q & anatomical)
                ),
            }
        )

    with (
        OUT / "significant_clusters.csv"
    ).open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=(
                rows[0].keys()
                if rows
                else ["cluster"]
            ),
        )

        w.writeheader()
        w.writerows(rows)

    # Anatomical ROI signals with initial-baseline normalization
    # and run-offset removal.
    raw = np.nanmean(
        seqs[..., anatomical],
        axis=-1,
    )

    cbv = np.empty_like(raw)

    base_roi = np.mean(
        base[:, anatomical],
        axis=1,
    )

    for rr, run in enumerate(RUN_LABELS):
        cbv[
            ids[:, 0] == rr
        ] = 100 * (
            raw[ids[:, 0] == rr] / base_roi[rr] - 1
        )

    offsets = {}
    adj = cbv.copy()

    for r, run in enumerate(RUN_LABELS):
        off = float(
            np.nanmean(
                cbv[ids[:, 0] == r]
            )
        )

        offsets[r] = off
        adj[ids[:, 0] == r] -= off

    mean = np.nanmean(adj, axis=0)

    n = np.sum(
        np.isfinite(adj),
        axis=0,
    )

    sem = (
        np.nanstd(
            adj,
            axis=0,
            ddof=1,
        )
        / np.sqrt(n)
    )

    c10 = np.nanmean(
        adj[:, :100],
        axis=1,
    )

    i10 = np.nanmean(
        adj[:, 100:],
        axis=1,
    )

    test = ttest_rel(i10, c10)

    # QC and primary figures.
    fig, ax = plt.subplots(
        n_runs,
        1,
        figsize=(13, 7),
        sharex=True,
        constrained_layout=True,
        squeeze=False,
    )

    ax = ax.ravel()

    for r, a in enumerate(ax):
        a.plot(
            np.arange(1, n_PARENT + 1),
            rz[r],
            color="black",
            lw=1,
        )

        a.axhline(
            5,
            color="red",
            ls="--",
        )

        a.set_ylabel(
            f"Run {RUN_LABELS[r]}\n|robust z|"
        )

        for tr in reject[r]:
            a.axvspan(
                11 + (tr - 1) * 20,
                10 + tr * 20,
                color="red",
                alpha=0.12,
            )

    ax[-1].set_xlabel(
        "Parent BF frame (1 based)"
    )

    fig.savefig(
        OUT / "artifact_censoring_qc.png",
        dpi=220,
    )

    plt.close(fig)

    fig, ax = plt.subplots(
        figsize=(12, 5),
        constrained_layout=True,
    )

    a = GRID < 10
    b = ~a

    ax.fill_between(
        GRID[a],
        mean[a] - sem[a],
        mean[a] + sem[a],
        color="blue",
        alpha=0.2,
    )

    ax.plot(
        GRID[a],
        mean[a],
        color="blue",
        lw=2.5,
        label="Contralateral",
    )

    ax.fill_between(
        GRID[b],
        mean[b] - sem[b],
        mean[b] + sem[b],
        color="red",
        alpha=0.2,
    )

    ax.plot(
        np.r_[GRID[a][-1], GRID[b]],
        np.r_[mean[a][-1], mean[b]],
        color="red",
        lw=2.5,
        label="Ipsilateral",
    )

    for x in np.arange(0, 8.21, 0.2):
        ax.axvspan(
            x,
            x + 0.02,
            color="#f2c94c",
            alpha=0.35,
        )

    for x in np.arange(10, 18.21, 0.2):
        ax.axvspan(
            x,
            x + 0.02,
            color="#f2c94c",
            alpha=0.35,
        )

    ax.axhline(
        0,
        color="#888",
        lw=0.8,
    )

    ax.axvline(
        10,
        color="#555",
        lw=1,
    )

    ax.set(
        xlim=(0, 20),
        xlabel="Time (s)",
        ylabel="Run-adjusted anatomical LGN CBV change (%)",
        title="Theta-burst fUSI response: Subject H",
    )

    ax.legend()

    fig.savefig(
        OUT / "anatomical_lgn_cbv_trace.png",
        dpi=240,
    )

    fig.savefig(
        OUT / "anatomical_lgn_cbv_trace.pdf"
    )

    plt.close(fig)

    mri, xo, zo = mri_slice()

    fig, ax = plt.subplots(
        figsize=(7, 10),
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

    im = ax.pcolormesh(
        xm,
        zm,
        obs,
        cmap="coolwarm",
        vmin=-6,
        vmax=6,
        alpha=0.78,
        shading="nearest",
    )

    for row in rows:
        q = siglabels == row["cluster"]

        ax.contour(
            xm,
            zm,
            q.astype(float),
            [0.5],
            colors=(
                "yellow"
                if "ipsilateral >" in row["direction"]
                else "cyan"
            ),
            linewidths=1.5,
        )

    ax.contour(
        xm,
        zm,
        anatomical.astype(float),
        [0.5],
        colors="lime",
        linewidths=1.5,
    )

    ax.set(
        xlim=(-5, 25),
        ylim=(50, -5),
        xlabel="Distance from midline (mm)",
        ylabel="MRI depth (mm)",
        title=(
            "Ipsilateral − contralateral: entire 10-second period\n"
            "cluster-mass FWER p < 0.05"
        ),
    )

    fig.colorbar(
        im,
        ax=ax,
        label="paired t statistic",
    )

    fig.savefig(
        OUT / "whole_10s_cluster_tmap_on_mri.png",
        dpi=240,
    )

    plt.close(fig)

    np.savez_compressed(
        OUT / "theta_analysis_results.npz",
        time_s=GRID,
        trial_run_and_number=ids,
        pulse_times_s=pulse_rel,
        power_trials=seqs,
        log_change_trials=logs,
        anatomical_lgn_mask=anatomical,
        anatomical_lgn_cbv_trials=adj,
        anatomical_lgn_raw_power_trials=raw,
        ipsilateral_log_maps=ipsi,
        contralateral_log_maps=contra,
        observed_t=obs,
        significant_cluster_labels=siglabels,
        null_max_cluster_mass=null,
        cluster_forming_t=thr,
        mri_x_mm=xm,
        mri_depth_mm=zm,
        available_n=n,
    )

    summary = {
        "selected_original_runs": RUN_LABELS,
        "rate_hz": 10,
        "retained_pairs": len(ids),
        "retained_by_run": {
            str(int(r)): int(
                np.sum(ids[:, 0] == r)
            )
            for r, run in enumerate(RUN_LABELS)
        },
        "artifact_threshold": (
            "parent maximum absolute whole-field "
            "robust z >= 5"
        ),
        "artifact_frames_1_based": artifact_frames,
        "removed_pairs": {
            str(k): v
            for k, v in reject.items()
        },
        "condition_order": [
            "contralateral/right/control",
            "ipsilateral/left/active",
        ],
        "run_offsets_percent": offsets,
        "anatomical_lgn_whole_10s": {
            "contralateral_mean_percent": float(
                c10.mean()
            ),
            "ipsilateral_mean_percent": float(
                i10.mean()
            ),
            "paired_t": float(
                test.statistic
            ),
            "df": len(i10) - 1,
            "paired_p": float(
                test.pvalue
            ),
        },
        "spatial_smoothing_fwhm_mm": 0.5,
        "cluster_forming_two_sided_p": CLUSTER_P,
        "cluster_forming_abs_t": thr,
        "permutations": N_PERM,
        "significant_clusters": rows,
    }

    (
        OUT / "summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(
        json.dumps(
            summary,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

