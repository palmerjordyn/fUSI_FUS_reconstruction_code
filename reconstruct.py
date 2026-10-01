#!/usr/bin/env python3

"""Stream BF data and reconstruct 10-Hz PD using SVD or Butterworth filtering."""

from pathlib import Path

import json
import re
import time
import zipfile
import config

import numpy as np

from scipy.linalg import svd
from scipy.signal import butter, sosfiltfilt


# ============================================================
# Select data source
# ============================================================

if config.ZIP:
    source = config.ZIP_PATH
else:
    source = config.DATA_PATH


# ============================================================
# Check reconstruction type
# ============================================================

if config.RECONSTRUCTION_TYPE not in ("SVD", "BUTTER"):
    raise ValueError(
        f"Invalid RECONSTRUCTION_TYPE: "
        f"{config.RECONSTRUCTION_TYPE!r}. "
        f"Expected 'SVD' or 'BUTTER'."
    )


# ============================================================
# Output
# ============================================================

OUT = (
    config.OUTPUT_DIR
    / config.RECONSTRUCTION_TYPE
    / "reconstruction_10hz"
)

RUNS = config.ACQS


# ============================================================
# Reconstruction parameters
# ============================================================

NZ, NX, NT = 339, 38, 200

PER_PARENT = 8
WINDOW = 25

numBF = max(config.BF)


# ============================================================
# Butterworth parameters
# ============================================================

IQ_RATE_HZ = 250.0
BUTTERWORTH_CUTOFF_HZ = 70.0
BUTTERWORTH_ORDER = 4


# ============================================================
# Find BF files
# ============================================================

def members(source, run, is_zip):

    pat = re.compile(
        r"P7-4__BF_(\d+)\.bin$"
    )

    found = []

    if is_zip:

        for name in source.namelist():

            if f"/{run}/" not in name:
                continue

            m = pat.search(name)

            if m:
                found.append(
                    (int(m.group(1)), name)
                )

    else:

        run_dir = source / run

        for path in run_dir.iterdir():

            if not path.is_file():
                continue

            m = pat.search(path.name)

            if m:
                found.append(
                    (int(m.group(1)), path)
                )

    return [
        n
        for _, n in sorted(found)
    ]


# ============================================================
# Main
# ============================================================

def main():

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    cache = (
        OUT /
        "power_doppler_10hz.npz"
    )

    if cache.exists():

        print(
            f"already exists: {cache}"
        )

        return


    # ========================================================
    # Design Butterworth filter
    # ========================================================

    if config.RECONSTRUCTION_TYPE == "BUTTER":

        sos = butter(
            BUTTERWORTH_ORDER,
            BUTTERWORTH_CUTOFF_HZ,
            btype="highpass",
            fs=IQ_RATE_HZ,
            output="sos"
        )

        print(
            f"Butterworth high-pass: "
            f"order={BUTTERWORTH_ORDER}, "
            f"cutoff={BUTTERWORTH_CUTOFF_HZ} Hz, "
            f"sampling={IQ_RATE_HZ} Hz"
        )

    else:

        sos = None


    # ========================================================
    # Allocate output arrays
    # ========================================================

    n_runs = len(RUNS)

    power = np.empty(
        (
            n_runs,
            numBF * PER_PARENT,
            NZ,
            NX
        ),
        np.float32
    )

    parent_bmode = np.empty(
        (
            n_runs,
            numBF,
            NZ,
            NX
        ),
        np.float32
    )


    started = time.time()

    inventory = {}


    # ========================================================
    # Open source
    # ========================================================

    if config.ZIP:

        source = zipfile.ZipFile(
            config.ZIP_PATH
        )

    else:

        source = config.DATA_PATH


    try:

        # ====================================================
        # Process each acquisition
        # ====================================================

        for rr, run in enumerate(RUNS):

            checkpoint = (
                OUT /
                f"checkpoint_run_{rr + 1}.npz"
            )


            # =================================================
            # Load checkpoint
            # =================================================

            if checkpoint.exists():

                q = np.load(
                    checkpoint
                )

                power[rr] = q[
                    "power_doppler"
                ]

                parent_bmode[rr] = q[
                    "parent_bmode"
                ]

                inventory[run] = {
                    "source_bf_files": len(
                        members(
                            source,
                            run,
                            config.ZIP
                        )
                    ),
                    "used_bf_files": config.BF[rr],
                    "loaded_checkpoint": True
                }

                print(
                    f"run {rr + 1}/{n_runs} "
                    f"{run}: "
                    f"loaded checkpoint",
                    flush=True
                )

                continue


            # =================================================
            # Find BF files
            # =================================================

            all_members = members(
                source,
                run,
                config.ZIP
            )


            inventory[run] = {
                "source_bf_files": len(
                    all_members
                ),
                "used_bf_files": config.BF[rr]
            }


            if len(all_members) < config.BF[rr]:

                raise RuntimeError(
                    f"{run}: "
                    f"only {len(all_members)} BF files"
                )


            # =================================================
            # Process BF files
            # =================================================

            for i, name in enumerate(
                all_members[:config.BF[rr]]
            ):

                # ---------------------------------------------
                # Read BF file
                # ---------------------------------------------

                if config.ZIP:

                    raw = source.read(
                        name
                    )

                else:

                    raw = name.read_bytes()


                # ---------------------------------------------
                # Decode complex BF data
                # ---------------------------------------------

                cube = np.frombuffer(
                    raw,
                    dtype="<c16"
                )


                if cube.size != NZ * NX * NT:

                    raise RuntimeError(
                        f"{name}: "
                        f"{cube.size} complex samples"
                    )


                cube = cube.reshape(
                    (
                        NZ,
                        NX,
                        NT
                    ),
                    order="F"
                )


                # =================================================
                # Parent B-mode
                # =================================================

                parent_bmode[
                    rr,
                    i
                ] = np.median(
                    np.abs(cube),
                    axis=2
                )


                # =================================================
                # SVD reconstruction
                # =================================================
                if config.RECONSTRUCTION_TYPE == "SVD":

                    # -----------------------------------------
                    # Convert to [spatial_pixels, time]
                    # -----------------------------------------

                    matrix = cube.reshape(
                        (
                            NZ * NX,
                            NT
                        ),
                        order="F"
                    )


                    # -----------------------------------------
                    # SVD
                    # -----------------------------------------

                    u, s, vh = svd(
                        matrix,
                        full_matrices=False,
                        check_finite=False,
                        lapack_driver="gesdd"
                    )


                    # -----------------------------------------
                    # Remove components 0 through 29
                    # -----------------------------------------

                    filtered = matrix.copy()

                    for rank in range(30):

                        filtered -= np.multiply.outer(
                            u[:, rank] * s[rank],
                            vh[rank]
                        )


                    # -----------------------------------------
                    # Restore [NZ, NX, NT]
                    # -----------------------------------------

                    filtered = filtered.reshape(
                        (
                            NZ,
                            NX,
                            NT
                        ),
                        order="F"
                    )


                # =================================================
                # Butterworth reconstruction
                # =================================================

                elif config.RECONSTRUCTION_TYPE == "BUTTER":
                    # -----------------------------------------
                    # Subtract first IQ frame
                    # -----------------------------------------

                    cube = (
                        cube
                        - cube[:, :, 0][:, :, np.newaxis]
                    )


                    # -----------------------------------------
                    # Convert to [spatial_pixels, time]
                    # -----------------------------------------

                    matrix = cube.reshape(
                        (
                            NZ * NX,
                            NT
                        ),
                        order="F"
                    )


                    # -----------------------------------------
                    # 70-Hz Butterworth high-pass filter
                    # -----------------------------------------

                    filtered = sosfiltfilt(
                        sos,
                        matrix,
                        axis=1
                    )


                    # -----------------------------------------
                    # Restore [NZ, NX, NT]
                    # -----------------------------------------

                    filtered = filtered.reshape(
                        (
                            NZ,
                            NX,
                            NT
                        ),
                        order="F"
                    )


                # =================================================
                # Divide 200 frames into 8 × 25-frame PD images
                # =================================================

                for j in range(PER_PARENT):

                    q = filtered[
                        :,
                        :,
                        j * WINDOW:(j + 1) * WINDOW
                    ]


                    # -----------------------------------------
                    # Power Doppler
                    # -----------------------------------------

                    power[
                        rr,
                        i * PER_PARENT + j
                    ] = np.sum(
                        np.abs(q) ** 2,
                        axis=2,
                        dtype=np.float64
                    )


                # =================================================
                # Progress
                # =================================================

                if (
                    i == 0
                    or (i + 1) % 20 == 0
                ):

                    print(
                        f"run {rr + 1}/{n_runs} "
                        f"{run}: "
                        f"{i + 1}/{config.BF[rr]} "
                        f"({time.time() - started:.1f}s)",
                        flush=True
                    )


            # =================================================
            # Save checkpoint
            # =================================================

            np.savez_compressed(
                checkpoint,
                power_doppler=power[rr],
                parent_bmode=parent_bmode[rr]
            )


    finally:

        if config.ZIP:

            source.close()


    # ========================================================
    # Across-run temporal median
    # ========================================================

    median = np.median(
        power,
        axis=(0, 1)
    )


    finite = np.isfinite(
        median
    )


    support = (
        finite
        &
        (
            median
            >
            np.percentile(
                median[finite],
                40
            )
        )
    )


    # ========================================================
    # Robust z-score
    # ========================================================

    robust_z = np.empty(
        (
            n_runs,
            numBF * PER_PARENT
        ),
        np.float32
    )


    for r in range(n_runs):

        g = np.median(
            np.log10(
                np.maximum(
                    power[r][:, support],
                    1e-30
                )
            ),
            axis=1
        )


        center = np.median(
            g
        )


        mad = np.median(
            np.abs(
                g - center
            )
        )


        robust_z[r] = (
            .67448975
            * (g - center)
            / max(
                mad,
                np.finfo(float).eps
            )
        )


    # ========================================================
    # B-mode
    # ========================================================

    bmode = np.median(
        parent_bmode,
        axis=(0, 1)
    )


    ref = np.percentile(
        bmode[bmode > 0],
        99.9
    )


    bmode_db = np.clip(
        20
        * np.log10(
            np.maximum(
                bmode / ref,
                1e-6
            )
        ),
        -60,
        0
    ).astype(
        np.float32
    )


    # ========================================================
    # Save reconstruction
    # ========================================================

    np.savez_compressed(
        cache,
        power_doppler=power,
        whole_field_robust_z=robust_z,
        support_mask=support,
        images_per_parent=PER_PARENT,
        iq_frames_per_image=WINDOW,
        subimage_center_offsets_s=(
            np.arange(PER_PARENT) + 0.5
        ) / 10,
        run_names=np.asarray(RUNS)
    )


    # ========================================================
    # Save B-mode
    # ========================================================

    np.savez_compressed(
        OUT /
        "temporal_median_bmode.npz",
        bmode_db=bmode_db,
        envelope=bmode.astype(
            np.float32
        )
    )


    # ========================================================
    # Summary
    # ========================================================

    summary = {

        "source": (
            str(config.ZIP_PATH)
            if config.ZIP
            else str(config.DATA_PATH)
        ),

        "runs": RUNS,

        "inventory": inventory,

        "reconstruction_type":
            config.RECONSTRUCTION_TYPE,

        "rate_hz": 10,

        "parent_frames_per_run":
            numBF,

        "images_per_parent":
            PER_PARENT,

        "iq_frames_per_image":
            WINDOW,

        "shape":
            list(power.shape),

        "elapsed_processing_seconds":
            time.time() - started
    }


    # ========================================================
    # SVD-specific summary information
    # ========================================================

    if config.RECONSTRUCTION_TYPE == "SVD":

        summary.update({

            "svd_removed_components_zero_based":
                [0, 29]

        })


    # ========================================================
    # Butterworth-specific summary information
    # ========================================================

    elif config.RECONSTRUCTION_TYPE == "Butter":

        summary.update({

            "iq_sampling_rate_hz":
                IQ_RATE_HZ,

            "butterworth_cutoff_hz":
                BUTTERWORTH_CUTOFF_HZ,

            "butterworth_order":
                BUTTERWORTH_ORDER,

            "filter":
                "70-Hz Butterworth high-pass",

            "filter_method":
                "sosfiltfilt",

            "first_frame_subtracted":
                True

        })


    # ========================================================
    # Write summary
    # ========================================================

    (
        OUT /
        "reconstruction_summary.json"
    ).write_text(
        json.dumps(
            summary,
            indent=2
        )
        + "\n"
    )


    print(
        json.dumps(
            summary,
            indent=2
        )
    )


if __name__ == "__main__":
    main()