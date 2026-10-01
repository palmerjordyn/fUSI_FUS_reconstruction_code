#!/usr/bin/env python3

from pathlib import Path
import json
import tempfile
import zipfile
import numpy as np
import sys
import matplotlib.pyplot as plt
import config


sys.path.insert(0, str(Path(__file__).parent))

from decode_intan_fusi_timing import read_trigger_channels, rising

if (config.ZIP):
    source = config.ZIP_PATH
else:
    source = config.DATA_PATH


OUT = config.OUTPUT_DIR / "design"
OUT.mkdir(parents=True, exist_ok=True)

LEN_GROUPS = 0

# def read_take(z, take):
#     names = sorted(
#         n
#         for n in z.namelist()
#         if f"/hobbesFusiFus20sblock15rep_take{take}_" in n and n.endswith(".rhs")
#     )

#     ticks = []
#     digital = []
#     header = None

#     for name in names:
#         with tempfile.NamedTemporaryFile(suffix=".rhs") as f:
#             f.write(z.read(name))
#             f.flush()

#             t, d, _, h = read_trigger_channels(Path(f.name))

#         ticks.append(t.astype(np.int64))
#         digital.append(d)
#         header = header or h

#     tick = np.concatenate(ticks)
#     dig = np.concatenate(digital, axis=1)

#     return (
#         (tick - tick[0]) / header["sample_rate"],
#         dig,
#         header,
#         names,
#     )

def read_take(source, acquisition_index, is_zip):
    identifier = config.RHS_IDENTIFIER[acquisition_index]

    if is_zip:
        names = sorted(
            n
            for n in source.namelist()
            if (
                f"Hobbes{config.DATE}/" in n
                and identifier in Path(n).name
                and n.endswith(".rhs")
            )
        )


    else:
        eeg_dir = source / f"Hobbes{config.DATE}"

        names = sorted(
            p
            for p in eeg_dir.rglob("*.rhs")
            if (
                p.is_file()
                and identifier in p.name
            )
        )

    if not names:
        raise RuntimeError(
            f"No RHS files found for identifier '{identifier}'"
        )

    print(
        f"\nAcquisition {acquisition_index + 1}: "
        f"identifier='{identifier}', found {len(names)} RHS files"
    )

    if config.VERBOSE:
        for name in names:
            print(f"  {Path(name).name}")

    ticks = []
    digital = []
    header = None

    for name in names:

        if is_zip:
            with tempfile.NamedTemporaryFile(suffix=".rhs") as f:
                f.write(source.read(name))
                f.flush()

                t, d, _, h = read_trigger_channels(Path(f.name))

        else:
            t, d, _, h = read_trigger_channels(name)

        ticks.append(t.astype(np.int64))
        digital.append(d)
        header = header or h

    tick = np.concatenate(ticks)
    dig = np.concatenate(digital, axis=1)

    return (
        (tick - tick[0]) / header["sample_rate"],
        dig,
        header,
        names,
    )
# def trains(edges):
#     return [
#         g
#         for g in np.split(
#             edges,
#             np.flatnonzero(np.diff(edges) > 1) + 1
#         )
#         if len(g) == 42
#     ]

def trains(edges):
    if ((edges[17]-edges[16])>2):
        wait = 3
    else:
        wait = 1

    groups = np.split(
        edges,
        np.flatnonzero(np.diff(edges) > wait) + 1
    )

    length = len(groups[2])
    LEN_GROUPS = length

    if config.VERBOSE:
        print("Group lengths:")
        for i, g in enumerate(groups):
            print(f"  Group {i}: {len(g)} edges")

    return [
        g
        for g in groups
        if len(g) == length
    ]


def main():

    all_fusi = []
    all_pulses = []
    summ = []

    # Open either the ZIP or the normal data directory
    if config.ZIP:
        source = zipfile.ZipFile(config.ZIP_PATH)
    else:
        source = config.DATA_PATH

    try:
        for i, (take, acq) in enumerate(zip(config.RHS_IDENTIFIER, config.ACQS)):

            # i = 0 -> RHS_IDENTIFIER[0]
            # i = 1 -> RHS_IDENTIFIER[1]
            # i = 2 -> RHS_IDENTIFIER[2]
            time, d, h, names = read_take(
                source,
                i,
                config.ZIP
            )

            e0 = time[rising(d[0])]
            fu = time[rising(d[1])]

            if config.VERBOSE:
                print("First 100 e0 intervals:")
                print(np.diff(e0[:100]))
                print(e0[:100])
                print(f"number of fUSI edges: {len(fu)}")
                print("first 100 fUSI intervals:")
                print(np.diff(fu[:101]))

                plt.figure(figsize=(12, 3))

                plt.vlines(
                    e0[(e0 >= 150) & (e0 <= 200)],
                    0,
                    1,
                    label="FUS"
                )

                plt.vlines(
                    fu[(fu >= 150) & (fu <= 200)],
                    1,
                    2,
                    label="fUSI"
                )

                plt.xlabel("Time (s)")
                plt.yticks([0.5, 1.5], ["FUS", "fUSI"])
                plt.legend()
                plt.tight_layout()
                plt.show()

            groups = trains(e0)

            if len(groups) != 30 or len(fu) != config.BF[i]:
                raise RuntimeError(
                    f"take {take}: groups={len(groups)}, "
                    f"fUSI={len(fu)}"
                )

            pulse = np.stack(groups)

            all_fusi.append(fu)
            all_pulses.append(pulse)

            starts = pulse[:, 0]
            within = np.diff(pulse, axis=1)

            summ.append(
                {
                    "take": take,
                    "acquisition": acq,
                    "rhs_identifier": config.RHS_IDENTIFIER[i],
                    "rhs_files": len(names),
                    "duration_s": float(time[-1]),
                    "sample_rate_hz": float(h["sample_rate"]),
                    "fusi_ttls": len(fu),
                    "stimulation_blocks": len(groups),
                    "recorded_edges_per_block": len(groups[2]),
                    "within_block_interval_mean_s": float(
                        within.mean()
                    ),
                    "within_block_duration_s": float(
                        np.median(
                            pulse[:, -1] - pulse[:, 0]
                        )
                    ),
                    "between_block_start_interval_s": float(
                        np.median(np.diff(starts))
                    ),
                    "first_block_eeg_s": float(starts[0]),
                    "first_fusi_eeg_s": float(fu[0]),
                }
            )

    finally:
        if config.ZIP:
            source.close()

    np.savez_compressed(
        OUT / "selected_run_timing.npz",
        take_numbers=np.asarray(config.RHS_IDENTIFIER),
        acquisition_names=np.asarray(config.ACQS),
        fusi_ttl_s=np.stack(all_fusi),
        block_pulse_ttl_s=np.stack(all_pulses),
    )

    summary = {
        "study": "Hobbes Oct 31 2026 theta burst",
        "selected_original_runs": config.RHS_IDENTIFIER,
        "condition_order": (
            "first 10 s right LGN/control/contralateral; "
            "second 10 s left LGN/active/ipsilateral"
        ),
        "pairs_per_run": 15,
        "blocks_per_run": 30,
        "recorded_stimulation_edges_per_block": LEN_GROUPS,
        "note": "no note",
        "runs": summ,
    }

    (
        OUT / "design_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

