#!/usr/bin/env python3

"""Decode fUSI, sonication, and Arduino target triggers from Intan RHS files."""

from __future__ import annotations

import csv
import json
import struct
from pathlib import Path

import numpy as np
import config
from scipy.io import loadmat


EEG_DIR = Path("/Users/jpalmer/Downloads/HobbesfusiFUSAug312026/Hobbes20260831/EEG")
HOBBES = Path("/System/Volumes/Data/v/raid1b/backup/jpalmer/FUS/fUSI_FUS_analysis/fUSI_FUS_code/outputs/hobbes20260831/hobbes20230207.mat")
OUT = config.OUTPUT_DIR
BLOCK = 128


def unpack(f, fmt):
    size = struct.calcsize("<" + fmt)
    data = f.read(size)

    if len(data) != size:
        raise EOFError

    return struct.unpack("<" + fmt, data)[0]


def qstring(f):
    n = unpack(f, "I")

    if n == 0xFFFFFFFF:
        return ""

    return f.read(n).decode(
        "utf-16-le",
        errors="replace",
    )


def read_header(f):
    if unpack(f, "I") != 0xD69127AC:
        raise ValueError("Not an Intan RHS file")

    major, minor = unpack(f, "h"), unpack(f, "h")
    sample_rate = unpack(f, "f")

    f.read(
        2
        + 4 * 8
        + 2
        + 4 * 2
        + 2 * 2
        + 4 * 3
    )

    notes = [qstring(f) for _ in range(3)]

    dc_saved, eval_mode = unpack(f, "h"), unpack(f, "h")
    reference = qstring(f)

    groups = unpack(f, "h")

    channels = {
        i: []
        for i in range(7)
    }

    for group in range(groups):
        group_name, prefix = (
            qstring(f),
            qstring(f),
        )

        enabled, nchan, namp = (
            unpack(f, "h"),
            unpack(f, "h"),
            unpack(f, "h"),
        )

        if nchan > 0 and enabled > 0:
            for _ in range(nchan):
                native, custom = (
                    qstring(f),
                    qstring(f),
                )

                native_order, custom_order = (
                    unpack(f, "h"),
                    unpack(f, "h"),
                )

                signal_type, channel_enabled = (
                    unpack(f, "h"),
                    unpack(f, "h"),
                )

                chip_channel = unpack(f, "h")

                command_stream, board_stream = (
                    unpack(f, "h"),
                    unpack(f, "h"),
                )

                f.read(
                    2 * 4
                    + 4 * 2
                )

                if channel_enabled:
                    channels[signal_type].append(
                        {
                            "native": native,
                            "custom": custom,
                            "native_order": native_order,
                        }
                    )

    return {
        "major": major,
        "minor": minor,
        "sample_rate": sample_rate,
        "dc_saved": dc_saved,
        "channels": channels,
        "header_bytes": f.tell(),
        "notes": notes,
        "reference": reference,
    }


def read_trigger_channels(path):
    with path.open("rb") as f:
        h = read_header(f)

        n_amp = len(h["channels"][0])
        n_adc = len(h["channels"][3])
        n_dac = len(h["channels"][4])

        n_din = len(h["channels"][5])
        n_dout = len(h["channels"][6])

        per = (
            BLOCK * 4
            + BLOCK * (6 if h["dc_saved"] else 4) * n_amp
        )

        per += (
            BLOCK * 2 * (n_adc + n_dac)
            + (BLOCK * 2 if n_din else 0)
            + (BLOCK * 2 if n_dout else 0)
        )

        remaining = (
            path.stat().st_size
            - h["header_bytes"]
        )

        if remaining % per:
            raise ValueError(
                f"Non-integral block count in {path.name}: "
                f"remainder {remaining % per}"
            )

        nblocks = remaining // per

        timestamps = np.empty(
            nblocks * BLOCK,
            np.int32,
        )

        adc = np.empty(
            (n_adc, nblocks * BLOCK),
            np.uint16,
        )

        dinraw = (
            np.empty(
                nblocks * BLOCK,
                np.uint16,
            )
            if n_din
            else None
        )

        pos = 0

        for _ in range(nblocks):
            timestamps[
                pos : pos + BLOCK
            ] = np.fromfile(
                f,
                "<i4",
                BLOCK,
            )

            f.seek(
                BLOCK
                * (6 if h["dc_saved"] else 4)
                * n_amp,
                1,
            )

            if n_adc:
                block = np.fromfile(
                    f,
                    "<u2",
                    BLOCK * n_adc,
                ).reshape(
                    BLOCK,
                    n_adc,
                )

                adc[
                    :,
                    pos : pos + BLOCK,
                ] = block.T

            f.seek(
                BLOCK * 2 * n_dac,
                1,
            )

            if n_din:
                dinraw[
                    pos : pos + BLOCK
                ] = np.fromfile(
                    f,
                    "<u2",
                    BLOCK,
                )

            if n_dout:
                f.seek(
                    BLOCK * 2,
                    1,
                )

            pos += BLOCK

    digital = np.vstack(
        [
            (
                dinraw
                & (1 << c["native_order"])
            )
            > 0
            for c in h["channels"][5]
        ]
    )

    analog = (
        312.5e-6
        * (adc.astype(np.float32) - 32768)
    )

    return timestamps, digital, analog, h


def rising(binary):
    return np.flatnonzero(
        np.diff(
            binary.astype(np.int8),
            prepend=0,
        )
        > 0
    )


def runs(binary):
    padded = np.r_[
        False,
        binary,
        False,
    ]

    return (
        np.flatnonzero(
            np.diff(
                padded.astype(np.int8)
            )
            == 1
        ),
        np.flatnonzero(
            np.diff(
                padded.astype(np.int8)
            )
            == -1
        ),
    )


def decode_arduino_messages(
    time,
    digital,
    bit_width=0.001,
):
    """Port of processArduino/dig2num, retaining message onset times."""

    on = rising(digital)
    dt = np.median(np.diff(time))

    window_n = (
        12
        * (
            round(bit_width / dt)
            + round(200e-6 / dt)
        )
        + 2
    )

    values, times = [], []
    cursor = 0

    while cursor < len(on):
        start = on[cursor]
        stop = min(
            len(digital),
            start + window_n,
        )

        segment = digital[start:stop]
        segtime = time[start:stop]

        rises, falls = runs(segment)

        if len(rises) == len(falls) and len(rises):
            bits = []

            for j in range(len(rises)):
                bits += [
                    1
                ] * int(
                    np.floor(
                        (
                            segtime[falls[j] - 1]
                            - segtime[rises[j]]
                        )
                        / bit_width
                        + 1e-6
                    )
                )

                if j + 1 < len(rises):
                    bits += [
                        0
                    ] * int(
                        np.floor(
                            (
                                segtime[rises[j + 1]]
                                - segtime[falls[j] - 1]
                            )
                            / bit_width
                            + 1e-6
                        )
                    )

            payload = bits[1:-1]

            if payload:
                values.append(
                    int(
                        "".join(
                            map(str, payload)
                        ),
                        2,
                    )
                )

                times.append(
                    time[start]
                )

        cursor = np.searchsorted(
            on,
            stop,
            side="right",
        )

    return (
        np.asarray(times),
        np.asarray(values, dtype=int),
    )


def main():
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    files = sorted(
        EEG_DIR.glob("hobbes*.rhs")
    )

    all_t = []
    all_d = []
    all_a = []
    headers = []

    for i, path in enumerate(files):
        ts, dig, alg, header = (
            read_trigger_channels(path)
        )

        all_t.append(
            ts.astype(np.int64)
        )

        all_d.append(dig)
        all_a.append(alg)
        headers.append(header)

        print(
            f"{i + 1:02d}/{len(files)} "
            f"{path.name}: {len(ts)} samples",
            flush=True,
        )

    ticks = np.concatenate(all_t)

    # Intan timestamps are continuous across one-minute files.
    time = (
        ticks - ticks[0]
    ) / headers[0]["sample_rate"]

    digital = np.concatenate(
        all_d,
        axis=1,
    )

    analog = np.concatenate(
        all_a,
        axis=1,
    )

    print("\n===== CHANNEL DIAGNOSTIC =====")
    print("Digital shape:", digital.shape)
    print("Analog shape:", analog.shape)
    print("Digital channels:", headers[0]["channels"][5])
    print("Analog channels:", headers[0]["channels"][3])
    print("==============================\n")

    fusi_edges = time[
        rising(digital[1])
    ]

    son_binary = analog[0] > 1.0
    son_edges = time[
        rising(son_binary)
    ]

    arduino_t, arduino_num = (
        decode_arduino_messages(
            time,
            digital[0],
        )
    )

    # The recording itself shows the second pulse about 10 s after the first
    # (despite the accompanying email describing a 5 s baseline).
    pair_gap = np.diff(son_edges)

    actual = son_edges[1:][
        (pair_gap > 9.5)
        & (pair_gap < 10.5)
    ]

    baseline = son_edges[:-1][
        (pair_gap > 9.5)
        & (pair_gap < 10.5)
    ]

    mat = loadmat(
        HOBBES,
        squeeze_me=True,
        struct_as_record=False,
    )

    table = [
        tuple(
            np.asarray(
                x.target,
                float,
            )
        )
        for x in np.atleast_1d(
            mat["sonicationTable"]
        )
    ]

    hobbes_target = np.array(
        [
            table.index(
                tuple(
                    np.asarray(
                        x.sonication.target,
                        float,
                    )
                )
            )
            + 1
            for x in np.atleast_1d(
                mat["logSys"]
            )
        ]
    )

    # Each trial begins at the baseline trigger. Its 20 BF acquisitions straddle
    # the actual sonication: approximately 10 before and 10 after.
    trial_for_fusi = (
        np.searchsorted(
            baseline,
            fusi_edges,
            side="right",
        )
        - 1
    )

    decoded_target = np.full(
        len(actual),
        -1,
        int,
    )

    decoded_messages = []

    for i, onset in enumerate(actual):
        prior = np.flatnonzero(
            (arduino_t < onset)
            & (arduino_t >= onset - 20)
        )

        vals = arduino_num[prior]
        vals = vals[vals != 0]

        decoded_messages.append(
            vals.tolist()
        )

        if len(vals):
            decoded_target[i] = vals[-1]

    rows = []

    for i, onset in enumerate(actual):
        idx = np.flatnonzero(
            trial_for_fusi == i
        )

        rows.append(
            {
                "trial": i + 1,
                "baseline_trigger_s": baseline[i],
                "sonication_s": onset,
                "decoded_target": decoded_target[i],
                "hobbes_target": (
                    hobbes_target[i]
                    if i < len(hobbes_target)
                    else -1
                ),
                "target_match": int(
                    i < len(hobbes_target)
                    and decoded_target[i]
                    == hobbes_target[i]
                ),
                "arduino_messages": ";".join(
                    map(
                        str,
                        decoded_messages[i],
                    )
                ),
                "first_bf_ordinal": (
                    int(idx[0] + 1)
                    if len(idx)
                    else -1
                ),
                "last_bf_ordinal": (
                    int(idx[-1] + 1)
                    if len(idx)
                    else -1
                ),
                "bf_count": len(idx),
                "first_bf_delay_s": (
                    fusi_edges[idx[0]] - onset
                    if len(idx)
                    else np.nan
                ),
                "last_bf_delay_s": (
                    fusi_edges[idx[-1]] - onset
                    if len(idx)
                    else np.nan
                ),
            }
        )

    with (
        OUT / "trial_timing.csv"
    ).open(
        "w",
        newline="",
    ) as f:
        w = csv.DictWriter(
            f,
            fieldnames=list(rows[0]),
        )

        w.writeheader()
        w.writerows(rows)

    np.savez_compressed(
        OUT / "trigger_times.npz",
        fusi_trigger_s=fusi_edges,
        sonication_trigger_s=actual,
        baseline_trigger_s=baseline,
        arduino_message_s=arduino_t,
        arduino_value=arduino_num,
        decoded_target=decoded_target,
        hobbes_target=hobbes_target,
    )

    summary = {
        "rhs_files": len(files),
        "sample_rate_hz": headers[0]["sample_rate"],
        "recording_duration_s": float(time[-1]),
        "digital_channels": headers[0]["channels"][5],
        "analog_channels": headers[0]["channels"][3],
        "all_sonication_edges": len(son_edges),
        "baseline_actual_pairs": len(actual),
        "fusi_rising_edges": len(fusi_edges),
        "arduino_messages": len(arduino_num),
        "arduino_values": arduino_num.tolist(),
        "decoded_targets": decoded_target.tolist(),
        "hobbes_targets": hobbes_target.tolist(),
        "target_matches": int(
            np.sum(
                decoded_target[:len(hobbes_target)]
                == hobbes_target[:len(decoded_target)]
            )
        ),
        "bf_counts_per_trial": [
            int(r["bf_count"])
            for r in rows
        ],
        "first_bf_delay_s": [
            float(r["first_bf_delay_s"])
            for r in rows
        ],
        "timestamp_gaps": int(
            np.sum(
                np.diff(ticks) != 1
            )
        ),
    }

    (
        OUT / "summary.json"
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
