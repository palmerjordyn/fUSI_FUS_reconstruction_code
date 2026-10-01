#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from pathlib import Path


CODE_DIR = Path(__file__).resolve().parent


SCRIPTS = [
    "parse_hobbes_design.py",
    "reconstruct.py",
    "analyze_hobbes.py",
    "make_biomarker.py",
    "make_anatomical_lgn_movie.py",
]


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--data-path",
        type=str,
        default=None
    )

    parser.add_argument(
        "--acq",
        type=str,
        default=None
    )

    parser.add_argument(
        "--date",
        type=str,
        default=None
    )

    parser.add_argument(
        "--rhs",
        type=str,
        default=None
    )

    parser.add_argument(
        "--reconstruction",
        type=str,
        default=None
    )

    parser.add_argument(
        "--verbose",
        action="store_true"
    )

    parser.add_argument(
        "--extra-identifier",
        type=str,
        default=None
    )

    return parser.parse_args()


def main():

    args = parse_arguments()


    # ============================================================
    # Set environment variables if MATLAB/command line supplied
    # them.
    #
    # If nothing was supplied, config.py uses its defaults.
    # ============================================================

    env = os.environ.copy()

    if args.data_path is not None:

        data_path = Path(args.data_path)

        if not data_path.exists():
            raise FileNotFoundError(
                f"Data path does not exist:\n{data_path}"
            )

        env["FUS_ZIP_PATH"] = str(data_path)
        env["FUS_DATA_PATH"] = str(data_path)

        if data_path.suffix.lower() == ".zip":
            env["FUS_ZIP"] = "true"
        else:
            env["FUS_ZIP"] = "false"


    if args.acq is not None:
        env["FUS_ACQ"] = args.acq

    if args.date is not None:
        env["FUS_DATE"] = args.date


    if args.rhs is not None:
        env["FUS_RHS"] = args.rhs


    if args.reconstruction is not None:
        env["FUS_RECONSTRUCTION"] = args.reconstruction


    if args.verbose:
        env["FUS_VERBOSE"] = "true"


    if args.extra_identifier is not None:
        env["FUS_EXTRA_IDENTIFIER"] = args.extra_identifier


    # ============================================================
    # Run pipeline
    # ============================================================

    print("=" * 60)
    print("Running fUSI/FUS analysis pipeline")
    print("=" * 60)


    for i, script_name in enumerate(SCRIPTS, start=1):

        script = CODE_DIR / script_name

        if not script.exists():
            raise FileNotFoundError(
                f"Script not found: {script}"
            )


        print()
        print("=" * 60)
        print(
            f"STEP {i}/{len(SCRIPTS)}: {script_name}"
        )
        print("=" * 60)
        print()


        subprocess.run(
            [sys.executable, str(script)],
            check=True,
            cwd=CODE_DIR,
            env=env
        )


        print()
        print(f"Finished: {script_name}")


    print()
    print("=" * 60)
    print("ALL STEPS COMPLETED SUCCESSFULLY")
    print("=" * 60)


if __name__ == "__main__":
    main()