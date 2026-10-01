"""
This script contains all configurations for the run used by other scripts. It should be edited 
before run_all.py is ran.
Follow instructions below, you should only have to edit the DEFAULT USER SETTINGS section. 
"""

import os
from pathlib import Path
import zipfile


# ================================================================
# DEFAULT USER SETTINGS
# These are used when run_all.py is run directly.
# ================================================================


DEFAULT_DATE = "20260831"

#If using a Zip file, set this to True and enter path in DEFAULT_ZIP_PATH. 
#If using an extracted folder, set this to False and enter path in DEFAULT_DATA_PATH. 
# You can leave the one you don't use empty.
DEFAULT_ZIP = True

DEFAULT_ZIP_PATH = Path(
    "/Users/jpalmer/Downloads/HobbesfusiFUSAug312026.zip"
)

DEFAULT_DATA_PATH = Path(
    "/Users/jpalmer/Downloads/HobbesfusiFUSAug312026"
)

#Enter the acquisition names you want to process. These should match the folder names in the zip or extracted folder.
#Ex. "Acq_2026_09_16_10_57_16" can be multiple separated by commas. Enter Acqs first to last of when they were taken. 
# DON'T MIX ACQ FROM DIFFERENT EXPERIMENT TYPES. ex. don't mix theta experiment acqs with normal 4 trigger experiments.
# Run the code twice once for each experiment type, see DEFAULT_EXTRA_IDENTIFIER
DEFAULT_ACQS = [
    "Acq_2026_08_31_15_07_01",
    "Acq_2026_08_31_14_58_14",
    "Acq_2026_08_31_14_46_09"
]


# The acquisitions have to be matched up to their associated RHS files. 
# In the same order as the Acq, find an identifier for the rhs files associated with that acq.
# DEFAULT_ACQS[0] matches up with DEFAULT_RHS_IDENTFIER[0] so make sure they are in order.
# The code will find all files with DEFUALT_RHS_IDENTIFIER[0] and .rhs and associate them with DEFAULT_ACQS[0]. 
DEFAULT_RHS_IDENTIFIER = [
    "take3",
    "take2",
    "hobbesFusiFus20sblock15rep_2"
]

# There are currently two options for reconstruction type "BUTTER" or "SVD". "BUTTER" performs a 30 Hz high pass 
# Butterworth filter to isolate the blood signal. "SVD" uses a SVD math model to isolate the blood signal. 
# You can run both separately to compare without anything getting overwritten. 
# The end results will go in different folders according to the filter you used. 
DEFAULT_RECONSTRUCTION_TYPE = "BUTTER"

# If DEFAULT_VERBOSE is true, extra graphs and data will appear to help trouble shoot. 
DEFAULT_VERBOSE = False

# This can be used if you want to rerun the code and change something but don't want any resultant files overwritten.
# ex. if you have theta runs and normal runs on the same day, you might add DEFAULT_EXTRA_IDENTIFIER = "theta" when you 
# run the theta acquisitions so that a new folder gets created and no data gets overwritten. If you are not doing
# multiple runs per date, this can stay = ""
DEFAULT_EXTRA_IDENTIFIER = ""


# ================================================================
# USE ENVIRONMENT OVERRIDES IF THEY EXIST
# ================================================================

DATE = os.environ.get(
    "FUS_DATE",
    DEFAULT_DATE
)

ZIP = os.environ.get(
    "FUS_ZIP",
    str(DEFAULT_ZIP)
).lower() == "true"


ZIP_PATH = Path(
    os.environ.get(
        "FUS_ZIP_PATH",
        str(DEFAULT_ZIP_PATH)
    )
)


DATA_PATH = Path(
    os.environ.get(
        "FUS_DATA_PATH",
        str(DEFAULT_DATA_PATH)
    )
)


_acq = os.environ.get("FUS_ACQ")

if _acq:
    ACQS = [_acq]
else:
    ACQS = DEFAULT_ACQS


_rhs = os.environ.get("FUS_RHS")

if _rhs:
    RHS_IDENTIFIER = [_rhs]
else:
    RHS_IDENTIFIER = DEFAULT_RHS_IDENTIFIER


RECONSTRUCTION_TYPE = os.environ.get(
    "FUS_RECONSTRUCTION",
    DEFAULT_RECONSTRUCTION_TYPE
)


VERBOSE = os.environ.get(
    "FUS_VERBOSE",
    str(DEFAULT_VERBOSE)
).lower() == "true"


EXTRA_IDENTIFIER = os.environ.get(
    "FUS_EXTRA_IDENTIFIER",
    DEFAULT_EXTRA_IDENTIFIER
)


# ================================================================
# AUTOMATIC PATH SETUP
# ================================================================

CODE_DIR = Path(__file__).resolve().parent

BASE_DIR = CODE_DIR.parent

OUTPUT_DIR = (
    BASE_DIR
    / "fUSI_FUS_output"
    / f"hobbes{DATE}{EXTRA_IDENTIFIER}"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ================================================================
# COUNT BF FILES
# ================================================================

BF = []

if ZIP:

    with zipfile.ZipFile(ZIP_PATH) as z:

        for acq in ACQS:

            num_bf = sum(
                1
                for n in z.namelist()
                if acq in n
                and Path(n).name.startswith("P7-4__BF_")
                and n.endswith(".bin")
            )

            BF.append(num_bf)

else:

    for acq in ACQS:

        acq_dir = DATA_PATH / acq

        num_bf = len([
            p for p in acq_dir.iterdir()
            if p.is_file()
            and p.name.startswith("P7-4__BF_")
            and p.suffix == ".bin"
        ])

        BF.append(num_bf)