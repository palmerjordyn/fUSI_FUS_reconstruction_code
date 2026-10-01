import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
from pathlib import Path
import zipfile


# ==================================================
# DATA SOURCE
#
# Set DATA_SOURCE to either:
#
#   1. A ZIP file
#   2. A folder containing the extracted data
# ==================================================

DATA_SOURCE = Path(
    "/Users/jpalmer/Desktop/FUSI/HobbesfusiFUSAug312026.zip"
)

# Example if using an extracted folder instead:
#
DATA_SOURCE = Path(
    "/Users/jpalmer/Downloads/HobbesFusiFUSAug312026"
)


# ==================================================
# ACQUISITION AND BF NUMBER
#
# Enter the acquisition and BF number here.
#
# The code will search the ENTIRE ZIP/folder for
# the first matching BF file.
# ==================================================

ACQUISITION = "Acq_2026_08_31_14_46_09"

BF_NUMBER = 305


# ==================================================
# Expected BF filename
# ==================================================

BF_FILENAME = f"P7-4__BF_{BF_NUMBER:03d}.bin"


# ==================================================
# Find BF file
#
# Searches recursively through the entire ZIP or
# folder. The exact surrounding folder structure
# does not matter.
# ==================================================

def find_bf_file(data_source, acquisition, bf_filename):

    # --------------------------------------------------
    # ZIP FILE
    # --------------------------------------------------

    if data_source.is_file() and data_source.suffix.lower() == ".zip":

        print(f"Reading from ZIP: {data_source}")
        print(f"Searching for:")
        print(f"  Acquisition: {acquisition}")
        print(f"  BF file:     {bf_filename}")

        with zipfile.ZipFile(data_source, "r") as z:

            matches = []

            for name in z.namelist():

                # Normalize path separators
                normalized_name = name.replace("\\", "/")

                # Check that this is a file with the desired
                # acquisition somewhere in its path and the
                # desired BF filename at the end.
                if (
                    acquisition in normalized_name
                    and normalized_name.endswith("/" + bf_filename)
                ):
                    matches.append(name)

            if len(matches) == 0:
                raise FileNotFoundError(
                    f"\nCould not find BF file.\n"
                    f"Acquisition: {acquisition}\n"
                    f"BF file:     {bf_filename}\n"
                    f"ZIP:         {data_source}"
                )

            # Use the first match
            bf_file_zip = matches[0]

            print("\nFound BF file:")
            print(f"  {bf_file_zip}")

            if len(matches) > 1:
                print(
                    f"\nWarning: found {len(matches)} matching files."
                )
                print("Using the first match.")

            raw = z.read(bf_file_zip)

            return raw, bf_file_zip


    # --------------------------------------------------
    # EXTRACTED FOLDER
    # --------------------------------------------------

    elif data_source.is_dir():

        print(f"Reading from folder: {data_source}")
        print(f"Searching for:")
        print(f"  Acquisition: {acquisition}")
        print(f"  BF file:     {bf_filename}")

        matches = []

        # Recursive search through everything
        # underneath DATA_SOURCE
        for path in data_source.rglob(bf_filename):

            # Make sure the acquisition appears somewhere
            # in the path
            if acquisition in str(path):

                matches.append(path)

        if len(matches) == 0:
            raise FileNotFoundError(
                f"\nCould not find BF file.\n"
                f"Acquisition: {acquisition}\n"
                f"BF file:     {bf_filename}\n"
                f"Folder:      {data_source}"
            )

        # Use the first match
        bf_path = matches[0]

        print("\nFound BF file:")
        print(f"  {bf_path}")

        if len(matches) > 1:
            print(
                f"\nWarning: found {len(matches)} matching files."
            )
            print("Using the first match.")

        raw = bf_path.read_bytes()

        return raw, bf_path


    # --------------------------------------------------
    # INVALID DATA SOURCE
    # --------------------------------------------------

    else:

        raise FileNotFoundError(
            f"DATA_SOURCE does not exist or is not a ZIP/folder:\n"
            f"{data_source}"
        )


# ==================================================
# Find and read requested BF file
# ==================================================

raw, found_file = find_bf_file(
    DATA_SOURCE,
    ACQUISITION,
    BF_FILENAME
)


# ==================================================
# Convert binary data to complex IQ
# ==================================================

data = np.frombuffer(
    raw,
    dtype="<c16"
)


# ==================================================
# Expected dimensions
# ==================================================

NZ = 339
NX = 38
NT = 200

expected = NZ * NX * NT


if data.size != expected:

    raise RuntimeError(
        f"Expected {expected} complex samples, "
        f"but found {data.size}"
    )


# ==================================================
# Reshape into Z × X × time
# ==================================================

cube = data.reshape(
    (NZ, NX, NT),
    order="F"
)

print("Loaded cube:", cube.shape)


# ==================================================
# Butterworth high-pass filter parameters
# ==================================================

fs = 250       # Sampling frequency [Hz]
cutoff = 30    # High-pass cutoff [Hz]
order = 4


# ==================================================
# Create Butterworth filter
# ==================================================

B, A = butter(
    order,
    cutoff / (fs / 2),
    btype="high"
)

print("Butterworth filter coefficients:")
print("B =", B)
print("A =", A)


# ==================================================
# Subtract first frame from all frames
# ==================================================

cube = cube - cube[:, :, 0][:, :, np.newaxis]


# ==================================================
# Apply Butterworth filter to every pixel
#
# axis=2 is the temporal/IQ dimension
# ==================================================

filtered = lfilter(
    B,
    A,
    cube,
    axis=2
)

print("Filtered cube:", filtered.shape)


# ==================================================
# Generate power-Doppler images
# ==================================================

WINDOW = 25

N_IMAGES = NT // WINDOW

power_doppler = np.empty(
    (N_IMAGES, NZ, NX),
    dtype=np.float32
)


for i in range(N_IMAGES):

    start = i * WINDOW
    end = (i + 1) * WINDOW

    q = filtered[:, :, start:end]

    # Sum squared magnitude over the 25 IQ frames
    power_doppler[i] = np.sum(
        np.abs(q) ** 2,
        axis=2
    )


print(
    "Power Doppler shape:",
    power_doppler.shape
)


# ==================================================
# Convert power Doppler to dB
# ==================================================

power_doppler_safe = np.maximum(
    power_doppler,
    1e-12
)

power_db = 10 * np.log10(
    power_doppler_safe / np.max(power_doppler_safe)
)


# ==================================================
# Display all 8 PD images
# ==================================================

fig, axes = plt.subplots(
    1,
    8,
    figsize=(14, 7),
    constrained_layout=True
)


# Same dB scale for all 8 images
vmin = -40
vmax = 0


for i, ax in enumerate(axes.flat):

    im = ax.imshow(
        power_db[i],
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
        aspect="equal"
    )

    ax.set_title(
        f"PD {i + 1}  "
        f"(IQ frames {i * 25 + 1}–{(i + 1) * 25})"
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Z")


# ==================================================
# Colorbar
# ==================================================

fig.colorbar(
    im,
    ax=axes,
    label="Power (dB)"
)


fig.suptitle(
    f"{BF_FILENAME} — Butterworth High-Pass Power Doppler\n"
    f"{ACQUISITION}",
    fontsize=16
)


# ==================================================
# Show all 8 PD images
# ==================================================

plt.show(block=False)


# ==================================================
# ROI ANALYSIS
# ==================================================

print("\n========================================")
print("ROI ANALYSIS")
print("========================================")


# --------------------------------------------------
# Ask for ROI radius
# --------------------------------------------------

radius = float(
    input("\nEnter ROI radius in pixels: ")
)

if radius <= 0:

    raise ValueError(
        "ROI radius must be greater than zero."
    )


# ==================================================
# Open second window showing PD image 1
# ==================================================

fig_roi, ax_roi = plt.subplots(
    figsize=(7, 7)
)

ax_roi.imshow(
    power_db[0],
    cmap="gray",
    vmin=vmin,
    vmax=vmax,
    aspect="equal"
)

ax_roi.set_xlabel("X")
ax_roi.set_ylabel("Z")

ax_roi.set_title(
    "PD 1 — Click center of SIGNAL ROI"
)

plt.show(block=False)


# ==================================================
# Select SIGNAL ROI
# ==================================================

print("\nClick the center of the SIGNAL ROI.")

signal_point = plt.ginput(
    1,
    timeout=-1
)

if len(signal_point) == 0:

    raise RuntimeError(
        "No signal ROI was selected."
    )

signal_x, signal_z = signal_point[0]

print(
    f"Signal ROI center: "
    f"X = {signal_x:.2f}, "
    f"Z = {signal_z:.2f}"
)


# ==================================================
# Draw SIGNAL ROI
# ==================================================

signal_circle = plt.Circle(
    (signal_x, signal_z),
    radius,
    fill=False,
    linewidth=2
)

ax_roi.add_patch(signal_circle)

fig_roi.canvas.draw()


# ==================================================
# Select BACKGROUND ROI
# ==================================================

ax_roi.set_title(
    "PD 1 — Click center of BACKGROUND ROI"
)

fig_roi.canvas.draw()

print("\nClick the center of the BACKGROUND ROI.")

background_point = plt.ginput(
    1,
    timeout=-1
)

if len(background_point) == 0:

    raise RuntimeError(
        "No background ROI was selected."
    )

background_x, background_z = background_point[0]

print(
    f"Background ROI center: "
    f"X = {background_x:.2f}, "
    f"Z = {background_z:.2f}"
)


# ==================================================
# Draw BACKGROUND ROI
# ==================================================

background_circle = plt.Circle(
    (background_x, background_z),
    radius,
    fill=False,
    linewidth=2
)

ax_roi.add_patch(background_circle)

ax_roi.set_title(
    "PD 1 — Signal and Background ROIs"
)

fig_roi.canvas.draw()


# ==================================================
# Create circular ROI masks
# ==================================================

z_grid, x_grid = np.ogrid[
    :NZ,
    :NX
]


signal_mask = (
    (x_grid - signal_x) ** 2
    +
    (z_grid - signal_z) ** 2
    <= radius ** 2
)


background_mask = (
    (x_grid - background_x) ** 2
    +
    (z_grid - background_z) ** 2
    <= radius ** 2
)


# ==================================================
# Extract signal and background values
# ==================================================

# Use LINEAR power values for the statistics,
# not the dB-converted image.

signal_values = power_doppler[0][signal_mask]

background_values = power_doppler[0][background_mask]


# ==================================================
# Compute ROI statistics
# ==================================================

usignal = np.mean(signal_values)
stdsignal = np.std(signal_values)

ubackground = np.mean(background_values)
stdb = np.std(background_values)


# ==================================================
# Compute SNR
# ==================================================

SNRdB = 20 * np.log10(
    usignal / stdb
)


# ==================================================
# Compute CNR
# ==================================================

CNRdB = 20 * np.log10(
    np.abs(usignal - ubackground)
    /
    np.sqrt(
        (stdsignal ** 2 + stdb ** 2) / 2
    )
)


# ==================================================
# Print ROI results
# ==================================================

print("\n========================================")
print("ROI RESULTS")
print("========================================")

print("\nSignal ROI")
print("----------------------------------------")
print(f"Center X:      {signal_x:.2f}")
print(f"Center Z:      {signal_z:.2f}")
print(f"Radius:        {radius:.2f} pixels")
print(f"Pixels:        {signal_values.size}")
print(f"Mean:          {usignal:.6g}")
print(f"Std:           {stdsignal:.6g}")

print("\nBackground ROI")
print("----------------------------------------")
print(f"Center X:      {background_x:.2f}")
print(f"Center Z:      {background_z:.2f}")
print(f"Radius:        {radius:.2f} pixels")
print(f"Pixels:        {background_values.size}")
print(f"Mean:          {ubackground:.6g}")
print(f"Std:           {stdb:.6g}")

print("\nImage quality")
print("----------------------------------------")
print(f"SNR:           {SNRdB:.3f} dB")
print(f"CNR:           {CNRdB:.3f} dB")

print("========================================")


# ==================================================
# Keep both windows open
# ==================================================

plt.show()