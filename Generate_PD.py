import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
from pathlib import Path
import zipfile
from matplotlib.widgets import Slider


# ==================================================
# DATA SOURCE
#
# Set DATA_SOURCE to either:
#
#   1. A ZIP file
#   2. A folder containing the extracted data
# ==================================================

# DATA_SOURCE = Path(
#     "/Users/jpalmer/Desktop/FUSI/HobbesfusiFUSAug312026.zip"
# )
# Example if using an extracted folder instead:
#
DATA_SOURCE = Path(
    "/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected"
)


# ==================================================
# ACQUISITION AND BF NUMBER
#
# Enter the acquisition and BF number here.
#
# The code will search the ENTIRE ZIP/folder for
# the first matching BF file.
# ==================================================

ACQUISITION = "Acq_2026_10_05_10_58_30"

BF_NUMBER = 31


# ==================================================
# Do you want to Save the images?
# If true, fill out date
# ==================================================
Save = True
DATE = "20261005"

# ==================================================
# Expected BF filename
# ==================================================

BF_FILENAME = f"P7-4__BF_{BF_NUMBER:03d}.bin"


# ==================================================
# OUTPUT DIRECTORY
# ==================================================

CODE_DIR = Path(__file__).resolve().parent
BASE_DIR = CODE_DIR.parent

OUTPUT_DIR = (
    BASE_DIR
    / "fUSI_FUS_output"
    / "PDs"
)

if Save:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ==================================================
# Find BF file
#
# Searches recursively through the entire ZIP or
# folder. The exact surrounding folder structure
# does not matter.
# ==================================================
# ==================================================
# FIND BF FILE
# ==================================================

def find_bf_file(data_source, acquisition, bf_filename):

    matches = []

    # --------------------------------------------------
    # ZIP FILE
    # --------------------------------------------------

    if data_source.is_file() and data_source.suffix.lower() == ".zip":

        print(f"Reading from ZIP: {data_source}")
        print(f"Searching acquisition: {acquisition}")
        print(f"Searching for BF file: {bf_filename}")

        with zipfile.ZipFile(data_source, "r") as z:

            for name in z.namelist():

                normalized = name.replace("\\", "/")
                parts = normalized.strip("/").split("/")

                # Find the acquisition folder anywhere in the ZIP.
                if acquisition not in parts:
                    continue

                # Match the BF filename at any nesting depth.
                if parts[-1] == bf_filename:
                    matches.append(name)

            if not matches:
                raise FileNotFoundError(
                    f"\nCould not find BF file.\n"
                    f"Acquisition: {acquisition}\n"
                    f"BF file:     {bf_filename}\n"
                    f"ZIP:         {data_source}"
                )

            matches.sort()

            print("\nFound BF file:")
            print(f"  {matches[0]}")

            if len(matches) > 1:
                print(
                    f"\nWarning: found {len(matches)} matching files. "
                    "Using the first match."
                )

            raw = z.read(matches[0])

            return raw, matches[0]

    # --------------------------------------------------
    # EXTRACTED FOLDER
    # --------------------------------------------------

    elif data_source.is_dir():

        print(f"Reading from folder: {data_source}")
        print(f"Searching acquisition: {acquisition}")
        print(f"Searching for BF file: {bf_filename}")

        # Locate acquisition directories anywhere under DATA_SOURCE.
        acquisition_dirs = [
            path for path in data_source.rglob(acquisition)
            if path.is_dir()
        ]

        # Also allow DATA_SOURCE itself to be the acquisition folder.
        if data_source.name == acquisition:
            acquisition_dirs.insert(0, data_source)

        if not acquisition_dirs:
            raise FileNotFoundError(
                f"\nCould not find acquisition folder.\n"
                f"Acquisition: {acquisition}\n"
                f"Folder:      {data_source}"
            )

        # Search recursively within every matching acquisition folder.
        for acquisition_dir in acquisition_dirs:

            for path in acquisition_dir.rglob(bf_filename):
                if path.is_file():
                    matches.append(path)

        if not matches:
            raise FileNotFoundError(
                f"\nCould not find BF file inside the acquisition folder.\n"
                f"Acquisition: {acquisition}\n"
                f"BF file:     {bf_filename}\n"
                f"Folder:      {data_source}"
            )

        matches.sort(key=lambda path: str(path))

        print("\nFound BF file:")
        print(f"  {matches[0]}")

        if len(matches) > 1:
            print(
                f"\nWarning: found {len(matches)} matching files. "
                "Using the first match."
            )

        raw = matches[0].read_bytes()

        return raw, matches[0]

    # --------------------------------------------------
    # INVALID DATA SOURCE
    # --------------------------------------------------

    else:

        raise FileNotFoundError(
            f"DATA_SOURCE does not exist or is not a ZIP/folder:\n"
            f"{data_source}"
        )



def select_point_on_axes(fig, ax, prompt):
    """Select a point only when clicking inside the image axes."""
    selected_point = []

    print(prompt)

    def on_click(event):
        # Ignore slider clicks and clicks outside the image.
        if event.inaxes is not ax:
            return

        if event.button == 1 and event.xdata is not None and event.ydata is not None:
            selected_point.append((event.xdata, event.ydata))

    callback_id = fig.canvas.mpl_connect("button_press_event", on_click)

    try:
        while not selected_point and plt.fignum_exists(fig.number):
            plt.pause(0.05)
    finally:
        fig.canvas.mpl_disconnect(callback_id)

    if not selected_point:
        raise RuntimeError("Figure closed before a point was selected.")

    return selected_point[0]


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
# Display all 8 PD images with adjustable dB scale
# ==================================================

fig, axes = plt.subplots(
    1,
    8,
    figsize=(14, 7)
)

# Leave room at the bottom for the slider
fig.subplots_adjust(
    bottom=0.22,
    top=0.82,
    wspace=0.35,
    right=0.92
)

# Initial brightness limits
vmin = -40
vmax = 0

# Display all images using the same dB scale
images = []

for i, ax in enumerate(axes.flat):

    im = ax.imshow(
        power_db[i],
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
        aspect="equal"
    )

    images.append(im)

    ax.set_title(
        f"PD {i + 1} "
        f"(IQ frames {i * 25 + 1}–{(i + 1) * 25})"
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Z")

# Add a colorbar
cbar = fig.colorbar(
    images[-1],
    ax=axes,
    label="Power (dB)"
)

# Add the brightness slider along the bottom of the figure
slider_ax = fig.add_axes([0.20, 0.08, 0.60, 0.035])

brightness_slider = Slider(
    ax=slider_ax,
    label="Lower dB bound",
    valmin=-60,
    valmax=-10,
    valinit=-40,
    valstep=1
)

# Update all images when the slider moves
def update_brightness(value):

    for im in images:
        im.set_clim(
            vmin=value,
            vmax=0
        )

    fig.canvas.draw_idle()

brightness_slider.on_changed(update_brightness)

# Figure title
fig.suptitle(
    f"{BF_FILENAME} — Butterworth High-Pass Power Doppler\n"
    f"{ACQUISITION}",
    fontsize=16
)

# Display the images without blocking the ROI workflow
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
# Individual PD image with its own brightness slider
# ==================================================

fig_roi, ax_roi = plt.subplots(figsize=(7, 7))

# Leave room below the image for the slider.
fig_roi.subplots_adjust(bottom=0.20)

roi_image_index = 2  # 0 = PD 1, 1 = PD 2, etc.

roi_image = ax_roi.imshow(
    power_db[roi_image_index],
    cmap="gray",
    vmin=vmin,
    vmax=0,
    aspect="equal"
)

ax_roi.set_xlabel("X")
ax_roi.set_ylabel("Z")
ax_roi.set_title(
    f"PD {roi_image_index + 1} — Click center of SIGNAL ROI"
)

# Slider axes
roi_slider_ax = fig_roi.add_axes([0.20, 0.07, 0.65, 0.035])

roi_brightness_slider = Slider(
    ax=roi_slider_ax,
    label="Lower dB bound",
    valmin=-60,
    valmax=-10,
    valinit=vmin,
    valstep=1
)

def update_roi_brightness(value):
    roi_image.set_clim(vmin=value, vmax=0)
    fig_roi.canvas.draw_idle()

roi_brightness_slider.on_changed(update_roi_brightness)

plt.show(block=False)



# ==================================================
# Helper: Select a point only inside the image axes
# ==================================================

def select_point_on_axes(fig, ax, prompt):
    selected_point = []

    print(prompt)

    def on_click(event):
        # Ignore clicks on the slider or anywhere outside the image.
        if event.inaxes is not ax:
            return

        # Only accept left-clicks with valid image coordinates.
        if (
            event.button == 1
            and event.xdata is not None
            and event.ydata is not None
        ):
            selected_point.append((event.xdata, event.ydata))

    callback_id = fig.canvas.mpl_connect(
        "button_press_event",
        on_click
    )

    try:
        while not selected_point and plt.fignum_exists(fig.number):
            plt.pause(0.05)
    finally:
        fig.canvas.mpl_disconnect(callback_id)

    if not selected_point:
        raise RuntimeError(
            "Figure was closed before an ROI point was selected."
        )

    return selected_point[0]


# ==================================================
# Select SIGNAL ROI
# ==================================================

signal_x, signal_z = select_point_on_axes(
    fig_roi,
    ax_roi,
    "\nClick the center of the SIGNAL ROI inside the image."
)

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
fig_roi.canvas.draw_idle()


# ==================================================
# Select BACKGROUND ROI
# ==================================================

ax_roi.set_title(
    f"PD {roi_image_index + 1} — Click center of BACKGROUND ROI"
)

fig_roi.canvas.draw_idle()

background_x, background_z = select_point_on_axes(
    fig_roi,
    ax_roi,
    "\nClick the center of the BACKGROUND ROI inside the image."
)

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
    f"PD {roi_image_index + 1} — Signal and Background ROIs"
)

fig_roi.canvas.draw_idle()


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



# ==================================================
# SAVE INDIVIDUAL PD IMAGE
# ==================================================

if Save:

    # Save the selected individual PD image at -60 to 0 dB.
    save_fig, save_ax = plt.subplots(figsize=(7, 7))

    save_ax.imshow(
        power_db[roi_image_index],
        cmap="gray",
        vmin=-60,
        vmax=0,
        aspect="equal"
    )

    # # Draw the signal and background ROIs.
    # save_ax.add_patch(
    #     plt.Circle(
    #         (signal_x, signal_z),
    #         radius,
    #         fill=False,
    #         color="yellow",
    #         linewidth=2,
    #         label="Signal ROI"
    #     )
    # )

    # save_ax.add_patch(
    #     plt.Circle(
    #         (background_x, background_z),
    #         radius,
    #         fill=False,
    #         color="cyan",
    #         linewidth=2,
    #         label="Background ROI"
    #     )
    # )

    save_ax.set_xlabel("X")
    save_ax.set_ylabel("Z")

    # Include the ROI statistics on the saved image.
    save_ax.text(
        0.2,
        0.97,
        f"SNR: {SNRdB:.3f} dB\nCNR: {CNRdB:.3f} dB",
        transform=save_ax.transAxes,
        va="top",
        ha="left",
        color="white",
        fontsize=12,
        bbox=dict(
            facecolor="black",
            alpha=0.75,
            edgecolor="white"
        )
    )

    save_ax.set_title(
        f"PD {roi_image_index + 1} — {ACQUISITION}"
    )

    save_ax.legend(loc="lower right")

    save_fig.tight_layout()

    save_path = (
        OUTPUT_DIR
        / f"{DATE}_{ACQUISITION}.png"
    )

    save_fig.savefig(
        save_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(save_fig)

    print(f"\nSaved PD image: {save_path}")


plt.show()