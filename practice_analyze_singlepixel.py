
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz
from pathlib import Path
import zipfile


# --------------------------------------------------
# Path to ZIP file and BF file inside ZIP
# --------------------------------------------------
ZIP_FILE = Path(
    "/Users/jpalmer/Desktop/FUSI/HobbesfusiFUSAug312026.zip"
)

BF_FILE = (
    "HobbesfusiFUSAug312026/"
    "Acq_2026_08_31_14_46_09/"
    "P7-4__BF_020.bin"
)


# --------------------------------------------------
# Expected dimensions
# --------------------------------------------------
NZ = 339
NX = 38
NT = 200

X = 3
Z = 150

X = 24
Z = 196
# --------------------------------------------------
# Read the binary file directly from the ZIP
# --------------------------------------------------
with zipfile.ZipFile(ZIP_FILE, "r") as z:

    if BF_FILE not in z.namelist():
        raise FileNotFoundError(
            f"Could not find {BF_FILE} in {ZIP_FILE}"
        )

    raw = z.read(BF_FILE)


# Interpret the bytes as little-endian complex128
data = np.frombuffer(raw, dtype="<c16")


# Check that the file has the expected number of samples
expected = NZ * NX * NT

if data.size != expected:
    raise RuntimeError(
        f"Expected {expected} complex samples, "
        f"but found {data.size}"
    )


# Reshape into Z × X × time
cube = data.reshape(
    (NZ, NX, NT),
    order="F"
)


print(f"Loaded: {BF_FILE}")
print(f"Cube shape: {cube.shape}")


# --------------------------------------------------
# Parameters
# --------------------------------------------------
fs = 250       # Sampling frequency [Hz]
N = 200         # Number of IQ samples
cutoff = 70     # High-pass cutoff [Hz]
order = 4


# --------------------------------------------------
# Get pixel
# --------------------------------------------------
pixel = cube[Z, X, :].copy()




pixel = cube[Z, X, :] - cube[Z, X, 0]


# --------------------------------------------------
# FFT and manually keep frequencies < -70 Hz
# --------------------------------------------------
F = np.fft.fft(pixel)
freq = np.fft.fftfreq(len(pixel), d=1/fs)


# --------------------------------------------------
# Subtract first sample
# --------------------------------------------------
pixel = pixel - pixel[0]


# --------------------------------------------------
# 1. FFT after subtraction, before filtering
# --------------------------------------------------
F_before = np.fft.fft(pixel)

# Frequency axis corresponding to the FFT
freq = np.fft.fftfreq(N, d=1/fs)

# Shift zero frequency to the center
freq_shifted = np.fft.fftshift(freq)
F_before_shifted = np.fft.fftshift(F_before)


# --------------------------------------------------
# 2. Butterworth high-pass filter
# --------------------------------------------------
B, A = butter(
    order,
    cutoff / (fs / 2),
    btype="high"
)

print(f"Filter coefficients B: {B}")
print(f"Filter coefficients A: {A}")


# Filter the pixel
sb = lfilter(B, A, pixel)


# --------------------------------------------------
# 3. FFT after filtering
# --------------------------------------------------
F_after = np.fft.fft(sb)

F_after_shifted = np.fft.fftshift(F_after)


# --------------------------------------------------
# Filter frequency response
# --------------------------------------------------
# Evaluate H(f) at the SAME frequencies as the FFT
w = 2 * np.pi * freq_shifted / fs

_, H = freqz(B, A, worN=w)

H_shifted = np.abs(H)


# --------------------------------------------------
# Plot 1
# --------------------------------------------------
plt.figure(figsize=(10, 5))

plt.plot(
    freq_shifted,
    np.abs(F_before_shifted)
)

plt.xlabel("Frequency (Hz)")
plt.ylabel("FFT Magnitude")
plt.title(
    "1. Fourier Transform After Subtraction — Before Filtering"
)
plt.xlim(-125, 125)
plt.grid(True)
plt.tight_layout()


# --------------------------------------------------
# Plot 2
# --------------------------------------------------
plt.figure(figsize=(10, 5))

plt.plot(
    freq_shifted,
    H_shifted
)

plt.xlabel("Frequency (Hz)")
plt.ylabel("|H(f)|")
plt.title(
    "2. Butterworth High-Pass Filter Magnitude Response"
)
plt.xlim(-125, 125)
plt.ylim(0, 1.05)
plt.grid(True)
plt.tight_layout()


# --------------------------------------------------
# Plot 3
# --------------------------------------------------
plt.figure(figsize=(10, 5))

plt.plot(
    freq_shifted,
    np.abs(F_after_shifted)
)

plt.xlabel("Frequency (Hz)")
plt.ylabel("FFT Magnitude")
plt.title(
    "3. Fourier Transform After Filtering"
)
plt.xlim(-125, 125)
plt.grid(True)
plt.tight_layout()


#pixel = cube[Z, X, :].copy()

##PLOT 4##
plt.figure(figsize=(10, 5))

plt.plot(
    np.arange(N),
    np.abs(pixel)
)

plt.xlabel("IQ Sample")
plt.ylabel("Magnitude")
plt.title(f"Pixel Signal — Z={Z}, X={X}")
plt.grid(True)
plt.tight_layout()
plt.show()