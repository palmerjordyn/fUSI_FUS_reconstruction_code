# Hobbes September 10, 2026 Theta-Burst fUSI Analysis SOP

## Scope

This code reconstructs BF files from Hobbes into visual graphs/images. 

Two main options, see more details for each below:
1) Obtain "theta_burst_10hz_anatomical_lgn_response_movie.mp4" and associated figures/data: run_all.py
2) Generate example PD image to check acquisition validity, SNR and CNR: Generate_PD.py

## Software and Optional Installation of Virtual Environment

The following is necessary to run the code: 

- Python 3.10 or newer
- NumPy
- SciPy
- Matplotlib
- OpenCV (`opencv-python`)
- `imageio-ffmpeg`

Either 
1) Install missing packages with:

```bash
python3 -m pip install numpy scipy matplotlib opencv-python imageio-ffmpeg
```
    - I found that I had to 'brew install ffmpeg'
OR

2) Create and activate a virtual environment OUTSIDE of the folder created by GitHub (if in folder created by GitHub in terminal cd .. then run the following)

```bash
python -m venv fUSI_FUS_venv

source fUSI_FUS_venv/bin/activate
```

If on MAC, may need to delete extra files made when copying from GitHub. If you don't, installing requirements won't work.
In directory outside of both GitHub and environment:

```bash
find . -name '._*' -delete
```

Install dependencies from github folder:

```bash
cd fUSI_FUS_code
pip install -r requirements.txt
```
After venv is created, it can be activated in a new terminal using source fUSI_FUS_venv/bin/activate

## Code 1 steps: generating LGN fUSI video
1) edit config.py 
    - input date, data path, etc. Follow instructions in config.py
2) run run_all.py
    - This will run 
        "parse_hobbes_design.py",
        "reconstruct.py",
        "analyze_hobbes.py",
        "make_biomarker.py",
        "make_anatomical_lgn_movie.py",
3) That should be all you have to do! The code will run and show updates in the terminal. 
    All outputs will be sent to a folder "fUSI_FUS_output" which is in the same folder as fUSI_FUS_code. Reconstruction arrays will be under "reconstruction_10hz", videos will be under "ttl_aligned_10hz". Both these folders are under a folder titled BUTTER or SVD depending on the reconstruction type used. 
4) If you want to know more about what is happening in each script run, see MORE INFORMATION section of this README.

## Code 2 steps: generating PD images and estimating SNR, CNR
1) If you are on MATLAB, you can run "Generate_PD.m", if using python you can run Generate_PD.py
2) You only need to edit DATA_SOURCE (.zip or folder path), ACQUISITION, and BF_NUMBER
    - ex. 
        DATA_SOURCE = Path(
            "/Users/jpalmer/Downloads/HobbesFusiFUSAug312026"
        )

        ACQUISITION = "Acq_2026_08_31_14_46_09"

        BF_NUMBER = 305
3) This will first show you the 8 PD images the BF file makes up (the 200 images are split up into sections of 25, each forming one PD image). Checking to see if the nerve is faintly visible around (300, 27) indicates if the acquisition is valid before further reconstruction is done. 
4) Next you will be prompted in the terminal: "Enter ROI radius in pixels: ". Enter the radius of signal area you want to analyze (for the nerve 3 to 5 works). Then hit enter. An image should pop up of the first PD image. First select the center of your signal, then select an area of background noise. The code will report the SNR and CNR in dB in the terminal. A good run was found to have a SNR of 20dB and CNR of 5dB. 

## MORE INFORMATION: WHAT IS THIS CODE DOING? AND OTHER FILES
1) parse_hobbes_design.py : Parses RHS files
    - produces `design/design_summary.json` which summarizes the fUSI and FUS triggers
    - This code ensures that the .rhs files match up with the BF data. 

2) "reconstruct.py" : Reconstruct 10 Hz power Doppler

    For each parent BF file, the code:

    - loads one 339 × 38 × 200 complex-IQ cube;
    - applies SVD or Butterworth filtering to all 200 IQ frames;
    - SVD filtering removes singular components 0–29;
    - Butteroworth filter is 30Hz high pass. 
    - divides the filtered cube into eight consecutive groups of 25 IQ frames;
    - sums squared magnitude within each group to form eight power-Doppler images.

This produces 10 Hz sampling during the approximately 800 ms acquisition burst. It is not continuous 10 Hz acquisition: the inter-acquisition/save gap is preserved as missing data rather than interpolated.

3) "analyze_hobbes.py" : Runs the primary z=5 analysis

    - Artifact detection is performed independently within each run. For every 10 Hz image, the median log10 power Doppler over the supported whole fUSI field is calculated and robustly z-scored using the run median and MAD. A parent BF file is flagged when any of its eight subimages has absolute robust z greater than or equal to 5. If a flagged parent file occurs within a paired 20-second trial, the complete pair is rejected.
    - This is to take out the pairs that occurred when Hobbes moved.
    - This censor metric is a whole-image flash metric, not an LGN signal metric. The anatomical LGN occupies too few pixels to drive the field median.

4) "make_biomarker.py" and "make_anatomical_lgn_movie.py"
    - Uses files in fUSI_FUS_code/outputs to locate the lgn and create complete anatomical movie.
    - make_biomarker.py looks for significant differences between the Contralateral and ipsilateral responses. If it can not find significant differences it will throw a warning but continue to make the video nonetheless. 

5) "Generate_PD.py"
    - This code uses Butterworth filtering and then sums squared magnitude within each group to form eight power-Doppler images.

6) "practice_analyze_singlepixel.py"
    - This code helps you analyze what is happening to a single pixel in one BF file (across all 200 images). You can select the pixel and BF file desired (perhaps use Generate_PD.py to determine which pixel you want to analyze). Graphs will demonstrate pixel signal in the time domain and frequency domain, the filter used and the post filtered frequency domain of the pixel. 

7) "Generate_PD.m"
    - Same thing as Generate_PD.py but works in MATLAB 

8) "decode_itan_fusi_timing.py"
    - Usually not ran, contains functions used by parse to understand the .rhs files. 

9) "make_10hz_uncensored_movie.py"
    - Allows you to change Z value that removes artifacts. For uncensored use CENSOR_Z = float("nan"), for Z = 4, use CENSOR_Z = float(4)

10) "make_hobbes_design.py"
    - Creates pdf and png traces that could be used in reports. 

11) "plot_hobbes_anatomical_lgn_bar.py"
    - anatomical-LGN whole-period bar plot.

12) FUS_recon.m
    - this calls run_all.py from matlab. If avoidable, I would not suggest this route. Many issues can crop up due to mixed architecture. If necessary, here are instructions. 
    - make sure to create environment following steps in "Software and Optional Installation of Virtual Environment" option 2 at the top of this readme. The MATLAB code takes an environment venv folder. If you use a different name for the venv, you will have to change this in the FUS_recon.m file, find venv_folder. 
    - open FUS_recon.m in matlab, edit the user inputs like date, path to file, etc.
    - run the code. 
    - If you are having issues with x86_64 vs arm64 then your matlab architecture does not match your python version architecture. Follow the instructions below in MATLAB and Python Architecture. 


## MATLAB and Python Architecture

### Why this matters

On Apple Silicon Macs, Python and MATLAB may run under different CPU architectures.

* Apple Silicon Python normally runs as **ARM64**
* Intel MATLAB runs as **x86_64**
* Python packages such as NumPy and SciPy contain compiled code, so the Python environment must match the architecture of the application that is running it.

For example:

```text
ARM64 MATLAB  →  ARM64 Python environment
x86_64 MATLAB →  x86_64 Python environment
```

If the architectures do not match, errors such as the following may occur:

```text
mach-o file, but is an incompatible architecture
have 'arm64', need 'x86_64'
```

---

### 1. Check the MATLAB architecture

Open MATLAB and run:

```matlab
computer('arch')
```

You should see one of:

```text
maca64
```

or:

```text
maci64
```

These mean:

| MATLAB output | Architecture          |
| ------------- | --------------------- |
| `maca64`      | ARM64 / Apple Silicon |
| `maci64`      | x86_64 / Intel        |

You can also check from MATLAB with:

```matlab
system('uname -m')
```

Expected output:

```text
arm64
```

for ARM64 MATLAB, or:

```text
x86_64
```

for Intel MATLAB.

---

### 2. Check the Python virtual environment architecture

Activate the virtual environment in Terminal:

```bash
source fUSI_FUS_venv/bin/activate
```

Then run:

```bash
python -c "import platform, sys; print('Architecture:', platform.machine()); print('Python:', sys.executable)"
```

You should see:

```text
Architecture: arm64
```

or:

```text
Architecture: x86_64
```

Compare this with the MATLAB architecture.

For example:

```text
MATLAB:  maca64
Python:  arm64
```

is a match.

And:

```text
MATLAB:  maci64
Python:  x86_64
```

is also a match.

---

### 3. If MATLAB and Python have different architectures

If you have:

```text
MATLAB: maci64
Python: arm64
```

or:

```text
MATLAB: x86_64
Python: arm64
```

you should create a separate x86_64 Python environment for MATLAB.

Do **not** delete the existing ARM64 environment. The ARM64 environment can still be used by native Apple Silicon applications such as VS Code.

A useful setup on an Apple Silicon Mac is:

```text
fUSI_FUS_analysis/
├── fUSI_FUS_venv/          ← ARM64 environment
├── fUSI_FUS_venv_x86/      ← x86_64 environment for Intel MATLAB
└── fUSI_FUS_code/
```

---

### 4. Create the x86_64 virtual environment

First remove an incorrectly created x86 environment if necessary:

```bash
rm -rf fUSI_FUS_venv_x86
```

Use the universal Python installation explicitly and force it to run as x86_64:

```bash
arch -x86_64 /Library/Frameworks/Python.framework/Versions/3.13/bin/python3 \
    -m venv fUSI_FUS_venv_x86
```

Activate it:

```bash
source fUSI_FUS_venv_x86/bin/activate
```

The normal Terminal on an Apple Silicon Mac may still report `arm64` because the Terminal itself is running natively as ARM64. Therefore, explicitly test the Python executable as x86_64:

```bash
arch -x86_64 fUSI_FUS_venv_x86/bin/python \
    -c "import platform, sys; print('Architecture:', platform.machine()); print('Python:', sys.executable)"
```

The expected result is:

```text
Architecture: x86_64
Python: .../fUSI_FUS_venv_x86/bin/python
```

---

### 5. Verify the Python executable contains x86_64

You can also check the executable directly:

```bash
file fUSI_FUS_venv_x86/bin/python
```

On an Apple Silicon Mac using the universal Python installation, it may report:

```text
Mach-O universal binary with 2 architectures:
[x86_64:Mach-O 64-bit executable x86_64]
[arm64:Mach-O 64-bit executable arm64]
```

This is okay.

The important test is that it can be executed as x86_64:

```bash
arch -x86_64 fUSI_FUS_venv_x86/bin/python \
    -c "import platform; print(platform.machine())"
```

which should return:

```text
x86_64
```

---

### 6. Install the project requirements into the x86 environment

With the x86 environment activated, install the requirements using the x86 Python explicitly:

```bash
arch -x86_64 fUSI_FUS_venv_x86/bin/python \
    -m pip install --upgrade pip
```

Then:

```bash
arch -x86_64 fUSI_FUS_venv_x86/bin/python \
    -m pip install -r fUSI_FUS_code/requirements.txt
```

Verify NumPy and SciPy:

```bash
arch -x86_64 fUSI_FUS_venv_x86/bin/python \
    -c "import platform, numpy, scipy; print('Architecture:', platform.machine()); print('NumPy:', numpy.__version__); print('SciPy:', scipy.__version__)"
```

The expected output should include:

```text
Architecture: x86_64
```

---

### 7. Configure MATLAB to use the correct environment

The MATLAB entry script should point to the x86 environment when MATLAB is running as `maci64`.

For example:

```matlab
venv_folder = fullfile(project_folder, 'fUSI_FUS_venv_x86');

python = fullfile( ...
    venv_folder, ...
    'bin', ...
    'python3');
```

For ARM64 MATLAB, use the normal ARM environment:

```matlab
venv_folder = fullfile(project_folder, 'fUSI_FUS_venv');

python = fullfile( ...
    venv_folder, ...
    'bin', ...
    'python3');
```

---

### 8. Verify the Python used by MATLAB

The MATLAB script should check the Python architecture before running the analysis.

From MATLAB:

```matlab
python_info_command = sprintf( ...
    '"%s" -c "import sys, platform; print(sys.version.split()[0]); print(platform.machine()); print(sys.executable)"', ...
    python);

system(python_info_command);
```

For Intel MATLAB, the output should contain:

```text
x86_64
```

and the executable should point to:

```text
fUSI_FUS_venv_x86/bin/python3
```

For ARM64 MATLAB, the output should contain:

```text
arm64
```

and point to:

```text
fUSI_FUS_venv/bin/python3
```

---

### 9. Important: do not commit virtual environments to Git

Virtual environments contain machine-specific binaries and should not be included in the repository.

Do not commit:

```text
fUSI_FUS_venv/
fUSI_FUS_venv_x86/
```

Instead, commit:

```text
requirements.txt
README.md
```

Each user should create their own virtual environment appropriate for their operating system and MATLAB architecture.

For example:

```text
Apple Silicon + ARM64 MATLAB
    ↓
fUSI_FUS_venv

Apple Silicon + Intel MATLAB
    ↓
fUSI_FUS_venv_x86

Intel Mac + Intel MATLAB
    ↓
fUSI_FUS_venv
```

The exact environment name is not important. What matters is that the Python interpreter and MATLAB process use the same architecture.
