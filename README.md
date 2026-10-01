# Hobbes September 10, 2026 Theta-Burst fUSI Analysis SOP

## Scope

This code reconstructs BF files from Hobbes into visual graphs/images. 

Two main options, see more details for each below:
1) Obtain "theta_burst_10hz_anatomical_lgn_response_movie.mp4" and associated figures/data: run_all.py
2) Generate example PD image to check acquisition validity: Generate_PD.py

## Optional Installation of Virtual Environment

Create and activate a virtual environment outside of folder created by GitHub (Requirements are fairly simple -see requirements.txt- so if they are already installed globally this isn't necessary):

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

## Option 1 steps: generating LGN fUSI video
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

## Option 2 steps: generating PD images and estimating SNR, CNR
1) If you are on MATLAB, you can run "butterworth_power_doppler_roi.m", if using python you can run Generate_PD.py
2) You only need to edit DATA_SOURCE (.zip or folder path), ACQUISITION, and BF_NUMBER
    - ex. 
        DATA_SOURCE = Path(
            "/Users/jpalmer/Downloads/HobbesFusiFUSAug312026"
        )

        ACQUISITION = "Acq_2026_08_31_14_46_09"

        BF_NUMBER = 305
3) This will first show you the 8 PD images the BF file makes up (the 200 images are split up into sections of 25, each forming one PD image). Checking to see if the nerve is faintly visible around (300, 27) indicates if the acquisition is valid before further reconstruction is done. 
4) Next you will be prompted in the terminal: "Enter ROI radius in pixels: ". Enter the radius of signal area you want to analyze (for the nerve 3 to 5 works). Then hit enter. An image should pop up of the first PD image. First select the center of your signal, then select an area of background noise. The code will report the SNR and CNR in dB in the terminal. A good run was found to have a SNR of 20dB and CNR of 5dB. 

## MORE INFORMATION: WHAT IS THIS CODE DOING?
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