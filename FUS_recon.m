
%% run_fUSI_analysis.m
%
% MATLAB entry point for the fUSI/FUS analysis pipeline.
%
% This script automatically uses the Python virtual environment:
%
%     fUSI_FUS_venv/
%
% which should be located next to:
%
%     fUSI_FUS_code/
%
% Expected project structure:
%
%     fUSI_FUS_analysis/
%     |
%     +-- fUSI_FUS_code/
%     |   +-- run_fUSI_analysis.m
%     |   +-- run_all.py
%     |   +-- config.py
%     |   +-- ...
%     |
%     +-- fUSI_FUS_venv/
%
%
% The Python environment should be created using the instructions in
% README.md and the packages in requirements.txt.
%
% Python can still be run directly from VS Code with:
%
%     python run_all.py
%
% When this MATLAB script is run, the values below override the
% corresponding values in config.py.
%
% ========================================================================


clear;
clc;


%% ========================================================================
% USER INPUT
% ========================================================================


date = "20260916";

% -------------------------------------------------------------------------
% Path to the ZIP file OR unzipped data directory
% -------------------------------------------------------------------------

path_to_data = ...
    '/Users/jpalmer/Downloads/HobbesFusiFusSeptember16th2026.zip';


% -------------------------------------------------------------------------
% Acquisition to analyze
% -------------------------------------------------------------------------

ACQ_to_run = ...
    'Acq_2026_09_16_10_22_12';


% -------------------------------------------------------------------------
% RHS identifier
% -------------------------------------------------------------------------

RHS_identifier = ...
    '100ms_take1';


% -------------------------------------------------------------------------
% Reconstruction type
%
% Examples:
%     'BUTTER'
%     'SVD'
% -------------------------------------------------------------------------

reconstruction_type = ...
    'BUTTER';


% -------------------------------------------------------------------------
% Optional extra identifier for the output folder
% -------------------------------------------------------------------------

extra_identifier = ...
    'test';


%% ========================================================================
% LOCATE PROJECT
% ========================================================================

% Folder containing this MATLAB script.
%
% fileparts(mfilename('fullpath')) means the user does NOT need to
% hard-code the location of fUSI_FUS_code.

code_folder = fileparts(mfilename('fullpath'));


% The virtual environment is expected to be located next to the
% fUSI_FUS_code folder.

project_folder = fileparts(code_folder);


% Expected virtual environment folder

venv_folder = fullfile(project_folder, 'fUSI_FUS_venv');


%% ========================================================================
% LOCATE PYTHON
% ========================================================================

% On macOS/Linux, the executable should be:
%
%     fUSI_FUS_venv/bin/python
%
% On Windows, it should be:
%
%     fUSI_FUS_venv/Scripts/python.exe

if ispc

    python = fullfile( ...
        venv_folder, ...
        'Scripts', ...
        'python.exe');

else

    python = fullfile( ...
        venv_folder, ...
        'bin', ...
        'python3');

end

% Find ffmpeg
[ffmpeg_status, ffmpeg_path] = system('which ffmpeg');

if ffmpeg_status ~= 0
    error(['FFmpeg was not found on your PATH.' newline newline ...
           'Please install FFmpeg and make sure the ffmpeg executable ' ...
           'is available from the terminal.']);
end

ffmpeg_path = strtrim(ffmpeg_path);

fprintf('FFmpeg found at: %s\n', ffmpeg_path);

% Path to run_all.py

run_all = fullfile( ...
    code_folder, ...
    'run_all.py');


%% ========================================================================
% CHECK REQUIRED FILES
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('fUSI/FUS ANALYSIS PIPELINE SETUP\n');
fprintf('============================================================\n');


% -------------------------------------------------------------------------
% Check run_all.py
% -------------------------------------------------------------------------

if ~isfile(run_all)

    error( ...
        ['Could not find run_all.py.\n\n' ...
         'Expected location:\n' ...
         '  %s\n\n' ...
         'Make sure run_fUSI_analysis.m is located inside the\n' ...
         'fUSI_FUS_code folder.'], ...
        run_all);

end


% -------------------------------------------------------------------------
% Check virtual environment
% -------------------------------------------------------------------------

if ~isfolder(venv_folder)

    error( ...
        ['Could not find the Python virtual environment.\n\n' ...
         'Expected location:\n' ...
         '  %s\n\n' ...
         'Please create the virtual environment using the instructions\n' ...
         'in README.md.\n\n' ...
         'The expected project structure is:\n\n' ...
         '  fUSI_FUS_analysis/\n' ...
         '  |-- fUSI_FUS_code/\n' ...
         '  |-- fUSI_FUS_venv/'], ...
        venv_folder);

end


% -------------------------------------------------------------------------
% Check Python executable
% -------------------------------------------------------------------------

if ~isfile(python)

    error( ...
        ['Could not find the Python executable inside the project\n' ...
         'virtual environment.\n\n' ...
         'Expected:\n' ...
         '  %s\n\n' ...
         'The virtual environment may not have been created correctly.\n\n' ...
         'Please follow the environment setup instructions in README.md.'], ...
        python);

end


% -------------------------------------------------------------------------
% Check data
% -------------------------------------------------------------------------

if ~isfile(path_to_data) && ~isfolder(path_to_data)

    error( ...
        ['Could not find the specified data.\n\n' ...
         'Data path:\n' ...
         '  %s\n\n' ...
         'Check path_to_data at the top of this script.'], ...
        path_to_data);

end


%% ========================================================================
% DISPLAY PROJECT INFORMATION
% ========================================================================

fprintf('\n');
fprintf('Project folder:\n');
fprintf('  %s\n', project_folder);

fprintf('\n');
fprintf('Code folder:\n');
fprintf('  %s\n', code_folder);

fprintf('\n');
fprintf('Virtual environment:\n');
fprintf('  %s\n', venv_folder);

fprintf('\n');
fprintf('Python executable:\n');
fprintf('  %s\n', python);

fprintf('\n');
fprintf('run_all.py:\n');
fprintf('  %s\n', run_all);


%% ========================================================================
% CHECK PYTHON
% ========================================================================

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('Checking Python environment...\n');
fprintf('------------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Check Python executable and retrieve version/architecture
% -------------------------------------------------------------------------

python_info_command = sprintf( ...
    '"%s" -c "import sys, platform; print(sys.version.split()[0]); print(platform.machine()); print(sys.executable)"', ...
    python);


[status_python, python_info] = system(python_info_command);


if status_python ~= 0

    error( ...
        ['Python could not be started from the project virtual environment.\n\n' ...
         'Python executable:\n' ...
         '  %s\n\n' ...
         'This usually means the virtual environment is broken or was\n' ...
         'created with a Python installation that is no longer available.\n\n' ...
         'Try recreating the environment using the instructions in README.md.'], ...
        python);

end


% Split returned information into lines

python_lines = strsplit(strtrim(python_info), newline);


if numel(python_lines) >= 1
    python_version = strtrim(python_lines{1});
else
    python_version = 'Unknown';
end


if numel(python_lines) >= 2
    python_architecture = strtrim(python_lines{2});
else
    python_architecture = 'Unknown';
end


fprintf('\n');
fprintf('Python version:\n');
fprintf('  %s\n', python_version);

fprintf('\n');
fprintf('Python architecture:\n');
fprintf('  %s\n', python_architecture);

fprintf('\n');
fprintf('Python executable:\n');
fprintf('  %s\n', python);


%% ========================================================================
% CHECK NUMPY
% ========================================================================

fprintf('\n');
fprintf('Checking NumPy...\n');


numpy_command = sprintf( ...
    '"%s" -c "import numpy as np; print(np.__version__)"', ...
    python);


[status_numpy, numpy_version] = system(numpy_command);


if status_numpy ~= 0

    error( ...
        ['NumPy could not be imported.\n\n' ...
         'Python:\n' ...
         '  %s\n\n' ...
         'This usually means NumPy is not installed in the project\n' ...
         'virtual environment, or that NumPy was installed for a\n' ...
         'different CPU architecture.\n\n' ...
         'Try:\n\n' ...
         '  source fUSI_FUS_venv/bin/activate\n' ...
         '  pip install -r fUSI_FUS_code/requirements.txt\n\n' ...
         'If you are on Apple Silicon, make sure Python and NumPy\n' ...
         'are both using arm64.'], ...
        python);

end


fprintf('  NumPy %s', numpy_version);


%% ========================================================================
% CHECK SCIPY
% ========================================================================

fprintf('Checking SciPy...\n');


scipy_command = sprintf( ...
    '"%s" -c "import scipy; print(scipy.__version__)"', ...
    python);


[status_scipy, scipy_version] = system(scipy_command);


if status_scipy ~= 0

    error( ...
        ['SciPy could not be imported.\n\n' ...
         'Python:\n' ...
         '  %s\n\n' ...
         'SciPy is either not installed in the project virtual\n' ...
         'environment or there is an architecture/package problem.\n\n' ...
         'Try:\n\n' ...
         '  source fUSI_FUS_venv/bin/activate\n' ...
         '  pip install -r fUSI_FUS_code/requirements.txt'], ...
        python);

end


fprintf('  SciPy %s', scipy_version);

%% ========================================================================
% CHECK FFMpeg
% ========================================================================

fprintf('\nChecking FFmpeg...\n');

[ffmpeg_status, ffmpeg_path] = system('which ffmpeg');

if ffmpeg_status ~= 0
    error(['FFmpeg could not be found by MATLAB.\n\n' ...
           'Make sure /opt/homebrew/bin is on MATLAB''s PATH.']);
end

fprintf('  FFmpeg: %s', ffmpeg_path);

%% ========================================================================
% CHECK PYTHON PACKAGE IMPORTS USED BY PIPELINE
%
% This section can be expanded if the pipeline gains additional required
% packages.
% ========================================================================

fprintf('\n');
fprintf('Checking core Python packages...\n');


package_check_command = sprintf( ...
    ['"%s" -c "import numpy; import scipy; ' ...
     'import matplotlib; import h5py; ' ...
     'print(''OK'')"'], ...
    python);


[status_packages, package_output] = system(package_check_command);


if status_packages ~= 0

    fprintf('\n');
    fprintf('Python package check failed.\n');
    fprintf('\n');
    fprintf('%s\n', package_output);

    error( ...
        ['One or more required Python packages could not be imported.\n\n' ...
         'Make sure all packages listed in requirements.txt are installed:\n\n' ...
         '  pip install -r fUSI_FUS_code/requirements.txt\n\n' ...
         'See README.md for the complete environment setup instructions.']);

end


fprintf('  NumPy      OK\n');
fprintf('  SciPy      OK\n');
fprintf('  Matplotlib OK\n');
fprintf('  h5py       OK\n');


%% ========================================================================
% FINISHED ENVIRONMENT CHECK
% ========================================================================

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('Python environment check PASSED\n');
fprintf('------------------------------------------------------------\n');


%% ========================================================================
% DISPLAY ANALYSIS SETTINGS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('ANALYSIS SETTINGS\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('Data:\n');
fprintf('  %s\n', path_to_data);

fprintf('\n');
fprintf('Acquisition:\n');
fprintf('  %s\n', ACQ_to_run);

fprintf('\n');
fprintf('RHS identifier:\n');
fprintf('  %s\n', RHS_identifier);

fprintf('\n');
fprintf('Reconstruction:\n');
fprintf('  %s\n', reconstruction_type);

fprintf('\n');
fprintf('Extra identifier:\n');
fprintf('  %s\n', extra_identifier);

fprintf('\n');
fprintf('============================================================\n');


%% ========================================================================
% BUILD PYTHON COMMAND
% ========================================================================

command = sprintf( ...
    '"%s" "%s" --date "%s" --data-path "%s" --acq "%s" --rhs "%s" --reconstruction "%s" --extra-identifier "%s"', ...
    python, ...
    run_all, ...
    date, ...
    path_to_data, ...
    ACQ_to_run, ...
    RHS_identifier, ...
    reconstruction_type, ...
    extra_identifier);


%% ========================================================================
% RUN PYTHON
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('STARTING PYTHON PIPELINE\n');
fprintf('============================================================\n');

fprintf('\n');

status = system(command);


%% ========================================================================
% CHECK RESULT
% ========================================================================

fprintf('\n');

if status == 0

    fprintf('============================================================\n');
    fprintf('fUSI/FUS ANALYSIS COMPLETED SUCCESSFULLY\n');
    fprintf('============================================================\n');

else

    fprintf('============================================================\n');
    fprintf('fUSI/FUS ANALYSIS FAILED\n');
    fprintf('============================================================\n');

    fprintf('\n');
    fprintf('Python returned error code: %d\n', status);

    fprintf('\n');
    fprintf('The Python environment itself passed the initial checks.\n');
    fprintf('The error above therefore occurred while running the\n');
    fprintf('analysis pipeline.\n');

    error( ...
        'Python pipeline returned error code %d.', ...
        status);

end

