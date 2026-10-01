
%% run_fUSI_analysis.m
%
% MATLAB entry point for the fUSI/FUS analysis pipeline.
%
% This script calls run_all.py.
%
% The Python pipeline can still be run directly from VS Code with:
%
%     python run_all.py
%
% In that case, the default values in config.py are used.
%
% When this MATLAB script is run, the values below override the
% corresponding values in config.py.

clear;
clc;


%% ========================================================================
% USER INPUT
% ========================================================================

% Path to the ZIP file OR unzipped data directory
path_to_data = ...
    '/Users/jpalmer/Downloads/HobbesFusiFusSeptember16th2026.zip';


% Acquisition to analyze
ACQ_to_run = ...
    'Acq_2026_09_16_10_57_16';


% RHS identifier
RHS_identifier = ...
    'TB_take2';


% Reconstruction type
%
% Examples:
%   'BUTTER'
%   'SVD'
%
reconstruction_type = ...
    'BUTTER';


% Optional extra identifier for the output folder
extra_identifier = ...
    'theta';


%% ========================================================================
% PYTHON SETUP
%
% These should normally only need to be configured once on each computer.
% ========================================================================

% Python executable
%
% IMPORTANT:
% Replace this with the Python executable from the environment that
% contains all of your required packages.

python = '/Library/Frameworks/Python.framework/Versions/3.13/bin/python3';


% Path to run_all.py

run_all = fullfile(fileparts(mfilename('fullpath')), 'run_all.py');

%% ========================================================================
% CHECK THAT EVERYTHING EXISTS
% ========================================================================

if ~isfile(run_all)
    error( ...
        'Could not find run_all.py:\n%s', ...
        run_all);
end


if ~isfile(python)
    error( ...
        'Could not find Python executable:\n%s', ...
        python);
end


if ~isfile(path_to_data) && ~isfolder(path_to_data)
    error( ...
        'Could not find data:\n%s', ...
        path_to_data);
end


%% ========================================================================
% DISPLAY SETTINGS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('fUSI/FUS ANALYSIS PIPELINE\n');
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
    '"%s" "%s" --data-path "%s" --acq "%s" --rhs "%s" --reconstruction "%s" --extra-identifier "%s"', ...
    python, ...
    run_all, ...
    path_to_data, ...
    ACQ_to_run, ...
    RHS_identifier, ...
    reconstruction_type, ...
    extra_identifier);


%% ========================================================================
% RUN PYTHON
% ========================================================================

fprintf('\n');
fprintf('Starting Python pipeline...\n');
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

    error( ...
        'Python pipeline returned error code %d.', ...
        status);

end
