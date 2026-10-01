%% Butterworth High-Pass Power Doppler + ROI Analysis
%
% Reads a BF file from a normal folder path.
%
% Processing:
%   1. Read complex BF data
%   2. Reshape into Z x X x time
%   3. Subtract first frame
%   4. Apply Butterworth high-pass filter
%   5. Generate 25-frame power-Doppler images
%   6. Display all 8 PD images
%   7. Select signal/background ROIs
%   8. Calculate SNR and CNR

clear;
clc;
close all;


%% ==================================================
% BF FILE PATH
%
% Set this to the complete path to the BF file.
%
% Example:
%
% bfPath = "/Users/jpalmer/Downloads/HobbesfusiFUSAug312026/Acq_2026_08_31_14_46_09/P7-4__BF_005.bin";
% ==================================================

bfPath = "/Users/jpalmer/Downloads/HobbesfusiFUSAug312026/Acq_2026_08_31_14_46_09/P7-4__BF_005.bin";


%% ==================================================
% Expected dimensions
% ==================================================

NZ = 339;
NX = 38;
NT = 200;


%% ==================================================
% Construct full BF file path
% ==================================================


if ~isfile(bfPath)

    error( ...
        "Could not find BF file:\n%s", ...
        bfPath);

end

fprintf("Reading BF file:\n%s\n", bfPath);


%% ==================================================
% Read binary BF data
% ==================================================

fid = fopen(bfPath, "r", "ieee-le");

if fid == -1
    error("Could not open BF file:\n%s", bfPath);
end

% Python dtype "<c16" consists of:
%   64-bit little-endian real
%   64-bit little-endian imaginary
%
% Read all values as doubles and reconstruct the complex samples.

raw = fread(fid, Inf, "double=>double");

fclose(fid);


%% ==================================================
% Convert binary data to complex IQ
% ==================================================

if mod(numel(raw), 2) ~= 0

    error( ...
        "Binary file contains an odd number of double values.");

end

data = complex( ...
    raw(1:2:end), ...
    raw(2:2:end));


expected = NZ * NX * NT;

if numel(data) ~= expected

    error( ...
        "Expected %d complex samples, but found %d.", ...
        expected, ...
        numel(data));

end


%% ==================================================
% Reshape into Z x X x time
%
% NumPy order='F' corresponds to MATLAB's native
% column-major reshape ordering.
% ==================================================

cube = reshape( ...
    data, ...
    [NZ, NX, NT]);

fprintf( ...
    "Loaded cube: %d x %d x %d\n", ...
    size(cube));


%% ==================================================
% Butterworth high-pass filter parameters
% ==================================================

fs = 250;       % Sampling frequency [Hz]
cutoff = 30;    % High-pass cutoff [Hz]
filterOrder = 4;


%% ==================================================
% Create Butterworth filter
% ==================================================

[B, A] = butter( ...
    filterOrder, ...
    cutoff / (fs / 2), ...
    "high");


fprintf("\nButterworth filter coefficients:\n");

disp("B =");
disp(B);

disp("A =");
disp(A);


%% ==================================================
% Subtract first frame from all frames
% ==================================================

firstFrame = cube(:, :, 1);

cube = cube - firstFrame;


%% ==================================================
% Apply Butterworth filter to every pixel
%
% Dimension 3 is the temporal/IQ dimension.
%
% This corresponds to:
%
% Python:
% filtered = lfilter(B, A, cube, axis=2)
%
% MATLAB:
% filtered = filter(B, A, cube, [], 3)
% ==================================================

filtered = filter( ...
    B, ...
    A, ...
    cube, ...
    [], ...
    3);


fprintf( ...
    "Filtered cube: %d x %d x %d\n", ...
    size(filtered));


%% ==================================================
% Generate power-Doppler images
% ==================================================

WINDOW = 25;

N_IMAGES = floor(NT / WINDOW);

power_doppler = zeros( ...
    N_IMAGES, ...
    NZ, ...
    NX, ...
    "single");


for i = 1:N_IMAGES

    startIdx = (i - 1) * WINDOW + 1;
    endIdx = i * WINDOW;

    q = filtered(:, :, startIdx:endIdx);

    % Sum squared magnitude over the 25 IQ frames.
    power_doppler(i, :, :) = single( ...
        sum(abs(q).^2, 3));

end


fprintf( ...
    "Power Doppler shape: %d x %d x %d\n", ...
    size(power_doppler));


%% ==================================================
% Convert power Doppler to dB
% ==================================================

power_doppler_safe = max( ...
    power_doppler, ...
    single(1e-12));


power_db = 10 * log10( ...
    power_doppler_safe ./ ...
    max(power_doppler_safe(:)));


%% ==================================================
% Display all 8 PD images
% ==================================================

figure( ...
    "Name", ...
    "Butterworth High-Pass Power Doppler");


tiledlayout( ...
    1, ...
    8, ...
    "TileSpacing", "compact", ...
    "Padding", "compact");


% Same dB scale for all 8 images
vmin = -40;
vmax = 0;


for i = 1:N_IMAGES

    nexttile;

    imagesc( ...
        squeeze(power_db(i, :, :)), ...
        [vmin vmax]);

    colormap(gray);

    axis image;

    title(sprintf( ...
        "PD %d\n(IQ frames %d-%d)", ...
        i, ...
        (i - 1) * WINDOW + 1, ...
        i * WINDOW));

    xlabel("X");
    ylabel("Z");

end


cb = colorbar;
cb.Layout.Tile = "east";
cb.Label.String = "Power (dB)";


sgtitle( ...
    "BF_005 — Butterworth High-Pass Power Doppler");


%% ==================================================
% ROI ANALYSIS
% ==================================================

fprintf("\n========================================\n");
fprintf("ROI ANALYSIS\n");
fprintf("========================================\n");


%% ==================================================
% Ask for ROI radius
% ==================================================

radius = input( ...
    "\nEnter ROI radius in pixels: ");


if isempty(radius) || radius <= 0

    error( ...
        "ROI radius must be greater than zero.");

end


%% ==================================================
% Open second window showing PD image 1
% ==================================================

figure( ...
    "Name", ...
    "ROI Analysis");


imagesc( ...
    squeeze(power_db(1, :, :)), ...
    [vmin vmax]);

colormap(gray);

axis image;

xlabel("X");
ylabel("Z");

title( ...
    "PD 1 — Click center of SIGNAL ROI");


%% ==================================================
% Select SIGNAL ROI
% ==================================================

fprintf( ...
    "\nClick the center of the SIGNAL ROI.\n");

[signal_x, signal_z] = ginput(1);


if isempty(signal_x)

    error( ...
        "No signal ROI was selected.");

end


fprintf( ...
    "Signal ROI center: X = %.2f, Z = %.2f\n", ...
    signal_x, ...
    signal_z);


%% ==================================================
% Draw SIGNAL ROI
% ==================================================

hold on;

theta = linspace(0, 2*pi, 200);

plot( ...
    signal_x + radius*cos(theta), ...
    signal_z + radius*sin(theta), ...
    "LineWidth", 2);


%% ==================================================
% Select BACKGROUND ROI
% ==================================================

title( ...
    "PD 1 — Click center of BACKGROUND ROI");

drawnow;

fprintf( ...
    "\nClick the center of the BACKGROUND ROI.\n");


[background_x, background_z] = ginput(1);


if isempty(background_x)

    error( ...
        "No background ROI was selected.");

end


fprintf( ...
    "Background ROI center: X = %.2f, Z = %.2f\n", ...
    background_x, ...
    background_z);


%% ==================================================
% Draw BACKGROUND ROI
% ==================================================

plot( ...
    background_x + radius*cos(theta), ...
    background_z + radius*sin(theta), ...
    "LineWidth", 2);


title( ...
    "PD 1 — Signal and Background ROIs");

drawnow;


%% ==================================================
% Create circular ROI masks
% ==================================================

[x_grid, z_grid] = meshgrid( ...
    1:NX, ...
    1:NZ);


signal_mask = ...
    (x_grid - signal_x).^2 + ...
    (z_grid - signal_z).^2 <= radius^2;


background_mask = ...
    (x_grid - background_x).^2 + ...
    (z_grid - background_z).^2 <= radius^2;


%% ==================================================
% Extract signal and background values
%
% Use LINEAR power values for the statistics,
% not the dB-converted image.
% ==================================================

PD1 = squeeze( ...
    power_doppler(1, :, :));


signal_values = ...
    PD1(signal_mask);


background_values = ...
    PD1(background_mask);


%% ==================================================
% Compute ROI statistics
% ==================================================

usignal = mean(signal_values);
stdsignal = std(signal_values);

ubackground = mean(background_values);
stdb = std(background_values);


%% ==================================================
% Compute SNR
% ==================================================

SNRdB = 20 * log10( ...
    usignal / stdb);


%% ==================================================
% Compute CNR
% ==================================================

CNRdB = 20 * log10( ...
    abs(usignal - ubackground) / ...
    sqrt( ...
        (stdsignal^2 + stdb^2) / 2));


%% ==================================================
% Print ROI results
% ==================================================

fprintf("\n========================================\n");
fprintf("ROI RESULTS\n");
fprintf("========================================\n");


fprintf("\nSignal ROI\n");
fprintf("----------------------------------------\n");

fprintf( ...
    "Center X:      %.2f\n", ...
    signal_x);

fprintf( ...
    "Center Z:      %.2f\n", ...
    signal_z);

fprintf( ...
    "Radius:        %.2f pixels\n", ...
    radius);

fprintf( ...
    "Pixels:        %d\n", ...
    numel(signal_values));

fprintf( ...
    "Mean:          %.6g\n", ...
    usignal);

fprintf( ...
    "Std:           %.6g\n", ...
    stdsignal);


fprintf("\nBackground ROI\n");
fprintf("----------------------------------------\n");

fprintf( ...
    "Center X:      %.2f\n", ...
    background_x);

fprintf( ...
    "Center Z:      %.2f\n", ...
    background_z);

fprintf( ...
    "Radius:        %.2f pixels\n", ...
    radius);

fprintf( ...
    "Pixels:        %d\n", ...
    numel(background_values));

fprintf( ...
    "Mean:          %.6g\n", ...
    ubackground);

fprintf( ...
    "Std:           %.6g\n", ...
    stdb);


fprintf("\nImage quality\n");
fprintf("----------------------------------------\n");

fprintf( ...
    "SNR:           %.3f dB\n", ...
    SNRdB);

fprintf( ...
    "CNR:           %.3f dB\n", ...
    CNRdB);

fprintf("========================================\n");


%% ==================================================
% Keep figures open
% ==================================================

drawnow;
