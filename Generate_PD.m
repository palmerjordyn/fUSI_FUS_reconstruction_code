
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
%   6. Display all 8 PD images with adjustable dB slider
%   7. Display an individual PD image with its own dB slider
%   8. Select signal/background ROIs
%   9. Calculate SNR and CNR
%
% The sliders change display contrast only.
% SNR/CNR calculations use linear power values.

clear;
clc;
close all;


%% ==================================================
% BF FILE PATH
% ==================================================

bfPath = "/Users/jpalmer/Downloads/Hobbes Oct 9th 2026/Acq_2026_10_09_12_12_53/P7-4__BF_214.bin";


%% ==================================================
% Expected dimensions
% ==================================================

NZ = 339;
NX = 38;
NT = 200;


%% ==================================================
% Check BF file
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

% Python dtype "<c16":
% 64-bit little-endian real + 64-bit imaginary

raw = fread(fid, Inf, "double=>double");

fclose(fid);


%% ==================================================
% Convert binary data to complex IQ
% ==================================================

if mod(numel(raw), 2) ~= 0

    error("Binary file contains an odd number of double values.");

end

data = complex( ...
    raw(1:2:end), ...
    raw(2:2:end));

expected = NZ * NX * NT;

if numel(data) ~= expected

    error( ...
        "Expected %d complex samples, but found %d.", ...
        expected, numel(data));

end


%% ==================================================
% Reshape into Z x X x time
% ==================================================

cube = reshape(data, [NZ, NX, NT]);

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
% Apply Butterworth filter along time dimension
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
    power_doppler_safe ./ max(power_doppler_safe(:)));


%% ==================================================
% Display all PD images with adjustable dB slider
% ==================================================

figAll = figure( ...
    "Name", "All Power Doppler Images", ...
    "Position", [50 200 1600 500]);

% Leave room at the bottom for the slider.
tl = tiledlayout(figAll, 1, N_IMAGES, ...
    "TileSpacing", "compact", ...
    "Padding", "compact");

tl.OuterPosition = [0.01 0.13 0.94 0.82];

% Initial display limits.
vmin = -60;
vmax = 0;

axAll = gobjects(1, N_IMAGES);

for i = 1:N_IMAGES

    axAll(i) = nexttile(tl);

    imagesc( ...
        axAll(i), ...
        squeeze(power_db(i, :, :)), ...
        [vmin vmax]);

    colormap(axAll(i), gray);
    axis(axAll(i), "image");

    title(axAll(i), sprintf( ...
        "PD %d\n(IQ frames %d-%d)", ...
        i, ...
        (i - 1) * WINDOW + 1, ...
        i * WINDOW));

    xlabel(axAll(i), "X");
    ylabel(axAll(i), "Z");

end

cb = colorbar(axAll(end));
cb.Label.String = "Power (dB)";
cb.Layout.Tile = "east";

sgtitle(tl, "Butterworth High-Pass Power Doppler");


%% ==================================================
% Slider for all 8 images
% ==================================================

dbLabelAll = uicontrol(figAll, ...
    "Style", "text", ...
    "Units", "normalized", ...
    "Position", [0.10 0.045 0.20 0.045], ...
    "String", sprintf("Lower dB limit: %d dB", vmin), ...
    "FontSize", 11, ...
    "BackgroundColor", get(figAll, "Color"));

sliderAll = uicontrol(figAll, ...
    "Style", "slider", ...
    "Units", "normalized", ...
    "Position", [0.32 0.05 0.48 0.04], ...
    "Min", -80, ...
    "Max", 0, ...
    "Value", vmin, ...
    "SliderStep", [1/80 5/80], ...
    "Callback", @(src, evt) ...
        updateDbLimit(src, axAll, dbLabelAll));

uicontrol(figAll, ...
    "Style", "text", ...
    "Units", "normalized", ...
    "Position", [0.25 0.045 0.07 0.045], ...
    "String", "-80 dB", ...
    "BackgroundColor", get(figAll, "Color"));

uicontrol(figAll, ...
    "Style", "text", ...
    "Units", "normalized", ...
    "Position", [0.80 0.045 0.06 0.045], ...
    "String", "0 dB", ...
    "BackgroundColor", get(figAll, "Color"));


%% ==================================================
% ROI ANALYSIS
% ==================================================

fprintf("\n========================================\n");
fprintf("ROI ANALYSIS\n");
fprintf("========================================\n");


%% ==================================================
% Select ROI radius
% ==================================================

radius = input("\nEnter ROI radius in pixels: ");

if isempty(radius) || ~isscalar(radius) || ...
        ~isfinite(radius) || radius <= 0

    error("ROI radius must be a positive number.");

end


%% ==================================================
% Individual PD image with adjustable dB slider
% ==================================================

% Choose which PD image to display for ROI selection.
roiImageIndex = 3;

if roiImageIndex < 1 || roiImageIndex > N_IMAGES

    error("roiImageIndex must be between 1 and %d.", N_IMAGES);

end

figROI = figure( ...
    "Name", "Individual PD Image - ROI Analysis", ...
    "Position", [250 100 750 750]);

axROI = axes(figROI, ...
    "Units", "normalized", ...
    "Position", [0.12 0.20 0.78 0.72]);

imagesc( ...
    axROI, ...
    squeeze(power_db(roiImageIndex, :, :)), ...
    [vmin vmax]);

colormap(axROI, gray);
axis(axROI, "image");

xlabel(axROI, "X");
ylabel(axROI, "Z");

title(axROI, sprintf( ...
    "PD %d — Click center of SIGNAL ROI", ...
    roiImageIndex));


%% ==================================================
% Slider for individual PD image
% ==================================================

dbLabelROI = uicontrol(figROI, ...
    "Style", "text", ...
    "Units", "normalized", ...
    "Position", [0.10 0.105 0.40 0.045], ...
    "String", sprintf("Lower dB limit: %d dB", vmin), ...
    "FontSize", 11, ...
    "BackgroundColor", get(figROI, "Color"));

sliderROI = uicontrol(figROI, ...
    "Style", "slider", ...
    "Units", "normalized", ...
    "Position", [0.12 0.065 0.76 0.04], ...
    "Min", -80, ...
    "Max", 0, ...
    "Value", vmin, ...
    "SliderStep", [1/80 5/80], ...
    "Callback", @(src, evt) ...
        updateDbLimit(src, axROI, dbLabelROI));

uicontrol(figROI, ...
    "Style", "text", ...
    "Units", "normalized", ...
    "Position", [0.02 0.06 0.09 0.045], ...
    "String", "-80 dB", ...
    "BackgroundColor", get(figROI, "Color"));

uicontrol(figROI, ...
    "Style", "text", ...
    "Units", "normalized", ...
    "Position", [0.88 0.06 0.09 0.045], ...
    "String", "0 dB", ...
    "BackgroundColor", get(figROI, "Color"));

drawnow;


%% ==================================================
% Select SIGNAL ROI
% ==================================================

fprintf("\nAdjust the individual image dB slider if desired.\n");
fprintf("Then click the center of the SIGNAL ROI.\n");

[signal_x, signal_z] = selectImagePoint(axROI, ...
    "Click the center of the SIGNAL ROI in the image.");

if isempty(signal_x)

    error("No signal ROI was selected.");

end

fprintf( ...
    "Signal ROI center: X = %.2f, Z = %.2f\n", ...
    signal_x, signal_z);


%% ==================================================
% Draw SIGNAL ROI
% ==================================================

hold(axROI, "on");

theta = linspace(0, 2*pi, 200);

plot( ...
    axROI, ...
    signal_x + radius*cos(theta), ...
    signal_z + radius*sin(theta), ...
    "r-", ...
    "LineWidth", 2);


%% ==================================================
% Select BACKGROUND ROI
% ==================================================

title(axROI, sprintf( ...
    "PD %d — Click center of BACKGROUND ROI", ...
    roiImageIndex));

drawnow;

fprintf("\nClick the center of the BACKGROUND ROI.\n");

[background_x, background_z] = selectImagePoint(axROI, ...
    "Click the center of the BACKGROUND ROI in the image.");

if isempty(background_x)

    error("No background ROI was selected.");

end

fprintf( ...
    "Background ROI center: X = %.2f, Z = %.2f\n", ...
    background_x, background_z);


%% ==================================================
% Draw BACKGROUND ROI
% ==================================================

plot( ...
    axROI, ...
    background_x + radius*cos(theta), ...
    background_z + radius*sin(theta), ...
    "y-", ...
    "LineWidth", 2);

title(axROI, sprintf( ...
    "PD %d — Signal (red) and Background (yellow)", ...
    roiImageIndex));

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
% Use linear power values, not dB values.
% Use the same PD image selected for ROI placement.
% ==================================================

PD1 = squeeze( ...
    power_doppler(roiImageIndex, :, :));

signal_values = PD1(signal_mask);

background_values = PD1(background_mask);


%% ==================================================
% Compute ROI statistics
% ==================================================

if isempty(signal_values) || isempty(background_values)

    error("One of the ROIs contains no pixels.");

end

usignal = mean(signal_values);
stdsignal = std(double(signal_values));

ubackground = mean(background_values);
stdb = std(double(background_values));


%% ==================================================
% Compute SNR
% ==================================================

if stdb > 0 && usignal > 0

    SNRdB = 20 * log10(usignal / stdb);

else

    SNRdB = NaN;

end


%% ==================================================
% Compute CNR
% ==================================================

cnrDenominator = sqrt( ...
    (stdsignal^2 + stdb^2) / 2);

if cnrDenominator > 0

    CNRdB = 20 * log10( ...
        abs(usignal - ubackground) / cnrDenominator);

else

    CNRdB = NaN;

end


%% ==================================================
% Print ROI results
% ==================================================

fprintf("\n========================================\n");
fprintf("ROI RESULTS\n");
fprintf("========================================\n");

fprintf("\nPD image used: %d\n", roiImageIndex);

fprintf("\nSignal ROI\n");
fprintf("----------------------------------------\n");

fprintf("Center X:      %.2f\n", signal_x);
fprintf("Center Z:      %.2f\n", signal_z);
fprintf("Radius:        %.2f pixels\n", radius);
fprintf("Pixels:        %d\n", numel(signal_values));
fprintf("Mean:          %.6g\n", usignal);
fprintf("Std:           %.6g\n", stdsignal);

fprintf("\nBackground ROI\n");
fprintf("----------------------------------------\n");

fprintf("Center X:      %.2f\n", background_x);
fprintf("Center Z:      %.2f\n", background_z);
fprintf("Radius:        %.2f pixels\n", radius);
fprintf("Pixels:        %d\n", numel(background_values));
fprintf("Mean:          %.6g\n", ubackground);
fprintf("Std:           %.6g\n", stdb);

fprintf("\nImage quality\n");
fprintf("----------------------------------------\n");

fprintf("SNR:           %.3f dB\n", SNRdB);
fprintf("CNR:           %.3f dB\n", CNRdB);

fprintf("========================================\n");


%% ==================================================
% Keep figures open
% ==================================================

drawnow;


%% ==================================================
% Local function: update dB display limit
% ==================================================

function updateDbLimit(src, ax, labelHandle)

    lowerDb = round(get(src, "Value"));

    % Update the display contrast without changing image data.
    for k = 1:numel(ax)

        if isgraphics(ax(k), "axes")

            set(ax(k), "CLim", [lowerDb 0]);

        end

    end

    % Update the text above the slider.
    set(labelHandle, "String", ...
        sprintf("Lower dB limit: %d dB", lowerDb));

    drawnow limitrate;

end




%% ==================================================
% Select a point only when the image axes are clicked
% ==================================================

function [x, z] = selectImagePoint(ax, prompt)

    fig = ancestor(ax, "figure");

    fprintf("%s\n", prompt);
    title(ax, prompt);

    x = [];
    z = [];

    % Wait for clicks until one lands inside the image axes.
    while isempty(x)

        % Wait for the next mouse click or key press.
        waitforbuttonpress;

        % Ignore keyboard presses.
        if ~strcmp(get(fig, "SelectionType"), "normal")
            continue;
        end

        % Identify the object under the mouse pointer.
        clickedObject = hittest(fig);

        % Ignore clicks on sliders, labels, buttons, etc.
        clickedAxes = ancestor(clickedObject, "axes");

        if isempty(clickedAxes) || ~isequal(clickedAxes, ax)
            continue;
        end

        % Get image coordinates.
        cp = get(ax, "CurrentPoint");

        xClick = cp(1, 1);
        zClick = cp(1, 2);

        % Only accept clicks inside the image boundaries.
        xLimits = xlim(ax);
        zLimits = ylim(ax);

        if xClick < min(xLimits) || xClick > max(xLimits) || ...
           zClick < min(zLimits) || zClick > max(zLimits)
            continue;
        end

        x = xClick;
        z = zClick;

    end

    fprintf("Selected X = %.2f, Z = %.2f\n", x, z);

end
