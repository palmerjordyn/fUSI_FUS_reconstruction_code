%% B-mode / MRI Overlay -- Interactive Multi-Angle Viewer + fUSI Analysis
%
% Overlays a phased-array B-mode image onto an oblique MRI slice
% extracted from a 3D MRI volume (sys structure, same as overlayUsFusi).
%
% Bottom row shows fUSI Doppler overlaid on the same MRI slice,
% positioned using the same B-mode adjustment controls.
%
% If a session folder has no Dop.bin, it is computed from the beamformed
% (_BF_) files before analysis, using the same power-Doppler steps as
% fUS_correlationZnorm.m.
%
% The Controls panel on the right runs analysis, switches session and
% animal, and opens the summary figures. Up/Down arrows change the beam
% angle and Left/Right arrows step through sagittal planes in the Beam Path
% panel. W/A/S/D move the ultrasound image (hold Shift for smaller steps).
% R and E rotate the image.
%
% Set fUSI = false in INPUTS to use this only for aligning the B-mode with
% the MRI. No fUSI data is read or computed, and the fUSI panels, analysis
% buttons and analysis hotkeys are turned off.
%
% COORDINATE SYSTEM (matches overlayUsFusi):
%   X-axis  = lateral position (mm), centered on transducer
%   Z-axis  = depth from transducer face along beam (mm), 0 = face

clear; close all force; clc;   % "force" skips old figures' CloseRequestFcn (its script functions are gone after clear)

%% ==================== INPUTS ============================
subject = {'Hobbes'};

% fUSI = true : full pipeline (load fUSI, activation maps, analysis, fUSI overlays)
% fUSI = false: B-mode / MRI alignment only. No fUSI data is read or computed and
%               sessionDefs(i).fusiDir is ignored (it can be left empty).
fUSI = false;

% sessionDefs(1).bMDir   = '/Users/jpalmer/Downloads/12-10-10-16-26_Boltzmann';
% sessionDefs(1).bMFile  = 'BMIQ_005';
% sessionDefs(1).fusiDir = '/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected/Acq_2026_10_05_10_58_30/Acq_2026_10_05_10_58_45';
sessionDefs = struct('bMDir', {}, 'bMFile', {}, 'fusiDir', {});
sessionDefs(1).bMDir   = '/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected/10-55-10-05-26';
sessionDefs(1).bMFile  = 'BMIQ_001';
sessionDefs(1).fusiDir = '/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected/Acq_2026_10_05_10_58_30/Acq_2026_10_05_10_58_45';

sessionDefs(2).bMDir = '/System/Volumes/Data/v/raid1b/backup/jpalmer/FUS/HobbesFusiFusSeptember16th2026/10-19-09-16-26';
sessionDefs(2).bMFile  = 'BMIQ_001';



%% ===================== SHARED SETTINGS =====================
scriptDir = fileparts(mfilename('fullpath'));
% --- Portable output location (no user-specific paths) ---
% Expected folder layout, relative to wherever this script lives:
%   <analysisRoot>/fUSI_FUS_reconstruction_code/LightsReconCode   <- this script
%   <analysisRoot>/fUSI_FUS_output/LightsRecon                    <- saveDir

reconCodeDir = fileparts(scriptDir);       % .../fUSI_FUS_reconstruction_code
analysisRoot = fileparts(reconCodeDir);    % folder that holds both code and output
[~, scriptFolderName] = fileparts(scriptDir);
[~, reconFolderName]  = fileparts(reconCodeDir);
if ~strcmp(scriptFolderName, 'LightsReconCode') || ~strcmp(reconFolderName, 'fUSI_FUS_reconstruction_code')
    warning(['FusiOverlay.m is not inside fUSI_FUS_reconstruction_code/LightsReconCode.\n' ...
             'Output folder will be created relative to: %s'], scriptDir);
end

saveDir = fullfile(analysisRoot, 'fUSI_FUS_output', 'LightsRecon');
cacheDir = fullfile(saveDir, 'cache');     % precomputed MRI slices are stored here
if ~isfolder(cacheDir), mkdir(cacheDir); end
forceRecomputeMRI = false;                 % set true to ignore the cache and recompute


% --- fUSI display (shared) ---
fusiAlpha    = 0.5;
fusiColormap = 'hot';

% --- Overlay display options ---
overlayAlpha     = 0.6;
overlayMode      = 'alpha';
fusionGamma      = 0.6;
overlayColormap  = 'hot';
checkerBlockSize = 16;

% --- Grid padding for shift accommodation ---
gridPad_mm = 20;           % Extra mm on each side so MRI fills view during shifts

%% ===================== LOAD BOTH ANIMALS =====================
%subjects = {'Calvin', 'Hobbes', 'Boltzmann'};
animalStore = cell(1, 1);

for iAnimal = 1:1


%% ===================== SUBJECT PROFILES =====================
subject = subject{iAnimal};
%sessionDefs = struct('bMDir', {}, 'bMFile', {}, 'fusiDir', {});  % reset per animal

switch subject
    case 'Calvin'

        % --- Per-session definitions (add rows for more sessions) ---
        %   Each session needs its own B-mode capture and fUSI directory.


        mrilogFile = fullfile(scriptDir, 'CalvinMRI.mat');

        % --- B-mode scan geometry (phased array) ---
        faceWidth_mm = 12;
        fovWidth_mm  = 125;
        fovDepth_mm  = 80;
        dynRange_dB  = 60;

        % --- Beam / transducer positioning ---
        yCenter    = 40.05e-3;    % +10mm Y from original 30.05e-3
        depthParam = 0.012;
        zDistAdj   = 0.008;       % +10mm Z offset applied to sys.zDist
        angles     = 40:0.5:57;

        % --- B-mode orientation corrections ---
        doRot180_US = false;
        flipLR_US   = true;
        flipUD_US   = false;

        % --- Initial B-mode overlay offset (WASD/R/E adjustments, default for all sessions) ---
        bm_dx_init   = 9.5;
        bm_dz_init   = 9.5;
        bm_rot_init  = 7.0;
        defaultAngle = 52;

        % --- fUSI display options ---
        fusiFrameIdx   = 100;
        fusiXOffset_mm = 0;

        % --- fUSI analysis parameters ---
        stim.zscore     = 1.65;
        %stim.zscore     = 1.28;

        stim.isSmooth   = false;
        stim.smoothN    = 2;
        stim.actPeriod  = 10;
        stim.restPeriod = 10;
        stim.numAct     = 15;
        stim.stimFirst  = false;

        % --- Transducer face marker ---
        txMarkerThick_mm = 1.5;
        txFaceOffset_mm  = 7;

        % --- Epoch plot labels ---
        epochTitlePrefix = 'fUSI Neurovascular Response';
        subjectLabel     = 'Subject C';

        % --- Intan EEG / Accelerometer (optional) ---
        intanPath     = '';   % No Intan data for Calvin yet
        intanBaseName = {};
        intanAccelCh  = 1;
        intanFs       = 20000;

    case 'Hobbes'
        % --- Paths ---
        % Old single-session paths (for reference):
        % bMDir      = 'D:\fusFusiData\April 14th Hobbes Pilot\15-54-04-14-26';
        % fusiDir    = 'D:\fusFusiData\April 14th Hobbes Pilot\Acq_2026_04_14_15_56_28\Acq_2026_04_14_15_57_16\Acq_2026_04_14_16_01_27';

        % --- Per-session definitions (add rows for more sessions) ---
        sessionDefs(1).bMDir   = '/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected/10-55-10-05-26';
        sessionDefs(1).bMFile  = 'BMIQ_001';
        sessionDefs(1).fusiDir = '/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected/Acq_2026_10_05_10_58_30/Acq_2026_10_05_10_58_45';
        % sessionDefs(2).bMDir   = 'D:\fUS Data\Hobbes July 8th 2026\11-30-07-08-26';
        % sessionDefs(2).bMFile  = 'BMIQ_001';
        % sessionDefs(2).fusiDir = 'D:\fUS Data\Hobbes July 8th 2026\Acq_2026_07_08_11_18_25\Acq_2026_07_08_11_18_59\Acq_2026_07_08_11_23_17';
        % sessionDefs(3).bMDir   = 'D:\fUS Data\Hobbes July 6th 2026\...';
        % sessionDefs(3).bMFile  = 'BMIQ_001';
        % sessionDefs(3).fusiDir = 'D:\fUS Data\Hobbes July 6th 2026\Acq_2026_07_06_11_01_47\Acq_2026_07_06_11_02_56\Acq_2026_07_06_11_07_08';

        mrilogFile = fullfile(scriptDir, 'HobbesMRI.mat');

        % --- B-mode scan geometry (phased array) ---
        faceWidth_mm = 12;
        fovWidth_mm  = 120;
        fovDepth_mm  = 80;
        dynRange_dB  = 60;

        % --- Beam / transducer positioning ---
        yCenter    = 40.05e-3;    % +10mm Y from original 30.05e-3
        depthParam = 0.012;
        zDistAdj   = 0.007;       % +10mm Z offset
        angles     = 34:0.5:56;

        % --- B-mode orientation corrections ---
        doRot180_US = false;
        flipLR_US   = true;
        flipUD_US   = false;

        % --- Initial B-mode overlay offset (WASD/R/E adjustments, default for all sessions) ---
        bm_dx_init   = 9.5;
        bm_dz_init   = 9.5;
        bm_rot_init  = 7.0;
        defaultAngle = 52;

        % --- fUSI display options ---
        fusiFrameIdx   = 0;
        fusiXOffset_mm = 0;

        % --- fUSI analysis parameters ---
        stim.zscore     = 1.65;
        stim.isSmooth   = false;
        stim.smoothN    = 2;
        stim.actPeriod  = 10;
        stim.restPeriod = 10;
        stim.numAct     = 15;
        %stim.numAct     = 11;

        stim.stimFirst  = false;

        % --- Transducer face marker ---
        txMarkerThick_mm = 1.5;
        txFaceOffset_mm  = 7;

        % --- Epoch plot labels ---
        epochTitlePrefix = 'fUSI Neurovascular Response';
        subjectLabel     = 'Subject H';

        % --- Intan EEG / Accelerometer (optional) ---
        %intanPath     = 'D:\fusFusiData\April 14th Hobbes Pilot\HobbesPilot300FramesLED\';
        intanPath = '';
        intanBaseName = {''};
        %intanBaseName = {'HobbesPilot300FrameLED_260414_'};
        intanAccelCh  = 1;
        intanFs       = 20000;

    case 'Boltzmann'
        % --- Paths ---

        % --- Per-session definitions (add rows for more sessions) ---
%          sessionDefs(1).bMDir   = '/Users/jpalmer/Downloads/12-10-10-16-26_Boltzmann';
%          sessionDefs(1).bMFile  = 'BMIQ_005';
%          sessionDefs(1).fusiDir = '/Users/jpalmer/Downloads/Hobbes Oct 5th 2026-selected/Acq_2026_10_05_10_58_30/Acq_2026_10_05_10_58_45';

        mrilogFile = fullfile(scriptDir, 'BoltzmannMRI.mat');
    
        % --- B-mode scan geometry (phased array) ---
        faceWidth_mm = 12;
        fovWidth_mm  = 120;
        fovDepth_mm  = 80;
        dynRange_dB  = 60;
    
        % --- Beam / transducer positioning ---
        yCenter    = 40.05e-3;    % +10mm Y from original 30.05e-3
        depthParam = 0.012;
        zDistAdj   = 0.007;       % +10mm Z offset
        angles     = 34:0.5:56;
    
        % --- B-mode orientation corrections ---
        doRot180_US = false;
        flipLR_US   = true;
        flipUD_US   = false;
    
        % --- Initial B-mode overlay offset (WASD/R/E adjustments, default for all sessions) ---
        bm_dx_init   = 9.5;
        bm_dz_init   = 9.5;
        bm_rot_init  = 7.0;
        defaultAngle = 52;
    
        % --- fUSI display options ---
        fusiFrameIdx   = 0;
        fusiXOffset_mm = 0;
    
        % --- fUSI analysis parameters ---
        stim.zscore     = 1.65;
        stim.isSmooth   = false;
        stim.smoothN    = 2;
        stim.actPeriod  = 10;
        stim.restPeriod = 10;
        stim.numAct     = 15;
        %stim.numAct     = 11;
    
        stim.stimFirst  = false;
    
        % --- Transducer face marker ---
        txMarkerThick_mm = 1.5;
        txFaceOffset_mm  = 7;
    
        % --- Epoch plot labels ---
        epochTitlePrefix = 'fUSI Neurovascular Response';
        subjectLabel     = 'Subject H';
    
        % --- Intan EEG / Accelerometer (optional) ---
        %intanPath     = 'D:\fusFusiData\April 14th Hobbes Pilot\HobbesPilot300FramesLED\';
        intanPath = '';
        intanBaseName = {''};
        %intanBaseName = {'HobbesPilot300FrameLED_260414_'};
        intanAccelCh  = 1;
        intanFs       = 20000;

    otherwise
        error('Unknown subject "%s". Use ''Calvin'' or ''Hobbes''.', subject);
end

nSessions = numel(sessionDefs);


fprintf('\n=== Loading subject: %s (%d session(s)) ===\n\n', subject, nSessions);

%% ===================== LOAD MRI VOLUME =====================
fprintf('Loading MRI volume from %s...\n', mrilogFile);
if ~isfile(mrilogFile)
    error('MRI file not found next to FusiOverlay.m:\n  %s', mrilogFile);
end
tmp = load(mrilogFile);
sys = tmp.sys;
clear tmp

% Apply per-subject Z offset to transducer reference distance
zDist_eff = sys.zDist + zDistAdj;
 
%% ===================== OUTPUT GRID =====================
dx_mri_mm = abs(median(diff(sys.ux))) * 1e3;
dz_mri_mm = abs(median(diff(sys.uz))) * 1e3;
pxSize_mm = min(dx_mri_mm, dz_mri_mm);
 
x_out_mm = -(fovWidth_mm/2 + gridPad_mm) : pxSize_mm : (fovWidth_mm/2 + gridPad_mm);
z_out_mm = -(txMarkerThick_mm + gridPad_mm) : pxSize_mm : (fovDepth_mm + gridPad_mm);
[X_OUT, Z_OUT] = meshgrid(x_out_mm, z_out_mm);

%% ===================== PRECOMPUTE MRI SLICES (CACHED) =====================
% The MRI oblique slices and the sagittal beam-path masks are slow to build
% (interp3 + 3-D masks for every angle). They are computed once, saved in
% <saveDir>/cache, and reloaded on later runs. The cache is rebuilt
% automatically if the MRI file or any geometry setting below changes.
nAngles = numel(angles);
hasLGN = isfield(sys, 'leftLgnRoi') && isfield(sys, 'rightLgnRoi');

% Sagittal beam column geometry (needed whether or not the cache is used)
txHeight_sag = 0.003;           % Beam column thickness (m)
txWidth_sag  = 0.015;           % Transducer lateral width (m)

mriInfo = dir(mrilogFile);
cacheSig = struct( ...
    'cacheVersion',     1, ...
    'subject',          subject, ...
    'mriFileBytes',     mriInfo.bytes, ...
    'mriFileDate',      mriInfo.datenum, ...
    'angles',           angles, ...
    'zDistAdj',         zDistAdj, ...
    'depthParam',       depthParam, ...
    'yCenter',          yCenter, ...
    'fovWidth_mm',      fovWidth_mm, ...
    'fovDepth_mm',      fovDepth_mm, ...
    'gridPad_mm',       gridPad_mm, ...
    'txMarkerThick_mm', txMarkerThick_mm, ...
    'txHeight_sag',     txHeight_sag, ...
    'txWidth_sag',      txWidth_sag, ...
    'gridSize',         [numel(z_out_mm), numel(x_out_mm)], ...
    'hasLGN',           hasLGN);
cacheFile = fullfile(cacheDir, sprintf('%s_MRIslices_cache.mat', subject));

cacheHit = false;
if ~forceRecomputeMRI && isfile(cacheFile)
    try
        hdr = load(cacheFile, 'cacheSig');
        if isequal(hdr.cacheSig, cacheSig)
            fprintf('Loading cached MRI slices from %s ...\n', cacheFile);
            C = load(cacheFile);
            mrSlices      = C.mrSlices;
            lgnSlices     = C.lgnSlices;
            faceY_mm      = C.faceY_mm;
            faceZ_mm      = C.faceZ_mm;
            tmpYimg       = C.tmpYimg;
            beamOverlay   = C.beamOverlay;
            sagFacePos    = C.sagFacePos;
            sagSysUy_mm   = C.sagSysUy_mm;
            sagSysUz_mm   = C.sagSysUz_mm;
            sagPivotY_mm  = C.sagPivotY_mm;
            sagPivotZ_mm  = C.sagPivotZ_mm;
            mrNorms = cellfun(@normalizeMriSlice, mrSlices, 'UniformOutput', false);
            clear C hdr
            cacheHit = true;
            fprintf('Done (cache hit).\n');
        else
            fprintf('MRI cache is out of date (MRI file or geometry changed). Recomputing...\n');
        end
    catch ME
        warning('Could not read MRI cache (%s). Recomputing...', ME.message);
        cacheHit = false;
    end
end

if ~cacheHit
    [Ymr, Xmr, Zmr] = meshgrid(sys.uy, sys.ux, sys.uz);
 
    if hasLGN, lgnVolume = sys.leftLgnRoi + sys.rightLgnRoi; end
 
    mrSlices  = cell(nAngles, 1);
    mrNorms   = cell(nAngles, 1);
    lgnSlices = cell(nAngles, 1);
    faceY_mm  = zeros(nAngles, 1);
    faceZ_mm  = zeros(nAngles, 1);
 
    X_mr_base = X_OUT * 1e-3;
 
    fprintf('Precomputing MRI slices for %d angles...\n', nAngles);
    for ai = 1:nAngles
        th = deg2rad(angles(ai));
        zF = zDist_eff + depthParam * cos(th);
        yF = yCenter   - depthParam * sin(th);
        faceY_mm(ai) = yF * 1e3;
        faceZ_mm(ai) = zF * 1e3;
 
        Z_mr_q = zF + Z_OUT * 1e-3 * cos(th);
        Y_mr_q = yF - Z_OUT * 1e-3 * sin(th);
 
        mr = interp3(Ymr, Xmr, Zmr, sys.aImg, Y_mr_q, X_mr_base, Z_mr_q, 'linear', NaN);
        mrSlices{ai} = mr;
 
        mrNorms{ai} = normalizeMriSlice(mr);
 
        if hasLGN
            lgnSlices{ai} = interp3(Ymr, Xmr, Zmr, lgnVolume, Y_mr_q, X_mr_base, Z_mr_q, 'linear', NaN);
        end

        % DEBUG: At theta=53 deg, find depths where the slice ACTUALLY contains LGN
        if angles(ai) == 53 && hasLGN
            lgn = lgnSlices{ai}; lgn(isnan(lgn)) = 0;
            [rows_with_lgn, ~] = find(lgn > 0.1);
            if ~isempty(rows_with_lgn)
                uniq_rows = unique(rows_with_lgn);
                fprintf('DEBUG theta=53: LGN appears at coronal depths:\n');
                for rr = uniq_rows(:)'
                    fprintf('  z_out=%.2fmm -> queries (uy=%.2f, uz=%.2f) mm\n', ...
                        z_out_mm(rr), Y_mr_q(rr,1)*1e3, Z_mr_q(rr,1)*1e3);
                end
            else
                fprintf('DEBUG theta=53: NO LGN voxels found in slice\n');
            end
        end

        % DEBUG: At theta=50 deg, print where the coronal queries at depth=50mm
        if angles(ai) == 50
            [~, depth_row] = min(abs(z_out_mm - 50));
            fprintf('DEBUG theta=50, depth=50mm: queries sys.uy=%.2fmm, sys.uz=%.2fmm (display Z=%.2fmm)\n', ...
                Y_mr_q(depth_row, 1)*1e3, Z_mr_q(depth_row, 1)*1e3, ...
                Z_mr_q(depth_row, 1)*1e3 - zDist_eff*1e3 - depthParam*cos(deg2rad(mean(angles)))*1e3);
        end
        if mod(ai, 5) == 0, fprintf('  %d/%d done\n', ai, nAngles); end
    end
    fprintf('Done.\n');

    %% ===================== PRECOMPUTE SAGITTAL BEAM PATH =====================
    % (Ported from overlayUsFusi subplot 1: sagittal MRI with beam overlay)
    fprintf('Precomputing sagittal beam path data...\n');

    % Transducer geometry for beam column visualization

    % Find midline x-index for sagittal slice extraction
    [~, centerIdxX] = min(abs(sys.ux));

    % Extract sagittal (Y-Z) slice at midline
    tmpYimg = squeeze(sys.aImg(centerIdxX,:,:)).';
    tmpYimg(end-5:end, end-9:end) = max(tmpYimg(:));  % Reference marker

    % Fixed transducer face Z-position (using mean angle as reference)
    refTheta_sag = deg2rad(mean(angles));
    zFace_sag    = depthParam * cos(refTheta_sag) + zDist_eff;
    [~, zFaceIdx] = min(abs(zFace_sag - sys.uz));

    % Precompute per-angle beam overlays and face positions
    beamOverlay = cell(nAngles, 1);
    sagFacePos  = zeros(nAngles, 2);   % [Y, Z] in mm relative to face

    for ai = 1:nAngles
        theta = deg2rad(angles(ai));

        % Beam path mask on 3D MR grid
        Zoff = Zmr - zDist_eff;
        Yoff = Ymr - yCenter;
        txColumn = zeros(size(sys.aImg));
        txColumn(Xmr < txWidth_sag/2 & Xmr > -txWidth_sag/2 & ...
            Zoff < -Yoff/tan(theta) + txHeight_sag/2 & Zoff > -Yoff/tan(theta) - txHeight_sag/2 & ...
            Yoff < -tan(theta)*Zoff + txHeight_sag/2 & Yoff > -tan(theta)*Zoff - txHeight_sag/2) = 1;

        % Extract midline sagittal slice of beam mask
        beamOverlay{ai} = squeeze(txColumn(centerIdxX,:,:)).';

        % Face center position in mm (Z relative to face)
        sagFacePos(ai,:) = [(yCenter - depthParam*sin(theta))*1e3, ...
                            (zDist_eff + depthParam*cos(theta))*1e3 - zFace_sag*1e3];
    end

    sagSysUy_mm   = sys.uy * 1e3;
    sagSysUz_mm   = sys.uz * 1e3 - zFace_sag * 1e3;
    sagPivotY_mm  = yCenter * 1e3;
    sagPivotZ_mm  = zDist_eff * 1e3 - zFace_sag * 1e3;
    fprintf('Done.\n');

    % ---- Save cache for next time ----
    fprintf('Saving MRI slice cache to %s ...\n', cacheFile);
    try
        save(cacheFile, 'cacheSig', 'mrSlices', 'lgnSlices', 'faceY_mm', 'faceZ_mm', ...
            'tmpYimg', 'beamOverlay', 'sagFacePos', 'sagSysUy_mm', 'sagSysUz_mm', ...
            'sagPivotY_mm', 'sagPivotZ_mm', '-v7.3');
        fprintf('Cache saved.\n');
    catch ME
        warning('Could not save MRI cache (%s). It will be recomputed next run.', ME.message);
    end
end  % if ~cacheHit

%% ===================== LOAD ALL SESSIONS (B-mode + fUSI + activation) =====================
% Per-session storage
sess = struct();
fusi_xaxis_mm = [];  fusi_zaxis_mm = [];   % filled in below when fUSI = true
for iSess = 1:nSessions
    sd = sessionDefs(iSess);
    fprintf('\n--- Session %d/%d ---\n', iSess, nSessions);
    fprintf('  B-mode: %s\n', sd.bMDir);
    if fUSI, fprintf('  fUSI:   %s\n', sd.fusiDir); end

    if ~isempty(sd.bMDir) && exist(sd.bMDir, 'dir'), addpath(sd.bMDir); end

    % ---- Load B-mode ----
    fprintf('  Loading B-mode data from %s...\n', sd.bMFile);
    S = load(fullfile(sd.bMDir, sd.bMFile));
    if isfield(S, 'IData') && isfield(S, 'QData')
        IQ = complex(S.IData, S.QData);
    elseif isfield(S, 'bmImage') && isfield(S.bmImage, 'IData') && isfield(S.bmImage, 'QData')
        IQ = complex(S.bmImage.IData, S.bmImage.QData);
    else
        error('Could not find IData/QData in %s.', sd.bMFile);
    end
    clear S

    if exist('bmode', 'file') == 2
        bmImg = bmode(IQ, dynRange_dB);
    else
        warning('bmode() not found. Using simple log-compression.');
        env = abs(IQ);
        env = env / max(env(:));
        env(env < 1e-10) = 1e-10;
        bmImg = 20*log10(env) + dynRange_dB;
        bmImg(bmImg < 0) = 0;
        bmImg = bmImg / dynRange_dB;
    end

    bmMaskRaw = abs(IQ) ~= 0;
    if exist('bwareaopen', 'file') == 2, bmMaskRaw = bwareaopen(bmMaskRaw, 50); end
    clear IQ

    if doRot180_US, bmImg = rot90(bmImg, 2); bmMaskRaw = rot90(bmMaskRaw, 2); end
    if flipLR_US,   bmImg = fliplr(bmImg);   bmMaskRaw = fliplr(bmMaskRaw);   end
    if flipUD_US,   bmImg = flipud(bmImg);   bmMaskRaw = flipud(bmMaskRaw);   end
    if ~isequal(size(bmMaskRaw), size(bmImg)), bmMaskRaw = imresize(bmMaskRaw, size(bmImg), 'nearest'); end

    [Nz_bm, Nx_bm] = size(bmImg);
    fprintf('  B-mode size: %d x %d pixels\n', Nz_bm, Nx_bm);

    % ---- Resample B-mode onto output grid ----
    fprintf('  Resampling B-mode onto output grid...\n');
    zClamp = max(Z_OUT, 0);
    col_frac = (X_OUT / fovWidth_mm + 0.5) * (Nx_bm - 1) + 1;
    row_frac = zClamp / fovDepth_mm * (Nz_bm - 1) + 1;

    bmResampled = interp2(double(bmImg), col_frac, row_frac, 'linear', NaN);
    bmMaskResampled = interp2(double(bmMaskRaw), col_frac, row_frac, 'nearest', 0) > 0.5;

    sectorMask = Z_OUT >= 0 & bmMaskResampled & ...
                 col_frac >= 1 & col_frac <= Nx_bm & ...
                 row_frac >= 1 & row_frac <= Nz_bm & ...
                 ~isnan(bmResampled);
    bmResampled(~sectorMask) = NaN;

    bmNorm = bmResampled;
    vp = sectorMask;
    if any(vp(:))
        lo = min(bmNorm(vp)); hi = max(bmNorm(vp));
        if hi > lo, bmNorm = (bmNorm - lo) / (hi - lo); else, bmNorm(:) = 0; end
    end
    bmNorm(~sectorMask) = 0;
    clear bmImg bmMaskRaw bmMaskResampled

    if fUSI
        % ---- Load fUSI ----
        fprintf('  Loading fUSI data from %s...\n', sd.fusiDir);
        SeqParamPath = dir(fullfile(sd.fusiDir, 'SequenceParam*'));
        SP = load(fullfile(sd.fusiDir, SeqParamPath(1).name));

        fusi_Nz_sess = SP.CParam.PData.Size(1);
        fusi_Nx_sess = SP.CParam.PData.Size(2);
        fusi_Nframes_sess = SP.UserParam.NumLoops + 1;

        if iSess == 1
            fusi_Nz = fusi_Nz_sess;
            fusi_Nx = fusi_Nx_sess;
            fusi_maxDepth_mm = fusi_Nz * SP.CParam.PData.PDelta(3) * SP.UserParam.CenterWavelength * 1e3;
            fusi_z0_mm = SP.CParam.PData.Origin(3) * SP.UserParam.CenterWavelength * 1e3;
            fusi_zaxis_mm = linspace(fusi_z0_mm, fusi_z0_mm + fusi_maxDepth_mm, fusi_Nz);
            fusi_x0_mm = SP.CParam.PData.Origin(1) * SP.UserParam.CenterWavelength * 1e3;
            fusi_xaxis_mm = (1:fusi_Nx) .* SP.CParam.PData.PDelta(1) .* SP.UserParam.CenterWavelength * 1e3;
            fusi_xaxis_mm = fusi_xaxis_mm + fusi_x0_mm + fusiXOffset_mm;
            fprintf('  fUSI grid: %d x %d (origin: x=%.1f mm, z=%.1f mm)\n', ...
                fusi_Nz, fusi_Nx, fusi_x0_mm, fusi_z0_mm);
        else
            if fusi_Nz_sess ~= fusi_Nz || fusi_Nx_sess ~= fusi_Nx
                warning('Session %d grid (%dx%d) differs from session 1 (%dx%d). Using session 1 grid.', ...
                    iSess, fusi_Nz_sess, fusi_Nx_sess, fusi_Nz, fusi_Nx);
            end
        end

        fprintf('  Session %d: %d frames\n', iSess, fusi_Nframes_sess);

        dopFile = fullfile(sd.fusiDir, 'Dop.bin');
        if ~isfile(dopFile)
            fprintf('  Dop.bin not found. Generating it from beamformed data...\n');
            generateDopBin(sd.fusiDir);
        end

        DopSess = readBinFile(dopFile, ...
            [fusi_Nz_sess, fusi_Nx_sess, fusi_Nframes_sess], 'double');

        nExpected = stim.numAct * (stim.restPeriod + stim.actPeriod);
        if fusi_Nframes_sess > nExpected
            fprintf('  Trimming %d -> %d frames (discarding %d extra)\n', ...
                fusi_Nframes_sess, nExpected, fusi_Nframes_sess - nExpected);
            DopSess = DopSess(:, :, 1:nExpected);
            fusi_Nframes_sess = nExpected;
        end

        if doRot180_US, DopSess = rot90(DopSess, 2); end
        if flipLR_US
            for fi = 1:fusi_Nframes_sess, DopSess(:,:,fi) = fliplr(DopSess(:,:,fi)); end
        end
        if flipUD_US
            for fi = 1:fusi_Nframes_sess, DopSess(:,:,fi) = flipud(DopSess(:,:,fi)); end
        end
        clear SP SeqParamPath

        % Build stim pattern for this session
        stimPattern_sess = -ones(1, fusi_Nframes_sess);
        stimOn_sess = (0:stim.numAct-1)' .* (stim.restPeriod + stim.actPeriod) + stim.restPeriod + (1:stim.actPeriod);
        stimOn_sess = stimOn_sess(:);
        stimOn_sess(stimOn_sess > fusi_Nframes_sess) = [];
        stimPattern_sess(stimOn_sess) = 1;
        if stim.stimFirst, stimPattern_sess = -stimPattern_sess; end
        stimPatternNorm_sess = (stimPattern_sess - mean(stimPattern_sess)) ./ std(stimPattern_sess, 1);

        % ---- Resample fUSI (display frame from this session) ----
        fprintf('  Resampling fUSI onto output grid...\n');
        fusi_Nframes_s = size(DopSess, 3);
        if fusiFrameIdx == 0
            dopFrame = mean(DopSess, 3);
        else
            dopFrame = DopSess(:, :, min(fusiFrameIdx, fusi_Nframes_s));
        end

        [Xf, Zf] = meshgrid(fusi_xaxis_mm, fusi_zaxis_mm);
        fusiResampled = interp2(Xf, Zf, dopFrame, X_OUT, Z_OUT, 'linear', NaN);
        fusiMask = ~isnan(fusiResampled) & Z_OUT >= 0;

        fusiNorm = fusiResampled;
        fusiNorm(~fusiMask) = 0;
        vpF = fusiMask;
        if any(vpF(:))
            fusiNorm = sqrt(abs(fusiNorm));
            lo = min(fusiNorm(vpF)); hi = max(fusiNorm(vpF));
            if hi > lo, fusiNorm = (fusiNorm - lo) / (hi - lo); end
        end
        fusiNorm(~fusiMask) = 0;

        % ---- Compute activation map for this session ----
        fprintf('  Computing activation map...\n');
        DopSmooth = DopSess;
        if stim.isSmooth
            fprintf('    Spatial smoothing (kernel=%d)...\n', stim.smoothN);
            for t = 1:fusi_Nframes_s
                DopSmooth(:,:,t) = imfilter(DopSmooth(:,:,t), ones(stim.smoothN)/stim.smoothN^2);
            end
        end

        dopMeanAll = mean(DopSmooth, 3, 'omitnan');
        stdDopAll  = std(DopSmooth, 1, 3, 'omitnan');
        stdDopAll(stdDopAll == 0) = 1;
        dopNormAll = (DopSmooth - dopMeanAll) ./ stdDopAll;

        ActMap_sess = zeros(fusi_Nz, fusi_Nx);
        for i = 1:fusi_Nz
            for ii = 1:fusi_Nx
                tmp = corrcoef(squeeze(dopNormAll(i, ii, :)), stimPatternNorm_sess(:));
                ActMap_sess(i, ii) = tmp(1, 2);
            end
        end
        actThreshold_sess = stim.zscore * std(ActMap_sess(:), 'omitnan') + mean(ActMap_sess(:), 'omitnan');
        actSigMaskNative_sess = ActMap_sess > actThreshold_sess;

        % Resample activation onto output grid
        actResampled = interp2(Xf, Zf, ActMap_sess, X_OUT, Z_OUT, 'linear', NaN);
        actSigResampled = interp2(Xf, Zf, double(actSigMaskNative_sess), X_OUT, Z_OUT, 'nearest', 0) > 0.5;
        actMaskGrid = ~isnan(actResampled) & Z_OUT >= 0;

        actNormGrid = actResampled;
        actNormGrid(actNormGrid < 0) = 0;
        actNormGrid(~actMaskGrid) = 0;
        hiAct = max(actNormGrid(:));
        if hiAct > 0, actNormGrid = actNormGrid / hiAct; end

        corrRawGrid = actResampled;
        corrRawGrid(~actMaskGrid) = 0;

        fprintf('  Activation map: %d significant pixels\n', sum(actSigMaskNative_sess(:)));
        clear DopSmooth dopNormAll dopMeanAll stdDopAll dopFrame fusiResampled actResampled
    else
        % B-mode / MRI alignment only: no fUSI data is read or computed
        fusiNorm = [];  fusiMask = [];  DopSess = [];
        actNormGrid = [];  corrRawGrid = [];  actSigResampled = [];  actMaskGrid = [];
        ActMap_sess = [];  actThreshold_sess = [];  stimPatternNorm_sess = [];
    end

    % ---- Store this session's data ----
    sess(iSess).bmResampled  = bmResampled;
    sess(iSess).bmNorm       = bmNorm;
    sess(iSess).sectorMask   = sectorMask;
    sess(iSess).fusiNorm     = fusiNorm;
    sess(iSess).fusiMask     = fusiMask;
    sess(iSess).DopAll       = DopSess;
    sess(iSess).actNormGrid  = actNormGrid;
    sess(iSess).corrRawGrid  = corrRawGrid;
    sess(iSess).actSigResampled = actSigResampled;
    sess(iSess).actMaskGrid  = actMaskGrid;
    sess(iSess).ActMap       = ActMap_sess;
    sess(iSess).actThreshold = actThreshold_sess;
    sess(iSess).stimPatternNorm = stimPatternNorm_sess;
    sess(iSess).bm_dx_init   = bm_dx_init;   % default; overridden below if saved
    sess(iSess).bm_dz_init   = bm_dz_init;
    sess(iSess).bm_rot_init  = bm_rot_init;
    sess(iSess).overlayState = 0;
    sess(iSess).roiResampled = [];
    sess(iSess).roiIsManual  = false;
    sess(iSess).savedCurrentIdx = [];

    % Alignment file. With fUSI on it lives in the fUSI folder, as before. With fUSI
    % off it lives next to the B-mode file (one file per B-mode capture). When fUSI is
    % on and the fUSI folder has none yet, an alignment saved by a B-mode-only run is used.
    bmAlignFile = fullfile(sd.bMDir, ['bmode_alignment_' sd.bMFile '.mat']);
    if fUSI
        alignFile = fullfile(sd.fusiDir, 'bmode_alignment.mat');
        if isfile(alignFile)
            loadAlignFile = alignFile;
        else
            loadAlignFile = bmAlignFile;
        end
    else
        alignFile = bmAlignFile;
        loadAlignFile = bmAlignFile;
    end
    sess(iSess).alignFile = alignFile;   % where this session's alignment is saved

    % Load saved alignment from disk if available
    if isfile(loadAlignFile)
        saved = load(loadAlignFile);
        if isfield(saved, 'bm_dx'),   sess(iSess).bm_dx_init  = saved.bm_dx;   end
        if isfield(saved, 'bm_dz'),   sess(iSess).bm_dz_init  = saved.bm_dz;   end
        if isfield(saved, 'bm_rot'),  sess(iSess).bm_rot_init = saved.bm_rot;  end
        if isfield(saved, 'angleIdx'), sess(iSess).savedCurrentIdx = saved.angleIdx; end
        if isfield(saved, 'overlayState'), sess(iSess).overlayState = saved.overlayState; end
        fprintf('  Loaded saved alignment: dx=%.1f dz=%.1f rot=%.1f\n', ...
            sess(iSess).bm_dx_init, sess(iSess).bm_dz_init, sess(iSess).bm_rot_init);
    end

    clear bmResampled bmNorm sectorMask fusiNorm fusiMask DopSess
    clear actNormGrid corrRawGrid actSigResampled actMaskGrid ActMap_sess
end
fprintf('\nAll %d sessions loaded.\n', nSessions);

% Session 1 frame count (used by Intan alignment below)
if fUSI
    fusi_Nframes = size(sess(1).DopAll, 3);
else
    fusi_Nframes = 0;
end

%% ===================== LOAD INTAN EEG + ACCEL (OPTIONAL) =====================
hasIntan = fUSI && ~isempty(intanPath) && ~isempty(intanBaseName);
if hasIntan
    fprintf('Loading Intan EEG + analog data from %s...\n', intanPath);
    [eegCollect, eegCollectRaw, digTrigFull, analogFull] = concatIntanFUSILights(intanPath, intanBaseName);

    % Extract accelerometer from full analog channels
    accelFull = analogFull(intanAccelCh, :);

    % --- Align Intan and fUSI time axes ---
    % Session is 300 s total. fUSI acquires 1 frame/s starting 10 s before
    % the first rising edge of the digital trigger.
    trigRiseIdx = find(diff(digTrigFull) == 1) + 1;
    fprintf('  Found %d trigger rising edges (fUSI has %d frames)\n', length(trigRiseIdx), fusi_Nframes);

    % Intan time axis (seconds) with t=0 at first trigger rising edge
    intanTimeAxis = ((1:length(digTrigFull)) - trigRiseIdx(1)) / intanFs;

    % fUSI frame times: 1 Hz, starting 10 s before first trigger
    fusiFrameTimes = -10 + (0:fusi_Nframes-1);
    nAligned = fusi_Nframes;

    fprintf('  Intan data: %.1f s total | fUSI session: %d to %d s\n', ...
        intanTimeAxis(end) - intanTimeAxis(1), fusiFrameTimes(1), fusiFrameTimes(end));

    % --- Quantify trigger epoch durations (on vs off) ---
    % The digital trigger oscillates during "on" epochs and stays low during
    % "off" epochs. Detect epoch boundaries by finding gaps between rising
    % edges that exceed a threshold (0.5 s covers any plausible pulse rate).
    trigRiseTimes = (trigRiseIdx(:) - trigRiseIdx(1)) / intanFs;  % seconds from first edge
    interPulseGaps = diff(trigRiseTimes);
    gapThreshold = 0.5;  % seconds; gaps longer than this mark epoch boundaries
    boundaryIdx = find(interPulseGaps > gapThreshold);
    boundaryIdx = boundaryIdx(:);  % ensure column vector

    % Build epoch table: each "on" epoch runs from one group of pulses to the next gap
    % boundaryIdx(k) = last pulse index before gap k
    onStarts  = [1; boundaryIdx + 1];            % first pulse index of each on-epoch
    onEnds    = [boundaryIdx; length(trigRiseTimes)];  % last pulse index of each on-epoch
    nEpochs   = length(onStarts);

    onDurations  = trigRiseTimes(onEnds) - trigRiseTimes(onStarts);  % on duration
    offDurations = zeros(nEpochs, 1);
    for ei = 1:nEpochs - 1
        offDurations(ei) = trigRiseTimes(onStarts(ei+1)) - trigRiseTimes(onEnds(ei));
    end
    offDurations(end) = NaN;  % no off period after last epoch

    fprintf('\n  --- Trigger Epoch Timing (expected: %.1f s on, %.1f s off) ---\n', ...
        stim.actPeriod, stim.restPeriod);
    fprintf('  %5s  %10s  %10s  %10s  %10s\n', 'Epoch', 'On (s)', 'Drift On', 'Off (s)', 'Drift Off');
    for ei = 1:nEpochs
        driftOn  = onDurations(ei) - stim.actPeriod;
        if ~isnan(offDurations(ei))
            driftOff = offDurations(ei) - stim.restPeriod;
            fprintf('  %5d  %10.3f  %+10.3f  %10.3f  %+10.3f\n', ...
                ei, onDurations(ei), driftOn, offDurations(ei), driftOff);
        else
            fprintf('  %5d  %10.3f  %+10.3f  %10s  %10s\n', ...
                ei, onDurations(ei), driftOn, '-', '-');
        end
    end

    cumDriftOn  = sum(onDurations) - nEpochs * stim.actPeriod;
    cumDriftOff = sum(offDurations(1:end-1)) - (nEpochs-1) * stim.restPeriod;
    fprintf('  -------------------------------------------------------\n');
    fprintf('  Cumulative drift:  On: %+.3f s   Off: %+.3f s\n', cumDriftOn, cumDriftOff);
    fprintf('  Mean epoch:        On: %.3f s    Off: %.3f s\n', ...
        mean(onDurations), mean(offDurations(1:end-1)));
    fprintf('  Total epochs: %d (expected: %d)\n\n', nEpochs, stim.numAct);

    % Store for potential plotting
    trigEpochInfo.onStarts     = trigRiseTimes(onStarts);
    trigEpochInfo.onEnds       = trigRiseTimes(onEnds);
    trigEpochInfo.onDurations  = onDurations;
    trigEpochInfo.offDurations = offDurations;
    trigEpochInfo.nEpochs      = nEpochs;

    % Check if EEG data was found
    hasEeg = ~isempty(eegCollect) && ~isempty(eegCollect{1});
else
    if fUSI
        fprintf('No Intan data path specified, skipping EEG/accel loading.\n');
    end
    eegCollect     = {};
    eegCollectRaw  = {};
    digTrigFull    = [];
    accelFull      = [];
    intanTimeAxis  = [];
    fusiFrameTimes = [];
    nAligned       = 0;
    hasEeg         = false;
    trigEpochInfo  = [];
end

%% ===================== PACK INTO animalStore =====================
as = struct();
as.subject        = subject;
as.angles         = angles;
as.defaultAngle   = defaultAngle;
as.nSessions      = nSessions;
as.currentSession = 1;
as.sessions       = sess;              % per-session B-mode, fUSI, activation, alignment
as.fusi_xaxis_mm  = fusi_xaxis_mm;
as.fusi_zaxis_mm  = fusi_zaxis_mm;
as.stim           = stim;
as.mrSlices       = mrSlices;
as.mrNorms        = mrNorms;
as.lgnSlices      = lgnSlices;
as.hasLGN         = hasLGN;
as.faceY_mm       = faceY_mm;
as.faceZ_mm       = faceZ_mm;
as.hasIntan       = hasIntan;
as.hasEeg         = hasEeg;
as.eegCollect     = eegCollect;
as.eegCollectRaw  = eegCollectRaw;
as.digTrigFull    = digTrigFull;
as.accelFull      = accelFull;
as.intanTimeAxis  = intanTimeAxis;
as.fusiFrameTimes = fusiFrameTimes;
as.nAligned       = nAligned;
as.intanFs        = intanFs;
as.trigEpochInfo   = trigEpochInfo;
as.tmpYimg        = tmpYimg;
% Sagittal plane browsing (Left/Right keys): keep the MRI volume so any plane can be shown
[~, sagCenterIdx] = min(abs(sys.ux));
as.sagVol         = sys.aImg;                  % [Nx x Ny x Nz]; plane ix is aImg(ix,:,:)
as.sagUx_mm       = sys.ux(:)' * 1e3;          % lateral position of each plane (mm)
as.sagIdx         = sagCenterIdx;              % start on the midline plane
as.sagStepVox     = max(1, round(1 / dx_mri_mm));   % planes per key press (about 1 mm); Shift = 1 plane
as.beamOverlay    = beamOverlay;
as.sagFacePos     = sagFacePos;
as.sagSysUy_mm    = sagSysUy_mm;
as.sagSysUz_mm    = sagSysUz_mm;
as.sagPivotY_mm   = sagPivotY_mm;
as.sagPivotZ_mm   = sagPivotZ_mm;
as.txWidth_sag_mm = txWidth_sag * 1e3;
as.angleRange     = [angles(1), angles(end)];
as.faceWidth_mm   = faceWidth_mm;
as.txMarkerThick_mm = txMarkerThick_mm;
as.txFaceOffset_mm  = txFaceOffset_mm;
as.fovWidth_mm    = fovWidth_mm;
as.fovDepth_mm    = fovDepth_mm;
as.x_out_mm       = x_out_mm;
as.z_out_mm       = z_out_mm;
as.pxSize_mm      = pxSize_mm;
as.epochData      = [];
as.epochTitlePrefix = epochTitlePrefix;
as.subjectLabel     = subjectLabel;
animalStore{iAnimal} = as;

fprintf('\n=== %s data loaded and stored ===\n\n', subject);
end  % for iAnimal

% Free workspace copies of large per-animal data (kept in animalStore)
clear sess mrSlices mrNorms lgnSlices sys
clear beamOverlay tmpYimg


%% ===================== INTERACTIVE DISPLAY =====================
ud = struct();
ud.animalStore   = animalStore;
ud.currentAnimal = 1;
ud.saveDir       = saveDir;   % local functions can't see script variables
ud.fUSI          = fUSI;      % false = B-mode / MRI alignment only

% Shared display settings
ud.fusiAlpha        = fusiAlpha;
ud.fusiColormap     = fusiColormap;
ud.overlayAlpha     = overlayAlpha;
ud.overlayMode      = overlayMode;
ud.overlayColormap  = overlayColormap;
ud.fusionGamma      = fusionGamma;
ud.checkerBlockSize = checkerBlockSize;

% Load first animal (session 1)
ud = loadAnimalIntoUd(ud, 1);

fig = figure('Name', figureTitle(ud), ...
    'NumberTitle', 'off', 'Color', 'k', 'Position', [20 40 1860 960], ...
    'WindowKeyPressFcn', @(src, evt) keyHandler(src, evt), ...
    'CloseRequestFcn', @(src, ~) closeAndSave(src));

fig.UserData = ud;
updateDisplay(fig);

%% ===================== LOAD ANIMAL INTO UD =====================
function ud = loadAnimalIntoUd(ud, animalIdx)
    as = ud.animalStore{animalIdx};
    ud.currentAnimal = animalIdx;
    ud.subject      = as.subject;
    ud.angles       = as.angles;
    ud.x_out_mm     = as.x_out_mm;
    ud.z_out_mm     = as.z_out_mm;
    ud.pxSize_mm    = as.pxSize_mm;
    ud.nSessions    = as.nSessions;
    ud.currentSession = as.currentSession;
    ud.mrSlices     = as.mrSlices;
    ud.mrNorms      = as.mrNorms;
    ud.lgnSlices    = as.lgnSlices;
    ud.hasLGN       = as.hasLGN;
    ud.faceY_mm     = as.faceY_mm;
    ud.faceZ_mm     = as.faceZ_mm;
    ud.faceWidth_mm = as.faceWidth_mm;
    ud.txMarkerThick_mm = as.txMarkerThick_mm;
    ud.txFaceOffset_mm  = as.txFaceOffset_mm;
    ud.fovWidth_mm  = as.fovWidth_mm;
    ud.fovDepth_mm  = as.fovDepth_mm;
    ud.fusi_xaxis_mm = as.fusi_xaxis_mm;
    ud.fusi_zaxis_mm = as.fusi_zaxis_mm;
    ud.stim         = as.stim;
    ud.hasIntan     = as.hasIntan;
    ud.hasEeg       = as.hasEeg;
    ud.eegCollect   = as.eegCollect;
    ud.eegCollectRaw = as.eegCollectRaw;
    ud.digTrigFull  = as.digTrigFull;
    ud.accelFull    = as.accelFull;
    ud.intanTimeAxis = as.intanTimeAxis;
    ud.fusiFrameTimes = as.fusiFrameTimes;
    ud.nAligned     = as.nAligned;
    ud.intanFs      = as.intanFs;
    ud.trigEpochInfo = as.trigEpochInfo;
    ud.tmpYimg      = as.tmpYimg;
    ud.sagVol       = as.sagVol;
    ud.sagUx_mm     = as.sagUx_mm;
    ud.sagIdx       = as.sagIdx;
    ud.sagStepVox   = as.sagStepVox;
    ud.beamOverlay  = as.beamOverlay;
    ud.sagFacePos   = as.sagFacePos;
    ud.sagSysUy_mm  = as.sagSysUy_mm;
    ud.sagSysUz_mm  = as.sagSysUz_mm;
    ud.sagPivotY_mm = as.sagPivotY_mm;
    ud.sagPivotZ_mm = as.sagPivotZ_mm;
    ud.txWidth_sag_mm    = as.txWidth_sag_mm;
    ud.angleRange        = as.angleRange;
    ud.epochTitlePrefix  = as.epochTitlePrefix;
    ud.subjectLabel      = as.subjectLabel;
    ud.sessions          = as.sessions;  % carry full session array for analysis

    % Load current session's display data
    ud = loadSessionIntoUd(ud, ud.currentSession);
end

%% ===================== LOAD SESSION INTO UD (display fields) =====================
function ud = loadSessionIntoUd(ud, sessIdx)
    ss = ud.sessions(sessIdx);
    ud.currentSession = sessIdx;
    ud.bmResampled  = ss.bmResampled;
    ud.bmNorm       = ss.bmNorm;
    ud.sectorMask   = ss.sectorMask;
    ud.fusiNorm     = ss.fusiNorm;
    ud.fusiMask     = ss.fusiMask;
    ud.DopAll       = ss.DopAll;
    ud.actNorm      = ss.actNormGrid;
    ud.corrRawGrid  = ss.corrRawGrid;
    ud.actSigMask   = ss.actSigResampled & ss.actMaskGrid;
    ud.actMask      = ss.actMaskGrid;
    ud.ActMap       = ss.ActMap;
    ud.actThreshold = ss.actThreshold;
    ud.stimPatternNorm = ss.stimPatternNorm;
    ud.roiResampled = ss.roiResampled;
    ud.roiIsManual  = ss.roiIsManual;
    ud.overlayState = ss.overlayState;

    % Load alignment offsets
    ud.bm_dx_mm     = ss.bm_dx_init;
    ud.bm_dz_mm     = ss.bm_dz_init;
    ud.bm_rot_deg   = ss.bm_rot_init;

    % Load angle index
    if ~isempty(ss.savedCurrentIdx) && ss.savedCurrentIdx <= numel(ud.angles)
        ud.currentIdx = ss.savedCurrentIdx;
    else
        ud.currentIdx = find(ud.angles == ud.animalStore{ud.currentAnimal}.defaultAngle, 1, 'first');
        if isempty(ud.currentIdx), ud.currentIdx = round(numel(ud.angles)/2); end
    end
end

%% ===================== SAVE UD STATE BACK TO animalStore =====================
function ud = saveAnimalState(ud)
    idx = ud.currentAnimal;
    si  = ud.currentSession;

    % Save per-session alignment and display state
    ud.animalStore{idx}.sessions(si).bm_dx_init     = ud.bm_dx_mm;
    ud.animalStore{idx}.sessions(si).bm_dz_init     = ud.bm_dz_mm;
    ud.animalStore{idx}.sessions(si).bm_rot_init    = ud.bm_rot_deg;
    ud.animalStore{idx}.sessions(si).overlayState    = ud.overlayState;
    ud.animalStore{idx}.sessions(si).roiResampled    = ud.roiResampled;
    ud.animalStore{idx}.sessions(si).roiIsManual     = ud.roiIsManual;
    ud.animalStore{idx}.sessions(si).savedCurrentIdx = ud.currentIdx;

    % Also update the session array in ud so it stays in sync
    ud.sessions(si) = ud.animalStore{idx}.sessions(si);

    % Save which session was active
    ud.animalStore{idx}.currentSession = ud.currentSession;
    ud.animalStore{idx}.defaultAngle   = ud.angles(ud.currentIdx);
    ud.animalStore{idx}.sagIdx         = ud.sagIdx;

    % Persist alignment to disk
    saveAlignmentToDisk(ud.animalStore{idx}.sessions(si));
end

%% ===================== PERSIST ALIGNMENT TO DISK =====================
function saveAlignmentToDisk(sessStruct)
    if ~isfield(sessStruct, 'alignFile') || isempty(sessStruct.alignFile), return; end
    alignFile = sessStruct.alignFile;
    bm_dx = sessStruct.bm_dx_init; %#ok<NASGU>
    bm_dz = sessStruct.bm_dz_init; %#ok<NASGU>
    bm_rot = sessStruct.bm_rot_init; %#ok<NASGU>
    angleIdx = sessStruct.savedCurrentIdx; %#ok<NASGU>
    overlayState = sessStruct.overlayState; %#ok<NASGU>
    try
        save(alignFile, 'bm_dx', 'bm_dz', 'bm_rot', 'angleIdx', 'overlayState');
    catch ME
        fprintf('  Warning: could not save alignment to %s: %s\n', alignFile, ME.message);
    end
end

%% ===================== CLOSE AND SAVE ALL ALIGNMENTS =====================
function closeAndSave(fig)
    ud = fig.UserData;
    if ~isempty(ud) && isstruct(ud) && isfield(ud, 'animalStore')
        % Save current session state
        ud = saveAnimalState(ud);
        % Save all sessions for all animals
        for iA = 1:length(ud.animalStore)
            as = ud.animalStore{iA};
            for iS = 1:as.nSessions
                saveAlignmentToDisk(as.sessions(iS));
            end
        end
        fprintf('All alignments saved to disk.\n');
    end
    delete(fig);
end

function updateDisplay(fig)
    ud  = fig.UserData;
    ai  = ud.currentIdx;
    ang = ud.angles(ai);
 
    % Redraw everything except the Controls panel. Deleting and rebuilding the buttons on
    % every redraw (clf) left the figure without keyboard focus after a button click.
    kids = fig.Children;
    for k = numel(kids):-1:1
        if ~strcmp(get(kids(k), 'Tag'), 'ControlsPanel')
            delete(kids(k));
        end
    end
    buildControlPanel(fig);   % creates the panel the first time, then only refreshes it

    tl = tiledlayout(fig, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl.OuterPosition = [0.005 0.01 0.805 0.98];

    whiteAxis = @(axH) set(axH, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');

    mrSlice = ud.mrSlices{ai};
    mrNorm  = ud.mrNorms{ai};

    % --- Panel 0: Sagittal beam path (from overlayUsFusi) ---
    ax0 = nexttile(tl, 1, [1 2]);
    % Compute face position first so we can offset everything
    theta_rad = deg2rad(ang);
    fY = ud.sagFacePos(ai,1);
    fZ = ud.sagFacePos(ai,2);
    % Offset: shift all Y coords so the transducer face is at x = 0
    yOff = fY;
    sagImg = squeeze(ud.sagVol(ud.sagIdx, :, :)).';           % [Nz x Ny] sagittal plane
    sagImg(end-5:end, end-9:end) = max(sagImg(:));            % reference marker (as before)
    imagesc(ax0, ud.sagSysUy_mm - yOff, ud.sagSysUz_mm, sagImg);
    colormap(ax0, gray); hold(ax0, 'on');
    % Transducer face patch (fY - yOff = 0)
    perpDir = [cos(theta_rad), sin(theta_rad)];
    beamDir = [-sin(theta_rad), cos(theta_rad)];
    txHalfW = 3;  txThick = 1.5;
    fY0 = 0;  % face is at origin after offset
    corners = [fY0 + txHalfW*perpDir(1),                        fZ + txHalfW*perpDir(2);
               fY0 + txHalfW*perpDir(1) - txThick*beamDir(1),  fZ + txHalfW*perpDir(2) - txThick*beamDir(2);
               fY0 - txHalfW*perpDir(1) - txThick*beamDir(1),  fZ - txHalfW*perpDir(2) - txThick*beamDir(2);
               fY0 - txHalfW*perpDir(1),                        fZ - txHalfW*perpDir(2)];
    patch(ax0, corners(:,1), corners(:,2), 'g', 'EdgeColor', 'g', 'FaceAlpha', 0.5);

    % --- Angle range wedge (translucent green, anchored at transducer face) ---
    wedgeR = 50;  % mm depth from face (matches beam line)
    % Near arc: face positions across all available angles (offset)
    faceArcY = ud.sagFacePos(:,1)' - yOff;
    faceArcZ = ud.sagFacePos(:,2)';
    % Far arc: 50mm along beam direction from each face position
    allThetas = deg2rad(ud.angles);
    farArcY = faceArcY + wedgeR * (-sin(allThetas));
    farArcZ = faceArcZ + wedgeR * cos(allThetas);
    % Polygon: face arc forward, far arc backward
    wedgeY = [faceArcY, fliplr(farArcY)];
    wedgeZ = [faceArcZ, fliplr(farArcZ)];
    fill(ax0, wedgeY, wedgeZ, 'g', 'FaceAlpha', 0.12, 'EdgeColor', 'g', ...
        'LineWidth', 1.2, 'LineStyle', '-');

    % --- Current beam path line (50mm from face) ---
    beamLineR = 50;  % mm
    lineEndY = fY0 + beamLineR * beamDir(1);
    lineEndZ = fZ + beamLineR * beamDir(2);
    plot(ax0, [fY0, lineEndY], [fZ, lineEndZ], 'y-', 'LineWidth', 1.8);

    hold(ax0, 'off');
    axis(ax0, 'image'); whiteAxis(ax0);
    ylim(ax0, [-40 60]);
    xlabel(ax0, 'Distance from face (mm)'); ylabel(ax0, 'Depth from face (mm)');
    title(ax0, {sprintf('Beam Path (%.1f\\circ)  [range: %.1f\\circ-%.1f\\circ]', ...
        ang, ud.angleRange(1), ud.angleRange(2)), ...
        sprintf('Sagittal plane x = %+.1f mm (beam at x = 0)', ud.sagUx_mm(ud.sagIdx))}, 'Color', 'w');
 
    if ud.bm_rot_deg ~= 0
        bmNorm   = imrotate(ud.bmNorm, ud.bm_rot_deg, 'bilinear', 'crop');
        mask     = imrotate(double(ud.sectorMask), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
        bmDisp   = imrotate(ud.bmResampled, ud.bm_rot_deg, 'bilinear', 'crop');
        bmDisp(isnan(bmDisp)) = 0;
        if ud.fUSI
            fusiNorm = imrotate(ud.fusiNorm, ud.bm_rot_deg, 'bilinear', 'crop');
            fusiMask = imrotate(double(ud.fusiMask), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
            actNorm    = imrotate(ud.actNorm, ud.bm_rot_deg, 'bilinear', 'crop');
            corrRaw    = imrotate(ud.corrRawGrid, ud.bm_rot_deg, 'bilinear', 'crop');
            actSigMask = imrotate(double(ud.actSigMask), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
            actMask    = imrotate(double(ud.actMask), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
            if ~isempty(ud.roiResampled)
                if ud.roiIsManual
                    roiResamp = ud.roiResampled;  % already in rotated frame
                else
                    roiResamp = imrotate(double(ud.roiResampled), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
                end
            else
                roiResamp = [];
            end
        end
    else
        bmNorm   = ud.bmNorm;
        mask     = ud.sectorMask;
        bmDisp   = ud.bmResampled; bmDisp(isnan(bmDisp)) = 0;
        if ud.fUSI
            fusiNorm = ud.fusiNorm;
            fusiMask = ud.fusiMask;
            actNorm    = ud.actNorm;
            corrRaw    = ud.corrRawGrid;
            actSigMask = ud.actSigMask;
            actMask    = ud.actMask;
            roiResamp  = ud.roiResampled;
        end
    end
 
    overlayState = ud.overlayState;  % 0=Doppler, 1=Activation, 2=Correlation
 
    bm_x = ud.x_out_mm + ud.bm_dx_mm;
    bm_z = ud.z_out_mm + ud.bm_dz_mm;
 
    % Axis limits: union of MRI original FOV and shifted B-mode FOV, plus margin
    margin = 3;  % mm
    mri_xRange = [-ud.fovWidth_mm/2, ud.fovWidth_mm/2];
    mri_zRange = [-ud.txMarkerThick_mm, ud.fovDepth_mm];
    bm_xRange  = [-ud.fovWidth_mm/2 + ud.bm_dx_mm, ud.fovWidth_mm/2 + ud.bm_dx_mm];
    bm_zRange  = [-ud.txMarkerThick_mm + ud.bm_dz_mm, ud.fovDepth_mm + ud.bm_dz_mm];
    xLim = [min(mri_xRange(1), bm_xRange(1)) - margin, max(mri_xRange(2), bm_xRange(2)) + margin];
    yLim = [min(mri_zRange(1), bm_zRange(1)) - margin, max(mri_zRange(2), bm_zRange(2)) + margin];
 
    % --- Panel 1: MRI + LGN ---
    ax1 = nexttile(tl, 3, [1 2]);
    imagesc(ax1, ud.x_out_mm, ud.z_out_mm, mrSlice);
    colormap(ax1, gray); axis(ax1, 'image'); whiteAxis(ax1);
    hold(ax1, 'on');
    if ud.hasLGN
        lgn = ud.lgnSlices{ai}; lgn(isnan(lgn)) = 0;
        lgnMax = max(lgn(:));
        if lgnMax > 0
            lgnNorm = lgn / lgnMax;
            lgnIdx = max(1, min(256, round(lgnNorm * 255) + 1));
            lgnRGB = ind2rgb(lgnIdx, winter(256));
            hLGN = image(ax1, ud.x_out_mm, ud.z_out_mm, lgnRGB);
            hLGN.AlphaData = double(lgnNorm > 0.1);
        end
        title(ax1, sprintf('MRI + LGN (%.1f\\circ)', ang), 'Color', 'w');
    else
        title(ax1, sprintf('MRI Slice (%.1f\\circ)', ang), 'Color', 'w');
    end
    hold(ax1, 'off');
    xlim(ax1, xLim); ylim(ax1, yLim);
    xlabel(ax1, 'Distance from midline (mm)'); ylabel(ax1, 'Depth (mm)');
 
    % --- Panel 2: B-mode ---
    ax2 = nexttile(tl, 5, [1 2]);
    imagesc(ax2, bm_x, bm_z, bmDisp);
    colormap(ax2, gray); axis(ax2, 'image'); whiteAxis(ax2);
    xlim(ax2, xLim); ylim(ax2, yLim);
    xlabel(ax2, 'Distance from midline (mm)'); ylabel(ax2, 'Depth (mm)');
    title(ax2, 'B-mode (adjusted)', 'Color', 'w');
 
    % --- Panel 3: B-mode / MRI overlay ---
    if ud.fUSI
        ax3 = nexttile(tl, 7, [1 2]);
    else
        ax3 = nexttile(tl, 8, [1 4]);   % no fUSI panels: give the B-mode / MRI overlay the whole bottom row
    end
    imagesc(ax3, ud.x_out_mm, ud.z_out_mm, mrNorm);
    colormap(ax3, gray); hold(ax3, 'on');
    switch ud.overlayMode
        case 'fusion'
            usG = bmNorm .^ ud.fusionGamma;
            magRGB = cat(3, usG, zeros(size(usG)), usG);
            hUS = image(ax3, bm_x, bm_z, magRGB);
            hUS.AlphaData = double(mask) .* min(usG * 1.5, 1);
        case 'checkerboard'
            [Nr, Nc] = size(bmNorm); [rr, cc] = ndgrid(1:Nr, 1:Nc);
            bs = ud.checkerBlockSize;
            checker = mod(floor((rr-1)/bs) + floor((cc-1)/bs), 2) == 1;
            hUS = image(ax3, bm_x, bm_z, repmat(bmNorm, [1 1 3]));
            hUS.AlphaData = double(mask & checker);
        case 'alpha'
            cmap = feval(ud.overlayColormap, 256);
            bmIdx = max(1, min(256, round(bmNorm * 255) + 1));
            hUS = image(ax3, bm_x, bm_z, ind2rgb(bmIdx, cmap));
            hUS.AlphaData = double(mask) .* ud.overlayAlpha .* (bmNorm .^ 0.8);
    end
    hold(ax3, 'off');
    axis(ax3, 'image'); whiteAxis(ax3);
    xlim(ax3, xLim);
    ylim(ax3, yLim);
    xlabel(ax3, 'Distance from midline (mm)'); ylabel(ax3, 'Depth (mm)');
    adjStr = sprintf('dx=%.1f  dz=%.1f mm  rot=%.1f', ud.bm_dx_mm, ud.bm_dz_mm, ud.bm_rot_deg);
    title(ax3, sprintf('B-mode overlay (%.1f\\circ) | %s', ang, adjStr), 'Color', 'w');
 
    if ud.fUSI
        % --- Panel 4: B-mode + fUSI/Activation overlay ---
        ax4 = nexttile(tl, 9, [1 2]);
        ax4.Tag = 'fusi_panel';
        imagesc(ax4, bm_x, bm_z, bmDisp);
        colormap(ax4, gray); hold(ax4, 'on');
 
        if overlayState == 1
            % Significant activation over B-mode
            act_cmap = hot(256);
            aIdx = max(1, min(256, round(actNorm * 255) + 1));
            aRGB = ind2rgb(aIdx, act_cmap);
            aAlpha = double(actSigMask) .* 0.8;
            hAct = image(ax4, bm_x, bm_z, aRGB);
            hAct.AlphaData = aAlpha;
            % ROI outline
            if any(roiResamp(:))
                [B_roi, ~] = bwboundaries(roiResamp);
                for k = 1:length(B_roi)
                    bnd = B_roi{k};
                    bnd_x = interp1(1:length(ud.x_out_mm), ud.x_out_mm, bnd(:,2), 'linear', 'extrap') + ud.bm_dx_mm;
                    bnd_z = interp1(1:length(ud.z_out_mm), ud.z_out_mm, bnd(:,1), 'linear', 'extrap') + ud.bm_dz_mm;
                    plot(ax4, bnd_x, bnd_z, 'g-', 'LineWidth', 1.5);
                end
            end
            title(ax4, 'B-mode + Activation', 'Color', 'w');
        elseif overlayState == 2
            % Full correlation map over B-mode (diverging blue-white-red)
            nC = 256;
            halfN = nC / 2;
            blue2white = [linspace(0,1,halfN)', linspace(0,1,halfN)', ones(halfN,1)];
            white2red  = [ones(halfN,1), linspace(1,0,halfN)', linspace(1,0,halfN)'];
            divCmap = [blue2white; white2red];
            cLim = max(abs(corrRaw(actMask)));
            if cLim == 0, cLim = 1; end
            corrScaled = (corrRaw / cLim + 1) / 2;  % map [-cLim,cLim] to [0,1]
            cIdx = max(1, min(nC, round(corrScaled * (nC - 1)) + 1));
            cRGB = ind2rgb(cIdx, divCmap);
            cAlpha = double(actMask) .* 0.7;
            hCorr = image(ax4, bm_x, bm_z, cRGB);
            hCorr.AlphaData = cAlpha;
            title(ax4, sprintf('B-mode + Correlation (r: %.2f to %.2f)', -cLim, cLim), 'Color', 'w');
        else
            % Doppler over B-mode
            fusi_cmap = feval(ud.fusiColormap, 256);
            fIdx = max(1, min(256, round(fusiNorm * 255) + 1));
            fRGB = ind2rgb(fIdx, fusi_cmap);
            fAlpha = double(fusiMask) .* ud.fusiAlpha .* (fusiNorm .^ 0.6);
            hFbm = image(ax4, bm_x, bm_z, fRGB);
            hFbm.AlphaData = fAlpha;
            title(ax4, 'B-mode + fUSI', 'Color', 'w');
        end
 
        hold(ax4, 'off');
        axis(ax4, 'image'); whiteAxis(ax4);
        xlim(ax4, xLim);
        ylim(ax4, yLim);
        xlabel(ax4, 'Distance from midline (mm)'); ylabel(ax4, 'Depth (mm)');
 
        % --- Panel 5: fUSI/Activation / MRI overlay ---
        ax5 = nexttile(tl, 11, [1 2]);
        imagesc(ax5, ud.x_out_mm, ud.z_out_mm, mrNorm);
        colormap(ax5, gray); hold(ax5, 'on');
 
        if overlayState == 1
            % Significant activation over MRI
            act_cmap = hot(256);
            aIdx = max(1, min(256, round(actNorm * 255) + 1));
            aRGB = ind2rgb(aIdx, act_cmap);
            aAlpha = double(actSigMask) .* 0.8;
            hAct2 = image(ax5, bm_x, bm_z, aRGB);
            hAct2.AlphaData = aAlpha;
            % ROI outline
            if any(roiResamp(:))
                [B_roi, ~] = bwboundaries(roiResamp);
                for k = 1:length(B_roi)
                    bnd = B_roi{k};
                    bnd_x = interp1(1:length(ud.x_out_mm), ud.x_out_mm, bnd(:,2), 'linear', 'extrap') + ud.bm_dx_mm;
                    bnd_z = interp1(1:length(ud.z_out_mm), ud.z_out_mm, bnd(:,1), 'linear', 'extrap') + ud.bm_dz_mm;
                    plot(ax5, bnd_x, bnd_z, 'g-', 'LineWidth', 1.5);
                end
            end
            overlayLabel = 'Activation';
        elseif overlayState == 2
            % Full correlation map over MRI (diverging blue-white-red)
            nC = 256;
            halfN = nC / 2;
            blue2white = [linspace(0,1,halfN)', linspace(0,1,halfN)', ones(halfN,1)];
            white2red  = [ones(halfN,1), linspace(1,0,halfN)', linspace(1,0,halfN)'];
            divCmap = [blue2white; white2red];
            cLim = max(abs(corrRaw(actMask)));
            if cLim == 0, cLim = 1; end
            corrScaled = (corrRaw / cLim + 1) / 2;
            cIdx = max(1, min(nC, round(corrScaled * (nC - 1)) + 1));
            cRGB = ind2rgb(cIdx, divCmap);
            cAlpha = double(actMask) .* 0.7;
            hCorr2 = image(ax5, bm_x, bm_z, cRGB);
            hCorr2.AlphaData = cAlpha;
            overlayLabel = sprintf('Correlation (r: %.2f to %.2f)', -cLim, cLim);
        else
            % Doppler over MRI
            fusi_cmap = feval(ud.fusiColormap, 256);
            fIdx = max(1, min(256, round(fusiNorm * 255) + 1));
            fRGB = ind2rgb(fIdx, fusi_cmap);
            fAlpha = double(fusiMask) .* ud.fusiAlpha .* (fusiNorm .^ 0.6);
            hF = image(ax5, bm_x, bm_z, fRGB);
            hF.AlphaData = fAlpha;
            overlayLabel = 'fUSI';
        end
 
        if ud.hasLGN
            lgn = ud.lgnSlices{ai}; lgn(isnan(lgn)) = 0;
            lgnMax = max(lgn(:));
            if lgnMax > 0
                lgnNorm = lgn / lgnMax;
                lgnIdx = max(1, min(256, round(lgnNorm * 255) + 1));
                hLGN2 = image(ax5, ud.x_out_mm, ud.z_out_mm, ind2rgb(lgnIdx, winter(256)));
                hLGN2.AlphaData = double(lgnNorm > 0.1) * 0.6;
            end
        end
        hold(ax5, 'off');
        axis(ax5, 'image'); whiteAxis(ax5);
        xlim(ax5, xLim);
        ylim(ax5, yLim);
        xlabel(ax5, 'Distance from midline (mm)'); ylabel(ax5, 'Depth (mm)');
        title(ax5, sprintf('%s / MRI (%.1f\\circ) | %s', overlayLabel, ang, adjStr), 'Color', 'w');
    end  % if ud.fUSI
 
    % Store ax4 handle for use by runFusiAnalysis
    if ud.fUSI, ud.ax4Handle = ax4; end
    fig.UserData = ud;
end
 
%% ===================== CYCLE OVERLAY CALLBACK =====================
function cycleOverlay(fig)
    ud = fig.UserData;
    ud.overlayState = mod(ud.overlayState + 1, 3);  % cycle 0 -> 1 -> 2 -> 0
    overlayNames = {'Doppler', 'Activation', 'Correlation'};
    fig.UserData = ud;
    fprintf('Overlay: %s (state=%d, sigPixels=%d)\n', ...
        overlayNames{ud.overlayState + 1}, ud.overlayState, sum(ud.actSigMask(:)));
    try
        updateDisplay(fig);
    catch ME
        fprintf(2, 'ERROR in updateDisplay: %s\n', ME.message);
        for si = 1:length(ME.stack)
            fprintf(2, '  %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
        end
    end
    figure(fig);
end

%% ===================== SAVE PANELS PDF (Save button) =====================
function savePanelsPDF(fig)
    ud = fig.UserData;
    ang = ud.angles(ud.currentIdx);

    baseDir = ud.saveDir;
    dateFolder = datestr(now, 'yyyy-mm-dd');
    outDir = fullfile(baseDir, dateFolder);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % --- Save the entire main figure as a single PDF + SVG ---
    fname_main = fullfile(outDir, sprintf('%s_%.1fdeg_main_panels.pdf', ud.subject, ang));
    exportPdfSvg(fig, fname_main, 'BackgroundColor', 'k');

    % --- Also export each tile row as a separate figure for cleaner output ---
    % Extract all tile axes from the figure
    allAxes = findobj(fig, 'Type', 'axes');
    % Sort by position (top to bottom, left to right)
    axPos = cell2mat(arrayfun(@(a) a.Position, allAxes, 'UniformOutput', false));
    [~, sortIdx] = sortrows(axPos, [-2, 1]);  % sort by Y descending, then X ascending
    allAxes = allAxes(sortIdx);

    panelNames = {};
    for k = 1:length(allAxes)
        ttl = get(get(allAxes(k), 'Title'), 'String');
        if isempty(ttl), ttl = sprintf('panel_%d', k); end
        % Clean title for filename
        ttl = regexprep(ttl, '[^a-zA-Z0-9_]', '_');
        ttl = regexprep(ttl, '_+', '_');
        panelNames{end+1} = ttl; %#ok<AGROW>

        fTmp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 500]);
        axCopy = copyobj(allAxes(k), fTmp);
        set(axCopy, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
        set(axCopy, 'XColor', 'k', 'YColor', 'k', 'Color', 'w');
        ttlObj = get(axCopy, 'Title');
        set(ttlObj, 'Color', 'k');

        panelFile = fullfile(outDir, sprintf('%s_%.1fdeg_%s.pdf', ud.subject, ang, ttl));
        exportPdfSvg(fTmp, panelFile);
        close(fTmp);
    end

    fprintf('=== All panels saved for %s at %.1f deg to %s ===\n\n', ud.subject, ang, outDir);

    if ~ud.fUSI
        fprintf('fUSI is off, so the correlation/activation figures and the analysis were skipped.\n');
        return
    end

    % --- Standalone MRI overlay figures ---
    ai = ud.currentIdx;
    mrNorm  = ud.mrNorms{ai};
    bm_x = ud.x_out_mm + ud.bm_dx_mm;
    bm_z = ud.z_out_mm + ud.bm_dz_mm;
    actNorm    = ud.actNorm;
    actSigMask = ud.actSigMask;
    corrRaw    = ud.corrRawGrid;
    actMask    = ud.actMask;

    % Shared axis limits
    margin = 3;
    mri_xR = [-ud.fovWidth_mm/2, ud.fovWidth_mm/2];
    mri_zR = [-ud.txMarkerThick_mm, ud.fovDepth_mm];
    bm_xR  = mri_xR + ud.bm_dx_mm;
    bm_zR  = mri_zR + ud.bm_dz_mm;
    xLim = [min(mri_xR(1), bm_xR(1)) - margin, max(mri_xR(2), bm_xR(2)) + margin];
    yLim = [min(mri_zR(1), bm_zR(1)) - margin, max(mri_zR(2), bm_zR(2)) + margin];

    % ---- 1) Correlation on MRI (no LGN) ----
    fTmp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 550]);
    ax = axes(fTmp);
    imagesc(ax, ud.x_out_mm, ud.z_out_mm, mrNorm);
    colormap(ax, gray); hold(ax, 'on');
    nC = 256; halfN = nC / 2;
    blue2white = [linspace(0,1,halfN)', linspace(0,1,halfN)', ones(halfN,1)];
    white2red  = [ones(halfN,1), linspace(1,0,halfN)', linspace(1,0,halfN)'];
    divCmap = [blue2white; white2red];
    cLim = max(abs(corrRaw(actMask)));
    if cLim == 0, cLim = 1; end
    corrScaled = (corrRaw / cLim + 1) / 2;
    cIdx = max(1, min(nC, round(corrScaled * (nC - 1)) + 1));
    cRGB = ind2rgb(cIdx, divCmap);
    cAlpha = double(actMask) .* 0.7;
    hC = image(ax, bm_x, bm_z, cRGB);
    hC.AlphaData = cAlpha;
    hold(ax, 'off'); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, sprintf('%s: Correlation on MRI (%.1f\\circ, r: %.2f to %.2f)', ...
        ud.subject, ang, -cLim, cLim));
    % Colorbar via hidden overlay axes (preserves gray colormap on main axes)
    axPos = get(ax, 'Position');
    axCb = axes(fTmp, 'Position', axPos, 'Visible', 'off');
    colormap(axCb, divCmap); clim(axCb, [-cLim cLim]);
    cb = colorbar(axCb, 'Location', 'eastoutside');
    ylabel(cb, 'Correlation (r)');
    % Resize main axes to match the colorbar-adjusted layout
    set(ax, 'Position', get(axCb, 'Position'));
    pdfName = fullfile(outDir, sprintf('%s_%.1fdeg_correlation_on_MRI.pdf', ud.subject, ang));
    exportPdfSvg(fTmp, pdfName);
    close(fTmp);

    % ---- 2) Activation on MRI (no LGN) ----
    fTmp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 550]);
    ax = axes(fTmp);
    imagesc(ax, ud.x_out_mm, ud.z_out_mm, mrNorm);
    colormap(ax, gray); hold(ax, 'on');
    act_cmap = hot(256);
    aIdx = max(1, min(256, round(actNorm * 255) + 1));
    aRGB = ind2rgb(aIdx, act_cmap);
    aAlpha = double(actSigMask) .* 0.8;
    hA = image(ax, bm_x, bm_z, aRGB);
    hA.AlphaData = aAlpha;
    hold(ax, 'off'); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, sprintf('%s: Activation on MRI (%.1f\\circ, z>%.2f)', ...
        ud.subject, ang, ud.stim.zscore));
    pdfName = fullfile(outDir, sprintf('%s_%.1fdeg_activation_on_MRI.pdf', ud.subject, ang));
    exportPdfSvg(fTmp, pdfName);
    close(fTmp);

    % ---- 3) Activation on MRI (with LGN) ----
    fTmp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 550]);
    ax = axes(fTmp);
    imagesc(ax, ud.x_out_mm, ud.z_out_mm, mrNorm);
    colormap(ax, gray); hold(ax, 'on');
    aIdx = max(1, min(256, round(actNorm * 255) + 1));
    aRGB = ind2rgb(aIdx, act_cmap);
    aAlpha = double(actSigMask) .* 0.8;
    hA = image(ax, bm_x, bm_z, aRGB);
    hA.AlphaData = aAlpha;
    if ud.hasLGN
        lgn = ud.lgnSlices{ai}; lgn(isnan(lgn)) = 0;
        lgnMax = max(lgn(:));
        if lgnMax > 0
            lgnNorm = lgn / lgnMax;
            lgnIdx = max(1, min(256, round(lgnNorm * 255) + 1));
            hLGN = image(ax, ud.x_out_mm, ud.z_out_mm, ind2rgb(lgnIdx, winter(256)));
            hLGN.AlphaData = double(lgnNorm > 0.1) * 0.6;
        end
    end
    hold(ax, 'off'); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, sprintf('%s: Activation + LGN on MRI (%.1f\\circ, z>%.2f)', ...
        ud.subject, ang, ud.stim.zscore));
    pdfName = fullfile(outDir, sprintf('%s_%.1fdeg_activation_LGN_on_MRI.pdf', ud.subject, ang));
    exportPdfSvg(fTmp, pdfName);
    close(fTmp);

    % --- Run fUSI analysis (equivalent to pressing P) and save time course + epoch ---
    fprintf('  Running fUSI analysis for PDF export...\n');
    try
        runFusiAnalysis(fig, false);
        % Find the analysis figure by name
        anFigName = sprintf('fUSI Analysis at %.1f deg', ang);
        anFig = findobj('Type', 'figure', 'Name', anFigName);
        if ~isempty(anFig)
            % Extract subplots by Tag into separate figures
            anAxes = findobj(anFig(1), 'Type', 'axes');
            for ia = 1:length(anAxes)
                axTag = get(anAxes(ia), 'Tag');

                saveThis = false;
                if strcmp(axTag, 'time_course')
                    pdfName = sprintf('%s_%.1fdeg_ROI_time_course.pdf', ud.subject, ang);
                    figW = 900; figH = 350;
                    saveThis = true;
                elseif strcmp(axTag, 'epoch_average_baseline_ref')
                    pdfName = sprintf('%s_%.1fdeg_epoch_average_baseline_ref.pdf', ud.subject, ang);
                    figW = 700; figH = 400;
                    saveThis = true;
                elseif strcmp(axTag, 'epoch_average')
                    pdfName = sprintf('%s_%.1fdeg_epoch_average.pdf', ud.subject, ang);
                    figW = 700; figH = 400;
                    saveThis = true;
                end

                if saveThis
                    fTmp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 figW figH]);
                    axCopy = copyobj(anAxes(ia), fTmp);
                    set(axCopy, 'Units', 'normalized', 'Position', [0.12 0.15 0.82 0.75]);
                    % Re-create legend by finding copied objects via DisplayName
                    lgdObj = get(anAxes(ia), 'Legend');
                    if ~isempty(lgdObj) && isvalid(lgdObj)
                        lgdLoc = lgdObj.Location;
                        lgdHandles = gobjects(0);
                        for si = 1:length(lgdObj.String)
                            h = findobj(axCopy, 'DisplayName', lgdObj.String{si});
                            if ~isempty(h)
                                lgdHandles(end+1) = h(1); %#ok<AGROW>
                            end
                        end
                        if ~isempty(lgdHandles)
                            legend(axCopy, lgdHandles, 'Location', lgdLoc);
                        end
                    end
                    outFile = fullfile(outDir, pdfName);
                    exportPdfSvg(fTmp, outFile);
                    close(fTmp);
                end
            end
            close(anFig(1));
        else
            fprintf('  Analysis figure not found, skipping.\n');
        end

        % Save bar plot if it exists
        barFigName = sprintf('Mean CBV: %s', ud.subjectLabel);
        barFig = findobj('Type', 'figure', 'Name', barFigName);
        if ~isempty(barFig)
            barPdfName = fullfile(outDir, sprintf('%s_%.1fdeg_mean_CBV_bar.pdf', ud.subject, ang));
            exportPdfSvg(barFig(1), barPdfName);
            close(barFig(1));
        end
    catch ME
        fprintf(2, '  Could not generate analysis figures: %s\n', ME.message);
    end

    fprintf('=== Save complete for %s at %.1f deg ===\n\n', ud.subject, ang, outDir);
end

%% ===================== EXPORT CROSS-ANIMAL EPOCH PDF =====================
function exportCrossAnimalPDF(fig)
    ud = fig.UserData;

    baseDir = ud.saveDir;
    dateFolder = datestr(now, 'yyyy-mm-dd');
    outDir = fullfile(baseDir, dateFolder);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fprintf('\n=== Exporting cross-animal epoch average ===\n');

    % Generate the cross-animal figure
    generateCrossAnimalEpochPlot(fig);

    % Find it by name and export
    crossFig = findobj('Type', 'figure', 'Name', 'Cross-Animal Epoch Average (Lights On vs Off)');
    if ~isempty(crossFig)
        fname = fullfile(outDir, 'epoch_average_across_animals.pdf');
        exportPdfSvg(crossFig(1), fname);
        close(crossFig(1));
    else
        fprintf('  Cross-animal figure not found.\n');
    end

    % Find and export combined bar plot
    nAnimals = length(ud.animalStore);
    barFigName = sprintf('Mean CBV Across Animals (N=%d)', nAnimals);
    barFig = findobj('Type', 'figure', 'Name', barFigName);
    if ~isempty(barFig)
        fname = fullfile(outDir, 'mean_CBV_bar_across_animals.pdf');
        exportPdfSvg(barFig(1), fname);
        close(barFig(1));
    end

    % Generate and export composite figure (all subjects + cross-animal)
    compFig = generateCompositeFigure(fig);
    if ~isempty(compFig)
        fname = fullfile(outDir, 'composite_epoch_summary.pdf');
        exportPdfSvg(compFig, fname);
        close(compFig);
    end

    fprintf('=== Cross-animal export complete ===\n\n');
end
 
%% ===================== COMPOSITE EPOCH FIGURE =====================
function compFig = generateCompositeFigure(fig)
    % Creates a single composite figure with per-subject epoch time courses
    % and bar plots in the top rows, and cross-animal pooled results on the
    % bottom row.  Layout: 3 rows x 2 columns (time course | bar plot).
    %
    % Row ordering matches the uploaded reference image:
    %   Row 1: Subject H  (Hobbes, animalStore index 2)
    %   Row 2: Subject C  (Calvin, animalStore index 1)
    %   Row 3: Cross-animal average

    ud = fig.UserData;
    nAnimals = length(ud.animalStore);

    % Verify all animals have epoch data
    for si = 1:nAnimals
        if isempty(ud.animalStore{si}.epochData)
            fprintf(2, '  Missing epoch data for %s. Cannot build composite.\n', ...
                ud.animalStore{si}.subject);
            compFig = [];
            return;
        end
    end

    % Display order: Hobbes (index 2) first, Calvin (index 1) second
    displayOrder = [2, 1];

    % Pool data for cross-animal row
    allDataStim = [];
    allDataCont = [];
    restPeriod  = [];
    for si = 1:nAnimals
        ed = ud.animalStore{si}.epochData;
        allDataStim = [allDataStim; ed.dataStim]; %#ok<AGROW>
        allDataCont = [allDataCont; ed.dataCont]; %#ok<AGROW>
        restPeriod  = ed.restPeriod;
    end
    nPooled = size(allDataStim, 1);

    pooledStimMean = mean(allDataStim, 1);
    pooledStimSEM  = std(allDataStim, 0, 1) / sqrt(nPooled);
    pooledContMean = mean(allDataCont, 1);
    pooledContSEM  = std(allDataCont, 0, 1) / sqrt(nPooled);
    % Prepend crossover point (t=0 of each condition = t=end of the other)
    stimEndVal = pooledStimMean(end);  contEndVal = pooledContMean(end);
    stimEndSEM = pooledStimSEM(end);   contEndSEM = pooledContSEM(end);
    pooledStimMean = [contEndVal, pooledStimMean];
    pooledStimSEM  = [contEndSEM, pooledStimSEM];
    pooledContMean = [stimEndVal, pooledContMean];
    pooledContSEM  = [stimEndSEM, pooledContSEM];

    pooledStimEpochMean = mean(allDataStim, 2);
    pooledContEpochMean = mean(allDataCont, 2);

    % --- Create figure ---
    compFig = figure('Name', 'Composite Epoch Summary', ...
        'NumberTitle', 'off', 'Color', 'w', ...
        'Position', [100 50 1200 900]);
    tl = tiledlayout(compFig, 3, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    % ---- Per-subject rows ----
    for ri = 1:nAnimals
        si = displayOrder(ri);
        ed = ud.animalStore{si}.epochData;
        subLabel = ud.animalStore{si}.subjectLabel;
        titlePrefix = ud.animalStore{si}.epochTitlePrefix;

        % -- Time course (left) --
        ax = nexttile(tl);
        hold(ax, 'on');
        t_stim = 0:ed.restPeriod;
        t_cont = ed.restPeriod:(2*ed.restPeriod);
        fill(ax, [t_stim fliplr(t_stim)], ...
            100*[ed.stimMean+ed.stimSEM fliplr(ed.stimMean-ed.stimSEM)], ...
            'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hS = plot(ax, t_stim, 100*ed.stimMean, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Lights On');
        fill(ax, [t_cont fliplr(t_cont)], ...
            100*[ed.contMean+ed.contSEM fliplr(ed.contMean-ed.contSEM)], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hC = plot(ax, t_cont, 100*ed.contMean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Lights Off');
        hold(ax, 'off');
        xlim(ax, [0, 2*ed.restPeriod]);
        ylim(ax, [-5 5]);
        xlabel(ax, 'Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel(ax, 'Change in LGN Blood Volume (%)', 'FontSize', 14, 'FontWeight', 'bold');
        title(ax, sprintf('%s: %s', titlePrefix, subLabel), ...
            'FontSize', 16, 'FontWeight', 'bold');
        legend(ax, [hS, hC], {'Lights On', 'Lights Off'}, ...
            'Location', 'northeast', 'FontSize', 11, 'FontWeight', 'bold');
        set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);

        % -- Bar plot (right) --
        ax = nexttile(tl);
        stimEpochMeans = mean(ed.dataStim, 2);
        contEpochMeans = mean(ed.dataCont, 2);
        plotEpochBar(ax, stimEpochMeans, contEpochMeans, ...
            sprintf('Mean CBV: %s', subLabel));
    end

    % ---- Cross-animal row ----
    % Time course (left)
    ax = nexttile(tl);
    hold(ax, 'on');
    t_stim = 0:restPeriod;
    t_cont = restPeriod:(2*restPeriod);
    fill(ax, [t_stim fliplr(t_stim)], ...
        100*[pooledStimMean+pooledStimSEM fliplr(pooledStimMean-pooledStimSEM)], ...
        'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hS = plot(ax, t_stim, 100*pooledStimMean, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Lights On');
    fill(ax, [t_cont fliplr(t_cont)], ...
        100*[pooledContMean+pooledContSEM fliplr(pooledContMean-pooledContSEM)], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hC = plot(ax, t_cont, 100*pooledContMean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Lights Off');
    hold(ax, 'off');
    xlim(ax, [0 2*restPeriod]);
    ylim(ax, [-5 5]);
    xlabel(ax, 'Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel(ax, 'Change in LGN Blood Volume (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title(ax, sprintf('Epoch Average Across Animals (N=%d)', nAnimals), ...
        'FontSize', 16, 'FontWeight', 'bold');
    legend(ax, [hS, hC], {'Lights On', 'Lights Off'}, ...
        'Location', 'northeast', 'FontSize', 11, 'FontWeight', 'bold');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);

    % Bar plot (right)
    ax = nexttile(tl);
    plotEpochBar(ax, pooledStimEpochMean, pooledContEpochMean, ...
        sprintf('Mean CBV Across Animals (N=%d)', nAnimals));

    fprintf('  Composite figure generated.\n');

    % --- Summary statistics for all rows ---
    nComparisons = nAnimals + 1;  % per-subject tests + pooled
    fprintf('\n  === Composite Figure Summary (Bonferroni k=%d) ===\n', nComparisons);
    fprintf('  %-12s  %8s  %10s  %10s  %10s  %8s\n', ...
        'Subject', 'r', 'p(uncorr)', 'p(Bonf)', 'Cohen''s d', 'N epochs');
    fprintf('  %s\n', repmat('-', 1, 68));
    for ri = 1:nAnimals
        si = displayOrder(ri);
        ed = ud.animalStore{si}.epochData;
        sStim = mean(ed.dataStim, 2);
        sCont = mean(ed.dataCont, 2);
        sdPool = sqrt((var(sStim) + var(sCont)) / 2);
        if sdPool > 0
            dVal = (mean(sStim) - mean(sCont)) / sdPool;
        else
            dVal = NaN;
        end
        [~, pSubj] = ttest(sCont, sStim);
        pBonf = min(pSubj * nComparisons, 1.0);
        fprintf('  %-12s  %8.3f  %10.4f  %10.4f  %10.3f  %8d\n', ...
            ud.animalStore{si}.subjectLabel, ed.rCorr, pSubj, pBonf, dVal, ed.numAct);
    end
    % Pooled row
    pooledSD_c = sqrt((var(pooledStimEpochMean) + var(pooledContEpochMean)) / 2);
    if pooledSD_c > 0
        dPooled = (mean(pooledStimEpochMean) - mean(pooledContEpochMean)) / pooledSD_c;
    else
        dPooled = NaN;
    end
    [~, pPooled] = ttest(pooledContEpochMean, pooledStimEpochMean);
    pBonfPooled = min(pPooled * nComparisons, 1.0);
    fprintf('  %-12s  %8s  %10.4f  %10.4f  %10.3f  %8d\n', ...
        sprintf('Pooled(N=%d)', nAnimals), '-', pPooled, pBonfPooled, dPooled, nPooled);
    fprintf('  %s\n', repmat('-', 1, 68));
end

%% ===================== INLINE BAR PLOT HELPER (for composite) =====================
function plotEpochBar(ax, stimEpochMeans, contEpochMeans, titleStr)
    % Draws a Lights On / Lights Off bar comparison into existing axes.
    stimGrandMean = mean(stimEpochMeans) * 100;
    contGrandMean = mean(contEpochMeans) * 100;
    nEpochs = length(stimEpochMeans);
    stimSEM = std(stimEpochMeans) / sqrt(nEpochs) * 100;
    contSEM = std(contEpochMeans) / sqrt(nEpochs) * 100;

    [~, pVal] = ttest(contEpochMeans, stimEpochMeans);

    hold(ax, 'on');
    bar(ax, 1, stimGrandMean, 0.6, 'FaceColor', [0.85 0.2 0.2], ...
        'EdgeColor', 'k', 'LineWidth', 1.5);
    bar(ax, 2, contGrandMean, 0.6, 'FaceColor', [0.2 0.4 0.85], ...
        'EdgeColor', 'k', 'LineWidth', 1.5);
    errorbar(ax, 1, stimGrandMean, stimSEM, 'k', 'LineWidth', 1.5, 'CapSize', 12);
    errorbar(ax, 2, contGrandMean, contSEM, 'k', 'LineWidth', 1.5, 'CapSize', 12);

    % Significance bracket
    allVals = [stimGrandMean + stimSEM, stimGrandMean - stimSEM, ...
               contGrandMean + contSEM, contGrandMean - contSEM];
    yMax = max(allVals);
    yRange = max(allVals) - min(allVals);
    if yRange == 0, yRange = 1; end
    yBar  = yMax + 0.12 * yRange;
    yTick = yMax + 0.08 * yRange;
    yText = yBar + 0.02 * yRange;

    plot(ax, [1 2], [yBar yBar], 'k-', 'LineWidth', 1.5);
    plot(ax, [1 1], [yTick yBar], 'k-', 'LineWidth', 1.5);
    plot(ax, [2 2], [yTick yBar], 'k-', 'LineWidth', 1.5);

    if pVal < 0.001,     starStr = '***';
    elseif pVal < 0.01,  starStr = '**';
    elseif pVal < 0.05,  starStr = '*';
    else,                starStr = 'n.s.';
    end
    text(ax, 1.5, yText, starStr, ...
        'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');

    hold(ax, 'off');
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Lights On', 'Lights Off'});
    ylabel(ax, 'Change in LGN Blood Volume (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title(ax, titleStr, 'FontSize', 16, 'FontWeight', 'bold');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
    xlim(ax, [0.3 2.7]);
    ylim(ax, [-5 5]);
end


%% ===================== CONTROL PANEL =====================
function name = figureTitle(ud)
    name = sprintf('%s  |  Session %d of %d', ud.subject, ud.currentSession, ud.nSessions);
end

function buildControlPanel(fig)
    ud = fig.UserData;
    %otherName = ud.animalStore{3 - ud.currentAnimal}.subject;
    if numel(ud.animalStore) > 1
        otherIdx = 3 - ud.currentAnimal;
        otherName = ud.animalStore{otherIdx}.subject;
    else
        otherName = '';
    end

    bothAnalyzed = true;
    for si = 1:length(ud.animalStore)
        if isempty(ud.animalStore{si}.epochData)
            bothAnalyzed = false;
            break;
        end
    end

    overlayLabels = {'Show: Doppler', 'Show: Activation', 'Show: Correlation'};
    overlayColors = {[0.82 0.82 0.82], [0.93 0.74 0.74], [0.74 0.84 0.93]};
    stateIdx = ud.overlayState + 1;

    if ud.nSessions > 1
        sessEnable = 'on';
        combineEnable = 'on';
    else
        sessEnable = 'off';
        combineEnable = 'off';
    end
    if bothAnalyzed
        exportEnable = 'on';
    else
        exportEnable = 'off';
    end

    pnl = findobj(fig, 'Type', 'uipanel', 'Tag', 'ControlsPanel');
    isNewPanel = isempty(pnl);
    if isNewPanel
        pnl = uipanel(fig, ...
            'Units', 'normalized', ...
            'Position', [0.815 0.012 0.175 0.976], ...
            'BackgroundColor', [0.14 0.14 0.14], ...
            'ForegroundColor', [0.92 0.92 0.92], ...
            'HighlightColor', [0.40 0.40 0.40], ...
            'BorderType', 'line', ...
            'Title', 'Controls', ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'Tag', 'ControlsPanel');
    end

    neutral = [0.90 0.90 0.90];
    rows = {
        'Session', [], [], 'label', 'on', ''
        'Previous session', @(~,~) stepSession(fig, -1), neutral, 'button', sessEnable, 'Go to the previous session for this animal'
        'Next session', @(~,~) stepSession(fig, 1), neutral, 'button', sessEnable, 'Go to the next session for this animal'
        sprintf('Switch to %s', otherName), @(~,~) switchAnimal(fig), [0.76 0.86 0.96], 'button', 'on', 'Show the other animal'
        'Analysis', [], [], 'label', 'on', ''
        'Run analysis', @(~,~) safeCall(@() runFusiAnalysis(fig, false), 'runFusiAnalysis'), [0.98 0.84 0.60], 'button', 'on', 'Activation analysis. Uses the LGN, or a region you already drew.'
        'Draw ROI and analyze', @(~,~) safeCall(@() runFusiAnalysis(fig, true), 'runFusiAnalysis'), [0.96 0.70 0.60], 'button', 'on', 'Draw a region on the fUS image, then run activation analysis'
        'View', [], [], 'label', 'on', ''
        overlayLabels{stateIdx}, @(~,~) cycleOverlay(fig), overlayColors{stateIdx}, 'button', 'on', 'Cycle the fUS panels through Doppler, activation, and correlation'
        'Figures', [], [], 'label', 'on', ''
        'Activation on MRI', @(~,~) safeCall(@() generateStandaloneActivationMap(fig), 'generateStandaloneActivationMap'), neutral, 'button', 'on', 'Open the activation map on the MRI slice'
        'Combine sessions', @(~,~) safeCall(@() generateCombinedSessionFigure(fig), 'generateCombinedSessionFigure'), neutral, 'button', combineEnable, 'One figure with every session for this animal'
        'Compare animals', @(~,~) safeCall(@() generateCrossAnimalEpochPlot(fig), 'generateCrossAnimalEpochPlot'), neutral, 'button', 'on', 'Epoch average across both animals. Run analysis on each animal first.'
        'Save PDF', @(~,~) savePanelsPDF(fig), [0.62 0.86 0.62], 'button', 'on', 'Save the current panels as a PDF'
        'Export both animals', @(~,~) exportCrossAnimalPDF(fig), [0.62 0.74 0.94], 'button', exportEnable, 'Save the cross-animal comparison. Enabled after both animals are analyzed.'
        'Alignment', [], [], 'label', 'on', ''
        'Reset alignment', @(~,~) resetAlignment(fig), neutral, 'button', 'on', 'Move the ultrasound image back to its starting position'
        };

    if ~ud.fUSI
        % B-mode / MRI alignment only: drop everything that needs fUSI data
        fusiOnly = {'Analysis', 'Run analysis', 'Draw ROI and analyze', 'View', ...
                    overlayLabels{stateIdx}, 'Activation on MRI', 'Combine sessions', ...
                    'Compare animals', 'Export both animals'};
        rows = rows(~ismember(rows(:, 1), fusiOnly), :);
    end

    labelH = 0.026;
    btnH = 0.044;
    gap = 0.005;
    y = 0.965;
    btnCount = 0;
    for iRow = 1:size(rows, 1)
        if strcmp(rows{iRow, 4}, 'label')
            y = y - labelH;
            if isNewPanel
                uicontrol(pnl, 'Style', 'text', ...
                    'String', upper(rows{iRow, 1}), ...
                    'Units', 'normalized', ...
                    'Position', [0.06 y 0.88 labelH], ...
                    'BackgroundColor', pnl.BackgroundColor, ...
                    'ForegroundColor', [0.62 0.74 0.86], ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold');
            end
        else
            y = y - btnH;
            btnCount = btnCount + 1;
            btnTag = sprintf('ctlBtn%d', btnCount);
            if isNewPanel
                cbk = rows{iRow, 2};
                uicontrol(pnl, 'Style', 'pushbutton', ...
                    'String', rows{iRow, 1}, ...
                    'Units', 'normalized', ...
                    'Position', [0.06 y 0.88 btnH], ...
                    'BackgroundColor', rows{iRow, 3}, ...
                    'FontSize', 10, ...
                    'FontWeight', 'bold', ...
                    'Enable', rows{iRow, 5}, ...
                    'TooltipString', rows{iRow, 6}, ...
                    'Tag', btnTag, ...
                    'Callback', @(btn, ~) runAndDefocus(btn, cbk));
            else
                % Existing button: refresh its label, colour and enabled state only
                set(findobj(pnl, 'Tag', btnTag), ...
                    'String', rows{iRow, 1}, ...
                    'BackgroundColor', rows{iRow, 3}, ...
                    'Enable', rows{iRow, 5}, ...
                    'TooltipString', rows{iRow, 6});
            end
        end
        y = y - gap;
    end

    if ~isNewPanel, return; end

    uicontrol(pnl, 'Style', 'text', ...
        'String', {'Up / Down: beam angle', 'Left / Right: sagittal plane', 'W A S D: move image', 'R / E: rotate', 'Shift: smaller steps'}, ...
        'Units', 'normalized', ...
        'Position', [0.06 0.012 0.88 0.19], ...
        'BackgroundColor', pnl.BackgroundColor, ...
        'ForegroundColor', [0.78 0.78 0.78], ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9);
end

function stepSession(fig, delta)
    ud = fig.UserData;
    if ud.nSessions < 2
        fprintf('  Only 1 session for %s.\n', ud.subject);
        return
    end
    ud = saveAnimalState(ud);
    newSess = mod(ud.currentSession - 1 + delta, ud.nSessions) + 1;
    ud = loadSessionIntoUd(ud, newSess);
    fig.UserData = ud;
    fig.Name = figureTitle(ud);
    fprintf('\n=== Switched to %s session %d/%d ===\n', ud.subject, ud.currentSession, ud.nSessions);
    updateDisplay(fig);
end

function switchAnimal(fig)
    ud = saveAnimalState(fig.UserData);
    newIdx = 3 - ud.currentAnimal;
    ud = loadAnimalIntoUd(ud, newIdx);
    fig.UserData = ud;
    fig.Name = figureTitle(ud);
    fprintf('\n=== Switched to %s (session %d/%d) ===\n', ud.subject, ud.currentSession, ud.nSessions);
    updateDisplay(fig);
end

function resetAlignment(fig)
    ud = fig.UserData;
    ud.bm_dx_mm = 0;
    ud.bm_dz_mm = ud.txFaceOffset_mm;
    ud.bm_rot_deg = 0;
    idx = ud.currentAnimal;
    si  = ud.currentSession;
    ud.animalStore{idx}.sessions(si).bm_dx_init = ud.bm_dx_mm;
    ud.animalStore{idx}.sessions(si).bm_dz_init = ud.bm_dz_mm;
    ud.animalStore{idx}.sessions(si).bm_rot_init = ud.bm_rot_deg;
    ud.animalStore{idx}.sessions(si).savedCurrentIdx = ud.currentIdx;
    saveAlignmentToDisk(ud.animalStore{idx}.sessions(si));
    fig.UserData = ud;
    updateDisplay(fig);
end

% Runs a button's action, then takes keyboard focus off the button. A pushbutton keeps
% focus after a click, and the arrow keys / W A S D then stop reaching the figure.
function runAndDefocus(btn, cb)
    cb([], []);
    if isvalid(btn) && strcmp(btn.Enable, 'on')
        btn.Enable = 'off';
        drawnow;
        btn.Enable = 'on';
    end
end

function safeCall(fcn, errLabel)
    try
        fcn();
    catch ME
        fprintf(2, 'ERROR in %s: %s\n', errLabel, ME.message);
        for si = 1:length(ME.stack)
            fprintf(2, '  %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
        end
    end
end

%% ===================== KEY HANDLER =====================
function keyHandler(src, evt)
    ud = src.UserData;
    isShift = any(strcmp(evt.Modifier, 'shift'));
    tStep = 0.1 + 0.4 * ~isShift;
    rStep = 0.2 + 0.8 * ~isShift;
 
    if ~ud.fUSI && any(strcmp(evt.Key, {'p', 'm', 'f', 'g', 'c'}))
        fprintf('fUSI is off (fUSI = false in INPUTS), so ''%s'' does nothing.\n', evt.Key);
        return
    end

    redraw = true;
    saveAlign = true;   % false for keys that don't change the alignment
    switch evt.Key
        case 'uparrow'          % beam angle +
            if ud.currentIdx < length(ud.angles), ud.currentIdx = ud.currentIdx + 1;
            else, redraw = false; end
        case 'downarrow'        % beam angle -
            if ud.currentIdx > 1, ud.currentIdx = ud.currentIdx - 1;
            else, redraw = false; end

        case {'rightarrow', 'leftarrow'}   % step through sagittal planes
            saveAlign = false;
            if strcmp(evt.Key, 'rightarrow'), dirX = 1; else, dirX = -1; end
            sgn = sign(ud.sagUx_mm(end) - ud.sagUx_mm(1));   % +1 if x increases with plane index
            if sgn == 0, sgn = 1; end
            if isShift, nStep = 1; else, nStep = ud.sagStepVox; end
            newSag = min(max(ud.sagIdx + dirX * sgn * nStep, 1), numel(ud.sagUx_mm));
            if newSag == ud.sagIdx
                redraw = false;
            else
                ud.sagIdx = newSag;
            end
 
        case 'w', ud.bm_dz_mm = ud.bm_dz_mm - tStep;
        case 's', ud.bm_dz_mm = ud.bm_dz_mm + tStep;
        case 'a', ud.bm_dx_mm = ud.bm_dx_mm - tStep;
        case 'd', ud.bm_dx_mm = ud.bm_dx_mm + tStep;
        case 'r', ud.bm_rot_deg = ud.bm_rot_deg + rStep;
        case 'e', ud.bm_rot_deg = ud.bm_rot_deg - rStep;
 
        case 'z'
            resetAlignment(src);
            redraw = false;

        case 'p'
            safeCall(@() runFusiAnalysis(src, false), 'runFusiAnalysis');
            redraw = false;

        case 'm'
            safeCall(@() runFusiAnalysis(src, true), 'runFusiAnalysis');
            redraw = false;

        case 'f'
            safeCall(@() generateStandaloneActivationMap(src), 'generateStandaloneActivationMap');
            redraw = false;

        case 'g'
            safeCall(@() generateCrossAnimalEpochPlot(src), 'generateCrossAnimalEpochPlot');
            redraw = false;

        case 'n'
            stepSession(src, 1);
            redraw = false;

        case 'c'
            safeCall(@() generateCombinedSessionFigure(src), 'generateCombinedSessionFigure');
            redraw = false;

        case 't'
            switchAnimal(src);
            redraw = false;

        otherwise
            redraw = false;
    end
 
    if redraw
        % Update animalStore and persist to disk immediately
        idx = ud.currentAnimal;
        si  = ud.currentSession;
        ud.animalStore{idx}.sessions(si).bm_dx_init     = ud.bm_dx_mm;
        ud.animalStore{idx}.sessions(si).bm_dz_init     = ud.bm_dz_mm;
        ud.animalStore{idx}.sessions(si).bm_rot_init    = ud.bm_rot_deg;
        ud.animalStore{idx}.sessions(si).savedCurrentIdx = ud.currentIdx;
        ud.animalStore{idx}.sagIdx = ud.sagIdx;
        if saveAlign
            saveAlignmentToDisk(ud.animalStore{idx}.sessions(si));
        end

        src.UserData = ud;
        updateDisplay(src);
    end
end

%% ===================== fUSI ANALYSIS =====================
function runFusiAnalysis(fig, forceManual)
    ud  = fig.UserData;
    ai  = ud.currentIdx;
    ang = ud.angles(ai);
    st  = ud.stim;
 
    if forceManual
        fprintf('\n=== fUSI Activation Analysis at %.1f deg (Manual ROI, %d session(s)) ===\n', ang, ud.nSessions);
    else
        fprintf('\n=== fUSI Activation Analysis at %.1f deg (%d session(s)) ===\n', ang, ud.nSessions);
    end
 
    % --- Determine ROI ---
    % Reuse existing ROI if available (from a previous P/M press),
    % unless the user explicitly requested a new manual ROI.
    roiInNative = [];
    roiMaskDisplay = [];   % Full-resolution mask on output grid (manual ROI only)
    reusingROI = false;

    if ~forceManual && ~isempty(ud.roiResampled)
        reusingROI = true;
        if ud.roiIsManual
            fprintf('  Reusing existing manual ROI.\n');
        else
            fprintf('  Reusing existing auto ROI.\n');
        end
        % Convert stored ROI (on output grid) back to native fUSI coords
        [Xf, Zf] = meshgrid(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm);
        [X_OUT_local, Z_OUT_local] = meshgrid(ud.x_out_mm, ud.z_out_mm);
        if ud.roiIsManual && ud.bm_rot_deg ~= 0
            % Manual ROI is stored in the rotated display frame
            cx = (ud.x_out_mm(1) + ud.x_out_mm(end)) / 2;
            cz = (ud.z_out_mm(1) + ud.z_out_mm(end)) / 2;
            theta = deg2rad(ud.bm_rot_deg);
            Xf_rot =  cos(theta)*(Xf - cx) + sin(theta)*(Zf - cz) + cx;
            Zf_rot = -sin(theta)*(Xf - cx) + cos(theta)*(Zf - cz) + cz;
            roiOnGrid = interp2(X_OUT_local, Z_OUT_local, double(ud.roiResampled), ...
                Xf_rot, Zf_rot, 'nearest', 0);
        else
            roiOnGrid = interp2(X_OUT_local, Z_OUT_local, double(ud.roiResampled), ...
                Xf, Zf, 'nearest', 0);
        end
        roiInNative = roiOnGrid > 0.5;
    end

    if isempty(roiInNative) && ud.hasLGN && ~forceManual
        lgn = ud.lgnSlices{ai};
        lgn(isnan(lgn)) = 0;
        lgnMask = lgn > 0.1 * max(lgn(:));
 
        % Shift LGN mask into fUSI pixel space to check overlap
        % The fUSI is displayed shifted by (bm_dx, bm_dz) relative to MRI.
        % So to find LGN in fUSI space, shift LGN by (-bm_dx, -bm_dz) in pixels.
        dx_px = round(ud.bm_dx_mm / ud.pxSize_mm);
        dz_px = round(ud.bm_dz_mm / ud.pxSize_mm);
 
        lgnShifted = circshift(lgnMask, [-dz_px, -dx_px]);
        % Zero out wrapped edges
        if dz_px > 0, lgnShifted(end-dz_px+1:end, :) = false; end
        if dz_px < 0, lgnShifted(1:-dz_px, :) = false; end
        if dx_px > 0, lgnShifted(:, end-dx_px+1:end) = false; end
        if dx_px < 0, lgnShifted(:, 1:-dx_px) = false; end
 
        overlap = lgnShifted & ud.fusiMask;
        overlapFrac = sum(overlap(:)) / max(sum(lgnMask(:)), 1);
        fprintf('  LGN overlap with fUSI: %.0f%%\n', overlapFrac * 100);
 
        if overlapFrac > 0.1  % At least 10% overlap
            fprintf('  Using LGN as ROI automatically.\n');
 
            % Map the overlap mask (on output grid) back to native fUSI coords
            [Xf, Zf] = meshgrid(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm);
            [X_OUT_local, Z_OUT_local] = meshgrid(ud.x_out_mm, ud.z_out_mm);
            roiOnGrid = interp2(X_OUT_local, Z_OUT_local, double(overlap), Xf, Zf, 'nearest', 0);
            roiInNative = roiOnGrid > 0.5;
        else
            fprintf('  LGN does not sufficiently overlap with fUSI.\n');
        end
    end
 
    if isempty(roiInNative)
        % Use the stored fUSI panel axes handle
        ax4 = ud.ax4Handle;
        axes(ax4);
 
        % Let user zoom in before drawing ROI
        fprintf('  Zoom into the region of interest on the figure.\n');
        zoom on;
        input('  Press ENTER here when ready to draw ROI...', 's');
        zoom off;
 
        fprintf('  Now draw your ROI polygon...\n');
        ud = fig.UserData;
        ax4 = ud.ax4Handle;
        axes(ax4);
        roi = drawpolygon('Color', 'g', 'LineWidth', 1.5);
        wait(roi);
 
        % Create mask from drawn ROI on the base B-mode image
        axChildren = get(ax4, 'Children');
        axChildren = axChildren(isvalid(axChildren));
        isImg = arrayfun(@(h) strcmp(get(h, 'Type'), 'image'), axChildren);
        imgHandles = axChildren(isImg);
        roiMaskDisplay = createMask(roi, imgHandles(end));  % end = first plotted (base image)
 
        % Map from display grid to native fUSI coordinates
        % The mask was drawn on the rotated display. Instead of un-rotating
        % the mask (which loses edge pixels via crop), rotate the fUSI query
        % coordinates into the rotated display frame for the lookup.
        [Xf, Zf] = meshgrid(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm);
        [X_OUT_local, Z_OUT_local] = meshgrid(ud.x_out_mm, ud.z_out_mm);
        if ud.bm_rot_deg ~= 0
            % Rotation center matches imrotate: center of the coordinate grid
            cx = (ud.x_out_mm(1) + ud.x_out_mm(end)) / 2;
            cz = (ud.z_out_mm(1) + ud.z_out_mm(end)) / 2;
            theta = deg2rad(ud.bm_rot_deg);
            % Rotate fUSI coords into the rotated display frame
            Xf_rot =  cos(theta)*(Xf - cx) + sin(theta)*(Zf - cz) + cx;
            Zf_rot = -sin(theta)*(Xf - cx) + cos(theta)*(Zf - cz) + cz;
        else
            Xf_rot = Xf;
            Zf_rot = Zf;
        end
        roiOnGrid = interp2(X_OUT_local, Z_OUT_local, double(roiMaskDisplay), ...
            Xf_rot, Zf_rot, 'nearest', 0);
        roiInNative = roiOnGrid > 0.5;
 
        fprintf('  ROI drawn: %d pixels in native fUSI space.\n', sum(roiInNative(:)));
    end
 
    if sum(roiInNative(:)) < 5
        warning('ROI too small (%d pixels). Aborting analysis.', sum(roiInNative(:)));
        return;
    end
 
    % --- Use precomputed activation map ---
    ActMap    = ud.ActMap;
    threshold = ud.actThreshold;
    sigMask   = ActMap > threshold;
    stimPatternNorm = ud.stimPatternNorm;
 
    % Use session 1 for display time course
    DopAll  = ud.DopAll;
    nFrames = size(DopAll, 3);
 
    % Store ROI on output grid for overlay drawing (skip if reusing)
    if ~reusingROI
        % For manual ROI, use the original full-resolution mask directly to avoid
        % shape distortion from round-tripping through the coarse fUSI grid.
        if ~isempty(roiMaskDisplay)
            % Manual ROI is already in the rotated display frame; store as-is
            % and flag it so updateDisplay skips the rotation step.
            ud.roiResampled = roiMaskDisplay;
            ud.roiIsManual  = true;
        else
            [Xf, Zf] = meshgrid(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm);
            [X_OUT_local, Z_OUT_local] = meshgrid(ud.x_out_mm, ud.z_out_mm);
            roiResampledGrid = interp2(Xf, Zf, double(roiInNative), X_OUT_local, Z_OUT_local, 'nearest', 0) > 0.5;
            ud.roiResampled = roiResampledGrid;
            ud.roiIsManual  = false;
        end
        fig.UserData = ud;
        fprintf('  ROI stored. Use toggle button to switch overlay.\n');
    end
 
    % --- Extract ROI time course (session 1 for display) ---
    fprintf('  Extracting ROI time course (session 1 for display)...\n');
    roiDop = DopAll .* roiInNative;
    roiTimeCourse = squeeze(mean(mean(roiDop, 1), 2));
    grandMean = mean(roiTimeCourse);
    rBV = (roiTimeCourse - grandMean) ./ grandMean;
    rBV = movmean(rBV, 3);
 
    [R, ~] = corrcoef(rBV, stimPatternNorm(:));
    rCorr = R(1,2);

    % --- Multi-session epoch extraction and pooling ---
    expLen = st.restPeriod * 2;
    nSessions = ud.nSessions;
    
    % Collect epoch data from all sessions
    allDataStim = [];
    allDataCont = [];
    hasEpoch = true;
    
    for iSess = 1:nSessions
        DopSess = ud.sessions(iSess).DopAll;
        nFramesSess = size(DopSess, 3);
        totalAcq = st.numAct * expLen;
        
        if totalAcq > nFramesSess
            fprintf('  Session %d: not enough frames for epoch averaging (%d < %d). Skipping.\n', ...
                iSess, nFramesSess, totalAcq);
            continue;
        end
        
        % Extract ROI time course for this session
        roiDopSess = DopSess .* roiInNative;
        roiTC_sess = squeeze(mean(mean(roiDopSess, 1), 2));
        gmSess = mean(roiTC_sess);
        rBV_sess = (roiTC_sess - gmSess) ./ gmSess;
        rBV_sess = movmean(rBV_sess, 3);
        
        % Segment into epochs
        dataFirst_sess  = zeros(st.numAct, st.restPeriod);
        dataSecond_sess = zeros(st.numAct, st.restPeriod);
        for ep = 1:st.numAct
            start = (ep-1) * expLen + 1;
            dataFirst_sess(ep, :) = rBV_sess(start : start + st.restPeriod - 1);
            secondStart = start + st.restPeriod;
            if secondStart + st.restPeriod - 1 <= length(rBV_sess)
                dataSecond_sess(ep, :) = rBV_sess(secondStart : secondStart + st.restPeriod - 1);
            end
        end
        
        if st.stimFirst
            sessStim = dataFirst_sess;
            sessCont = dataSecond_sess;
        else
            sessCont = dataFirst_sess;
            sessStim = dataSecond_sess;
        end
        
        % Per-session correlation
        stimPN_sess = ud.sessions(iSess).stimPatternNorm;
        [R_sess, ~] = corrcoef(rBV_sess(1:min(end,length(stimPN_sess))), ...
            stimPN_sess(1:min(end,length(rBV_sess))));
        rCorr_sess = R_sess(1,2);
        
        fprintf('  Session %d: %d epochs, r=%.3f\n', iSess, st.numAct, rCorr_sess);
        
        allDataStim = [allDataStim; sessStim];  %#ok<AGROW>
        allDataCont = [allDataCont; sessCont];  %#ok<AGROW>
    end
    
    nTotalEpochs = size(allDataStim, 1);
    hasEpoch = nTotalEpochs > 0;
    pTtest = NaN;
    
    if hasEpoch
        dataStim = allDataStim;
        dataCont = allDataCont;
        
        stimMean = mean(dataStim, 1);
        stimSEM  = std(dataStim, 0, 1) / sqrt(nTotalEpochs);
        contMean = mean(dataCont, 1);
        contSEM  = std(dataCont, 0, 1) / sqrt(nTotalEpochs);
        t_epoch = 0:st.restPeriod;
        % t=0 of each condition = t=10 of the opposite condition
        stimEndVal = stimMean(end);  contEndVal = contMean(end);
        stimEndSEM = stimSEM(end);   contEndSEM = contSEM(end);
        stimMean = [contEndVal, stimMean];
        stimSEM  = [contEndSEM, stimSEM];
        contMean = [stimEndVal, contMean];
        contSEM  = [stimEndSEM, contSEM];

        [~, pTtest] = ttest(mean(dataCont, 2), mean(dataStim, 2));
        stimEpochMeans_tmp = mean(dataStim, 2);
        contEpochMeans_tmp = mean(dataCont, 2);
        pooledSD = sqrt((var(stimEpochMeans_tmp) + var(contEpochMeans_tmp)) / 2);
        if pooledSD > 0
            cohensD = (mean(stimEpochMeans_tmp) - mean(contEpochMeans_tmp)) / pooledSD;
        else
            cohensD = NaN;
        end
        
        if nSessions > 1
            fprintf('  Pooled %d epochs across %d sessions\n', nTotalEpochs, nSessions);
        end
        fprintf('  ROI correlation (sess 1): r=%.3f | Epoch t-test: p=%.4f | Cohen''s d=%.3f\n', rCorr, pTtest, cohensD);
    else
        fprintf('  ROI correlation: r=%.3f | Not enough frames for epoch averaging.\n', rCorr);
    end

    % --- Build shared time axis for fUSI frames ---
    if ud.hasIntan && ud.nAligned > 0
        fusiTimeAxis = ud.fusiFrameTimes(1:nFrames);
        timeLabel = 'Time (s)';
    else
        fusiTimeAxis = 1:nFrames;
        timeLabel = 'Frame';
    end

    % --- Plot results ---
    nRows = 3;  % default: activation maps, time course, epoch average
    if hasEpoch, nRows = nRows + 1; end  % baseline-referenced epoch plot
    if ud.hasIntan, nRows = nRows + 1; end
    if ud.hasEeg,   nRows = nRows + 1; end

    fAn = figure('Name', sprintf('fUSI Analysis at %.1f deg', ang), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1200 200*nRows]);
 
    % Activation map
    subplot(nRows, 2, 1);
    imagesc(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm, ActMap);
    axis image; colormap(gca, 'hot'); colorbar;
    title('Activation Map'); xlabel('mm'); ylabel('mm');
 
    % Significant activation overlaid on Doppler
    subplot(nRows, 2, 2);
    avgDop = mean(DopAll, 3);
    interpDopNorm = 1.5 * sqrt(avgDop) / max(sqrt(avgDop(:)));
    imagesc(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm, interpDopNorm(:,:,[1 1 1]));
    hold on;
    TransparencyMap = zeros(size(ActMap));
    TransparencyMap(sigMask) = 1;
    ih = imagesc(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm, ActMap);
    ih.AlphaData = TransparencyMap;
    colormap(gca, 'hot'); clim([0 1]);
    axis image; colorbar;
    [B_roi, ~] = bwboundaries(roiInNative);
    for k = 1:length(B_roi)
        bnd = B_roi{k};
        bnd_x = interp1(1:length(ud.fusi_xaxis_mm), ud.fusi_xaxis_mm, bnd(:,2), 'linear', 'extrap');
        bnd_z = interp1(1:length(ud.fusi_zaxis_mm), ud.fusi_zaxis_mm, bnd(:,1), 'linear', 'extrap');
        plot(bnd_x, bnd_z, 'g-', 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Significant Activation (z>%.2f)', st.zscore));
    xlabel('mm'); ylabel('mm');

    % Time course
    subplot(nRows, 2, [3 4]);
    set(gca, 'Tag', 'time_course');
    hTC = plot(fusiTimeAxis, rBV, 'Color', [0.635 0.078 0.184], 'LineWidth', 1.5, 'DisplayName', 'rBV'); hold on;
    stimScaled = stimPatternNorm / max(abs(stimPatternNorm)) * max(abs(rBV)) * 0.5;
    hSP = plot(fusiTimeAxis, stimScaled, 'k', 'LineWidth', 1, 'DisplayName', 'Stim Pattern');
    hold off;
    xlabel('Time (s)'); ylabel('rBV (A.U.)');
    title(sprintf('ROI Time Course | r=%.3f, t-test p=%.4f', rCorr, pTtest));
    legend([hTC, hSP], 'Location', 'best');
    xlim([fusiTimeAxis(1) fusiTimeAxis(end)]);

    % --- Intan EEG and Acceleration panels ---
    if ud.hasIntan && ud.nAligned > 0
        intanTime = ud.intanTimeAxis;
        tLimits   = [fusiTimeAxis(1), fusiTimeAxis(end)];
        nextRow = 2;

        trigRiseIntan = find(diff(ud.digTrigFull) == 1) + 1;
        trigFallIntan = find(diff(ud.digTrigFull) == -1);

        if ud.hasEeg
            nextRow = nextRow + 1;
            subplot(nRows, 2, [nextRow*2-1, nextRow*2]);
            eegPow = ud.eegCollect{1};
            eegWindowDur = 0.2;
            eegOverlapFrac = 0.5;
            eegTimeStep = eegWindowDur * (1 - eegOverlapFrac);
            nEegBins = size(eegPow, 2);
            eegTimeAxis = (0:nEegBins-1) * eegTimeStep;
            eegTimeAxis = eegTimeAxis + intanTime(1);

            for ch = 1:size(eegPow, 1)
                plot(eegTimeAxis, eegPow(ch,:), 'LineWidth', 1.2); hold on;
            end
            yl = ylim;
            for ti = 1:min(length(trigRiseIntan), length(trigFallIntan))
                tOn  = intanTime(trigRiseIntan(ti));
                tOff = intanTime(trigFallIntan(ti));
                patch([tOn tOff tOff tOn], [yl(1) yl(1) yl(2) yl(2)], ...
                    [0.8 0.8 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            end
            hold off;
            xlim(tLimits); xlabel(timeLabel); ylabel('Power (A.U.)');
            title('EEG Gamma Power (30-80 Hz)');
            if size(eegPow,1) == 2
                legend('Ch 1', 'Ch 2', 'Location', 'best');
            end
        end

        nextRow = nextRow + 1;
        subplot(nRows, 2, [nextRow*2-1, nextRow*2]);
        plot(intanTime, ud.accelFull, 'Color', [0.2 0.5 0.2], 'LineWidth', 0.8);
        hold on;
        accelRange = max(ud.accelFull) - min(ud.accelFull);
        trigScaled = double(ud.digTrigFull) * accelRange * 0.4 + min(ud.accelFull);
        plot(intanTime, trigScaled, 'Color', [0.1 0.1 0.8], 'LineWidth', 1);
        hold off;
        xlim(tLimits); xlabel(timeLabel); ylabel('Accel (V)');
        title('Accelerometer + Digital Trigger');
        legend('Accel', 'Digital Trigger', 'Location', 'best');
    end

    % --- Epoch average plot ---
    if hasEpoch
        epochRow = nRows - 1;  % baseline-referenced plot takes the last row
        subplot(nRows, 2, [epochRow*2-1, epochRow*2]);
        set(gca, 'Tag', 'epoch_average');
        hold on;
        t_stim = 0:st.restPeriod;
        t_cont = st.restPeriod:(2*st.restPeriod);
        fill([t_stim fliplr(t_stim)], 100*[stimMean+stimSEM fliplr(stimMean-stimSEM)], ...
            'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hStim = plot(t_stim, 100*stimMean, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Lights On');
        fill([t_cont fliplr(t_cont)], 100*[contMean+contSEM fliplr(contMean-contSEM)], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hCont = plot(t_cont, 100*contMean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Lights Off');
        hold off;
        xlim([0 2*st.restPeriod]);
        ylim([-5 5]);
        legend([hStim, hCont], 'Location', 'northeast', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Change in LGN Blood Volume (%)', 'FontSize', 18, 'FontWeight', 'bold');
        if nSessions > 1
            title(sprintf('%s: %s (%d sessions, %d epochs)', ud.epochTitlePrefix, ud.subjectLabel, nSessions, nTotalEpochs), ...
                'FontSize', 20, 'FontWeight', 'bold');
        else
            title(sprintf('%s: %s', ud.epochTitlePrefix, ud.subjectLabel), ...
                'FontSize', 20, 'FontWeight', 'bold');
        end
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
        fprintf('  Epoch t-test: p=%.4f\n', pTtest);

        % --- Baseline-referenced epoch plot (pooled across sessions) ---
        allDataStim_bl = [];
        allDataCont_bl = [];
        for iSess_bl = 1:nSessions
            DopSess_bl = ud.sessions(iSess_bl).DopAll;
            nFramesSess_bl = size(DopSess_bl, 3);
            if st.numAct * expLen > nFramesSess_bl, continue; end

            roiDopSess_bl = DopSess_bl .* roiInNative;
            roiTC_bl = squeeze(mean(mean(roiDopSess_bl, 1), 2));

            baselineFrames_bl = min(10, length(roiTC_bl));
            if ud.hasIntan && ud.nAligned > 0 && iSess_bl == 1
                t0 = ud.fusiFrameTimes(1);
                baselineIdx_bl = find(ud.fusiFrameTimes(1:nFramesSess_bl) <= t0 + 10);
                if ~isempty(baselineIdx_bl)
                    baselineFrames_bl = length(baselineIdx_bl);
                end
            end
            blMean = mean(roiTC_bl(1:baselineFrames_bl));
            rBV_bl_sess = (roiTC_bl - blMean) ./ blMean;
            rBV_bl_sess = movmean(rBV_bl_sess, 3);

            dataFirst_bl  = zeros(st.numAct, st.restPeriod);
            dataSecond_bl = zeros(st.numAct, st.restPeriod);
            for ep = 1:st.numAct
                start = (ep-1) * expLen + 1;
                dataFirst_bl(ep, :) = rBV_bl_sess(start : start + st.restPeriod - 1);
                secondStart = start + st.restPeriod;
                if secondStart + st.restPeriod - 1 <= length(rBV_bl_sess)
                    dataSecond_bl(ep, :) = rBV_bl_sess(secondStart : secondStart + st.restPeriod - 1);
                end
            end
            if st.stimFirst
                allDataStim_bl = [allDataStim_bl; dataFirst_bl];   %#ok<AGROW>
                allDataCont_bl = [allDataCont_bl; dataSecond_bl];  %#ok<AGROW>
            else
                allDataCont_bl = [allDataCont_bl; dataFirst_bl];   %#ok<AGROW>
                allDataStim_bl = [allDataStim_bl; dataSecond_bl];  %#ok<AGROW>
            end
        end
        dataStim_bl = allDataStim_bl;
        dataCont_bl = allDataCont_bl;
        nTotalEpochs_bl = size(dataStim_bl, 1);
        stimMean_bl = mean(dataStim_bl, 1);
        stimSEM_bl  = std(dataStim_bl, 0, 1) / sqrt(nTotalEpochs_bl);
        contMean_bl = mean(dataCont_bl, 1);
        contSEM_bl  = std(dataCont_bl, 0, 1) / sqrt(nTotalEpochs_bl);
        % t=0 of each condition = t=10 of the opposite condition
        stimEndVal_bl = stimMean_bl(end);  contEndVal_bl = contMean_bl(end);
        stimEndSEM_bl = stimSEM_bl(end);   contEndSEM_bl = contSEM_bl(end);
        stimMean_bl = [contEndVal_bl, stimMean_bl];
        stimSEM_bl  = [contEndSEM_bl, stimSEM_bl];
        contMean_bl = [stimEndVal_bl, contMean_bl];
        contSEM_bl  = [stimEndSEM_bl, contSEM_bl];

        [~, pTtest_bl] = ttest(mean(dataCont_bl, 2), mean(dataStim_bl, 2));

        blRow = epochRow + 1;
        subplot(nRows, 2, [blRow*2-1, blRow*2]);
        set(gca, 'Tag', 'epoch_average_baseline_ref');
        hold on;
        t_stim_bl = 0:st.restPeriod;
        t_cont_bl = st.restPeriod:(2*st.restPeriod);
        fill([t_stim_bl fliplr(t_stim_bl)], 100*[stimMean_bl+stimSEM_bl fliplr(stimMean_bl-stimSEM_bl)], ...
            'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hStim_bl = plot(t_stim_bl, 100*stimMean_bl, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Lights On');
        fill([t_cont_bl fliplr(t_cont_bl)], 100*[contMean_bl+contSEM_bl fliplr(contMean_bl-contSEM_bl)], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hCont_bl = plot(t_cont_bl, 100*contMean_bl, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Lights Off');
        hold off;
        xlim([0 2*st.restPeriod]);
        ylim([-5 5]);
        legend([hStim_bl, hCont_bl], 'Location', 'northeast', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Change in LGN Blood Volume (%)', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('%s: %s (rel. to first 10s)', ud.epochTitlePrefix, ud.subjectLabel), ...
            'FontSize', 20, 'FontWeight', 'bold');
        set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
        fprintf('  Baseline-ref epoch t-test: p=%.4f\n', pTtest_bl);

        % --- Store epoch data for cross-animal averaging (key G) ---
        epochData = struct();
        epochData.dataStim    = dataStim;       % all sessions pooled
        epochData.dataCont    = dataCont;        % all sessions pooled
        epochData.stimMean    = stimMean;
        epochData.stimSEM     = stimSEM;
        epochData.contMean    = contMean;
        epochData.contSEM     = contSEM;
        epochData.restPeriod  = st.restPeriod;
        epochData.numAct      = nTotalEpochs;   % total epochs across sessions
        epochData.numActPerSession = st.numAct;
        epochData.nSessions   = nSessions;
        epochData.subject     = ud.subject;
        epochData.angle       = ang;
        epochData.pTtest      = pTtest;
        epochData.rCorr       = rCorr;
        epochData.pCorr       = pTtest;  % use t-test p as primary p-value
        ud.animalStore{ud.currentAnimal}.epochData = epochData;
        ud.animalStore{ud.currentAnimal}.sessions = ud.sessions;
        fig.UserData = ud;
        if nSessions > 1
            fprintf('  Epoch data stored for %s (%d sessions, %d total epochs).\n', ...
                ud.subject, nSessions, nTotalEpochs);
        else
            fprintf('  Epoch data stored for %s. Use Compare animals after analyzing both.\n', ud.subject);
        end

        % --- Generate bar plot for this animal ---
        barTitle = sprintf('Mean CBV: %s', ud.subjectLabel);
        generateEpochBarPlot(mean(dataStim, 2), mean(dataCont, 2), barTitle);
    else
        fprintf('  Not enough frames for epoch averaging.\n');
    end
 
    fprintf('=== Analysis complete ===\n\n');
 
    % Refresh main figure to enable toggle button
    updateDisplay(fig);
end

%% ===================== STANDALONE FIGURES + PDF EXPORT (F key) =====================
function generateStandaloneActivationMap(fig)
    ud  = fig.UserData;
    ai  = ud.currentIdx;
    ang = ud.angles(ai);
    mrNorm  = ud.mrNorms{ai};
    mrSlice = ud.mrSlices{ai};

    fprintf('\n=== Generating standalone figures (%s, %.1f deg) ===\n', ud.subject, ang);

    % --- Create output folder ---
    baseDir = ud.saveDir;
    dateFolder = datestr(now, 'yyyy-mm-dd');
    outDir = fullfile(baseDir, dateFolder);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % --- Shared geometry ---
    if ud.bm_rot_deg ~= 0
        bmNorm   = imrotate(ud.bmNorm, ud.bm_rot_deg, 'bilinear', 'crop');
        mask     = imrotate(double(ud.sectorMask), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
        bmDisp   = imrotate(ud.bmResampled, ud.bm_rot_deg, 'bilinear', 'crop');
        bmDisp(isnan(bmDisp)) = 0;
        actNorm    = imrotate(ud.actNorm, ud.bm_rot_deg, 'bilinear', 'crop');
        actSigMask = imrotate(double(ud.actSigMask), ud.bm_rot_deg, 'nearest', 'crop') > 0.5;
    else
        bmNorm   = ud.bmNorm;
        mask     = ud.sectorMask;
        bmDisp   = ud.bmResampled; bmDisp(isnan(bmDisp)) = 0;
        actNorm    = ud.actNorm;
        actSigMask = ud.actSigMask;
    end

    bm_x = ud.x_out_mm + ud.bm_dx_mm;
    bm_z = ud.z_out_mm + ud.bm_dz_mm;

    margin = 3;
    mri_xRange = [-ud.fovWidth_mm/2, ud.fovWidth_mm/2];
    mri_zRange = [-ud.txMarkerThick_mm, ud.fovDepth_mm];
    bm_xRange  = [-ud.fovWidth_mm/2 + ud.bm_dx_mm, ud.fovWidth_mm/2 + ud.bm_dx_mm];
    bm_zRange  = [-ud.txMarkerThick_mm + ud.bm_dz_mm, ud.fovDepth_mm + ud.bm_dz_mm];
    xLim = [min(mri_xRange(1), bm_xRange(1)) - margin, max(mri_xRange(2), bm_xRange(2)) + margin];
    yLim = [min(mri_zRange(1), bm_zRange(1)) - margin, max(mri_zRange(2), bm_zRange(2)) + margin];

    figPos = [150 150 800 650];
    figs = gobjects(4, 1);  % track all figures for PDF export

    % ---- Figure 1: Significant Activation on MRI ----
    figs(1) = figure('Name', sprintf('[%s] Activation on MRI (%.1f deg)', ud.subject, ang), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', figPos);
    ax = axes(figs(1));
    imagesc(ax, ud.x_out_mm, ud.z_out_mm, mrNorm);
    colormap(ax, gray); hold(ax, 'on');
    act_cmap = hot(256);
    aIdx = max(1, min(256, round(actNorm * 255) + 1));
    aRGB = ind2rgb(aIdx, act_cmap);
    hAct = image(ax, bm_x, bm_z, aRGB);
    hAct.AlphaData = double(actSigMask) .* 0.8;
    hold(ax, 'off'); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, sprintf('%s: Significant Activation on MRI (%.1f\\circ, z>%.2f)', ...
        ud.subject, ang, ud.stim.zscore));
    fprintf('  1/4  Activation on MRI\n');

    % ---- Figure 2: MRI + LGN ----
    figs(2) = figure('Name', sprintf('[%s] MRI + LGN (%.1f deg)', ud.subject, ang), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', figPos + [30 -30 0 0]);
    ax = axes(figs(2));
    imagesc(ax, ud.x_out_mm, ud.z_out_mm, mrSlice);
    colormap(ax, gray); hold(ax, 'on');
    if ud.hasLGN
        lgn = ud.lgnSlices{ai}; lgn(isnan(lgn)) = 0;
        lgnMax = max(lgn(:));
        if lgnMax > 0
            lgnNorm = lgn / lgnMax;
            lgnIdx = max(1, min(256, round(lgnNorm * 255) + 1));
            hLGN = image(ax, ud.x_out_mm, ud.z_out_mm, ind2rgb(lgnIdx, winter(256)));
            hLGN.AlphaData = double(lgnNorm > 0.1);
        end
        ttl = sprintf('%s: MRI + LGN (%.1f\\circ)', ud.subject, ang);
    else
        ttl = sprintf('%s: MRI Slice (%.1f\\circ)', ud.subject, ang);
    end
    hold(ax, 'off'); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, ttl);
    fprintf('  2/4  MRI + LGN\n');

    % ---- Figure 3: Plain B-mode ----
    figs(3) = figure('Name', sprintf('[%s] B-mode (%.1f deg)', ud.subject, ang), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', figPos + [60 -60 0 0]);
    ax = axes(figs(3));
    imagesc(ax, bm_x, bm_z, bmDisp);
    colormap(ax, gray); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, sprintf('%s: B-mode (%.1f\\circ)', ud.subject, ang));
    fprintf('  3/4  B-mode\n');

    % ---- Figure 4: MRI + B-mode overlay + LGN ----
    figs(4) = figure('Name', sprintf('[%s] MRI + B-mode + LGN (%.1f deg)', ud.subject, ang), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', figPos + [90 -90 0 0]);
    ax = axes(figs(4));
    imagesc(ax, ud.x_out_mm, ud.z_out_mm, mrNorm);
    colormap(ax, gray); hold(ax, 'on');
    % B-mode overlay
    switch ud.overlayMode
        case 'fusion'
            usG = bmNorm .^ ud.fusionGamma;
            magRGB = cat(3, usG, zeros(size(usG)), usG);
            hUS = image(ax, bm_x, bm_z, magRGB);
            hUS.AlphaData = double(mask) .* min(usG * 1.5, 1);
        case 'checkerboard'
            [Nr, Nc] = size(bmNorm); [rr, cc] = ndgrid(1:Nr, 1:Nc);
            bs = ud.checkerBlockSize;
            checker = mod(floor((rr-1)/bs) + floor((cc-1)/bs), 2) == 1;
            hUS = image(ax, bm_x, bm_z, repmat(bmNorm, [1 1 3]));
            hUS.AlphaData = double(mask & checker);
        case 'alpha'
            cmap = feval(ud.overlayColormap, 256);
            bmIdx = max(1, min(256, round(bmNorm * 255) + 1));
            hUS = image(ax, bm_x, bm_z, ind2rgb(bmIdx, cmap));
            hUS.AlphaData = double(mask) .* ud.overlayAlpha .* (bmNorm .^ 0.8);
    end
    % LGN annotation on top
    if ud.hasLGN
        lgn = ud.lgnSlices{ai}; lgn(isnan(lgn)) = 0;
        lgnMax = max(lgn(:));
        if lgnMax > 0
            lgnNorm = lgn / lgnMax;
            lgnIdx = max(1, min(256, round(lgnNorm * 255) + 1));
            hLGN = image(ax, ud.x_out_mm, ud.z_out_mm, ind2rgb(lgnIdx, winter(256)));
            hLGN.AlphaData = double(lgnNorm > 0.1) * 0.6;
        end
    end
    hold(ax, 'off'); axis(ax, 'image');
    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Distance from midline (mm)'); ylabel(ax, 'Depth (mm)');
    title(ax, sprintf('%s: MRI + B-mode + LGN (%.1f\\circ)', ud.subject, ang));
    fprintf('  4/4  MRI + B-mode + LGN\n');

    % ---- Save all as PDFs + SVGs ----
    names = {'activation_on_MRI', 'MRI_LGN', 'Bmode', 'MRI_Bmode_LGN'};
    for fi = 1:4
        fname = fullfile(outDir, sprintf('%s_%.1fdeg_%s.pdf', ud.subject, ang, names{fi}));
        exportPdfSvg(figs(fi), fname);
    end

    % ---- Save cross-animal epoch average if available ----
    allAnalyzed = true;
    for si = 1:length(ud.animalStore)
        if isempty(ud.animalStore{si}.epochData)
            allAnalyzed = false;
            break;
        end
    end
    if allAnalyzed
        generateCrossAnimalEpochPlot(fig);
        % Find the figures that were just created
        epochFig = findobj('Type', 'figure', 'Name', 'Cross-Animal Epoch Average (Lights On vs Off)');
        if ~isempty(epochFig)
            fname = fullfile(outDir, 'epoch_average_across_animals.pdf');
            exportPdfSvg(epochFig(1), fname);
            close(epochFig(1));
        end
        nAnimals = length(ud.animalStore);
        barFigName = sprintf('Mean CBV Across Animals (N=%d)', nAnimals);
        barFig = findobj('Type', 'figure', 'Name', barFigName);
        if ~isempty(barFig)
            fname = fullfile(outDir, 'mean_CBV_bar_across_animals.pdf');
            exportPdfSvg(barFig(1), fname);
            close(barFig(1));
        end
    else
        fprintf('  Epoch average not saved (run P on both animals first).\n');
    end

    fprintf('=== All figures saved to %s ===\n\n', outDir);
end

%% ===================== COMBINED SESSION FIGURE (C key) =====================
function generateCombinedSessionFigure(fig)
    ud  = fig.UserData;

    % Save current session state so sessArr has the latest alignment values
    ud = saveAnimalState(ud);
    fig.UserData = ud;

    ai  = ud.currentIdx;
    ang = ud.angles(ai);
    st  = ud.stim;
    nSess = ud.nSessions;
    sessArr = ud.animalStore{ud.currentAnimal}.sessions;  % freshly saved

    if nSess < 2
        fprintf('Only 1 session loaded for %s. Nothing to combine.\n', ud.subject);
        return;
    end

    fprintf('\n=== Combined session figure for %s (%d sessions, %.1f deg) ===\n', ...
        ud.subject, nSess, ang);

    % --- Colormaps ---
    nC = 256; halfN = nC / 2;
    blue2white = [linspace(0,1,halfN)', linspace(0,1,halfN)', ones(halfN,1)];
    white2red  = [ones(halfN,1), linspace(1,0,halfN)', linspace(1,0,halfN)'];
    divCmap = [blue2white; white2red];
    actCmap = hot(256);
    whiteAxis = @(axH) set(axH, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');

    % ====================================================================
    %  Prepare Session 1 as the reference frame
    % ====================================================================
    ss1 = sessArr(1);
    if ss1.bm_rot_init ~= 0
        bmD_ref    = imrotate(ss1.bmResampled, ss1.bm_rot_init, 'bilinear', 'crop');
        bmD_ref(isnan(bmD_ref)) = 0;
        crr_ref    = imrotate(ss1.corrRawGrid, ss1.bm_rot_init, 'bilinear', 'crop');
        actMsk_ref = imrotate(double(ss1.actMaskGrid), ss1.bm_rot_init, 'nearest', 'crop') > 0.5;
        actN_ref   = imrotate(ss1.actNormGrid, ss1.bm_rot_init, 'bilinear', 'crop');
        sigM_ref   = imrotate(double(ss1.actSigResampled & ss1.actMaskGrid), ss1.bm_rot_init, 'nearest', 'crop') > 0.5;
    else
        bmD_ref    = ss1.bmResampled; bmD_ref(isnan(bmD_ref)) = 0;
        crr_ref    = ss1.corrRawGrid;
        actMsk_ref = ss1.actMaskGrid;
        actN_ref   = ss1.actNormGrid;
        sigM_ref   = ss1.actSigResampled & ss1.actMaskGrid;
    end

    bm_x_ref = ud.x_out_mm + ss1.bm_dx_init;
    bm_z_ref = ud.z_out_mm + ss1.bm_dz_init;
    margin = 3;
    bm_xR = [-ud.fovWidth_mm/2 + ss1.bm_dx_init, ud.fovWidth_mm/2 + ss1.bm_dx_init];
    bm_zR = [-ud.txMarkerThick_mm + ss1.bm_dz_init, ud.fovDepth_mm + ss1.bm_dz_init];
    xLim = [bm_xR(1) - margin, bm_xR(2) + margin];
    yLim = [bm_zR(1) - margin, bm_zR(2) + margin];

    % ====================================================================
    %  Co-register all sessions in UN-ROTATED space (translation only)
    % ====================================================================
    combCorr  = ss1.corrRawGrid;
    combAct   = ss1.actNormGrid;
    combSig   = ss1.actSigResampled & ss1.actMaskGrid;
    combMask  = ss1.actMaskGrid;
    combCount = double(ss1.actMaskGrid);

    fprintf('  Session 1 (reference): dx=%.1f dz=%.1f rot=%.1f\n', ...
        ss1.bm_dx_init, ss1.bm_dz_init, ss1.bm_rot_init);

    % Per-session panel data (each rotated by its own bm_rot)
    sessPanels = cell(nSess, 1);
    sessPanels{1} = struct('bmD', bmD_ref, 'crr', crr_ref, 'actMsk', actMsk_ref, ...
                           'actN', actN_ref, 'sigM', sigM_ref, ...
                           'bm_x', bm_x_ref, 'bm_z', bm_z_ref);

    for iSess = 2:nSess
        ss = sessArr(iSess);
        dx_rel = ss.bm_dx_init - ss1.bm_dx_init;
        dz_rel = ss.bm_dz_init - ss1.bm_dz_init;

        fprintf('  Session %d: dx=%.1f dz=%.1f rot=%.1f (rel: dx=%.1f dz=%.1f)\n', ...
            iSess, ss.bm_dx_init, ss.bm_dz_init, ss.bm_rot_init, dx_rel, dz_rel);

        % --- Per-session panel data (rotated by its own bm_rot) ---
        if ss.bm_rot_init ~= 0
            crr_disp     = imrotate(ss.corrRawGrid, ss.bm_rot_init, 'bilinear', 'crop');
            actMsk_disp  = imrotate(double(ss.actMaskGrid), ss.bm_rot_init, 'nearest', 'crop') > 0.5;
            bmD_disp     = imrotate(ss.bmResampled, ss.bm_rot_init, 'bilinear', 'crop');
            bmD_disp(isnan(bmD_disp)) = 0;
            actN_disp    = imrotate(ss.actNormGrid, ss.bm_rot_init, 'bilinear', 'crop');
            sigM_disp    = imrotate(double(ss.actSigResampled & ss.actMaskGrid), ss.bm_rot_init, 'nearest', 'crop') > 0.5;
        else
            crr_disp     = ss.corrRawGrid;
            actMsk_disp  = ss.actMaskGrid;
            bmD_disp     = ss.bmResampled; bmD_disp(isnan(bmD_disp)) = 0;
            actN_disp    = ss.actNormGrid;
            sigM_disp    = ss.actSigResampled & ss.actMaskGrid;
        end
        bm_x_s = ud.x_out_mm + ss.bm_dx_init;
        bm_z_s = ud.z_out_mm + ss.bm_dz_init;
        sessPanels{iSess} = struct('bmD', bmD_disp, 'crr', crr_disp, 'actMsk', actMsk_disp, ...
                                   'actN', actN_disp, 'sigM', sigM_disp, ...
                                   'bm_x', bm_x_s, 'bm_z', bm_z_s);

        % --- Combining: un-rotated, translation only ---
        crr_raw = ss.corrRawGrid;
        act_raw = ss.actNormGrid;
        sig_raw = double(ss.actSigResampled & ss.actMaskGrid);
        msk_raw = ss.actMaskGrid;

        [X_grid, Z_grid] = meshgrid(1:size(crr_raw, 2), 1:size(crr_raw, 1));
        dx_px = dx_rel / ud.pxSize_mm;
        dz_px = dz_rel / ud.pxSize_mm;

        crr_shifted = interp2(X_grid, Z_grid, crr_raw, X_grid - dx_px, Z_grid - dz_px, 'linear', 0);
        act_shifted = interp2(X_grid, Z_grid, act_raw, X_grid - dx_px, Z_grid - dz_px, 'linear', 0);
        sig_shifted = interp2(X_grid, Z_grid, sig_raw, X_grid - dx_px, Z_grid - dz_px, 'nearest', 0) > 0.5;
        msk_shifted = interp2(X_grid, Z_grid, double(msk_raw), X_grid - dx_px, Z_grid - dz_px, 'nearest', 0) > 0.5;

        combCorr(msk_shifted)  = combCorr(msk_shifted) + crr_shifted(msk_shifted);
        combAct(msk_shifted)   = combAct(msk_shifted) + act_shifted(msk_shifted);
        combCount(msk_shifted) = combCount(msk_shifted) + 1;
        combSig  = combSig | sig_shifted;
        combMask = combMask | msk_shifted;
    end
    combCount(combCount == 0) = 1;
    combCorr = combCorr ./ combCount;
    combAct  = combAct  ./ combCount;
    hiAct = max(combAct(:)); if hiAct > 0, combAct = combAct / hiAct; end

    % Apply Session 1's rotation once for display
    if ss1.bm_rot_init ~= 0
        combCorr = imrotate(combCorr, ss1.bm_rot_init, 'bilinear', 'crop');
        combAct  = imrotate(combAct, ss1.bm_rot_init, 'bilinear', 'crop');
        combSig  = imrotate(double(combSig), ss1.bm_rot_init, 'nearest', 'crop') > 0.5;
        combMask = imrotate(double(combMask), ss1.bm_rot_init, 'nearest', 'crop') > 0.5;
    end

    % ====================================================================
    %  Pre-compute RGB/Alpha for both views (all panels)
    % ====================================================================
    % Helper: build correlation RGB + alpha
    buildCorrOverlay = @(crr, msk) deal( ...
        ind2rgb(max(1, min(nC, round((crr / max(max(abs(crr(msk)),[], 'all'), 1e-6) + 1) / 2 * (nC-1)) + 1)), divCmap), ...
        double(msk) .* 0.7);

    % Helper: build activation RGB + alpha
    buildActOverlay = @(actN, sig) deal( ...
        ind2rgb(max(1, min(256, round(actN * 255) + 1)), actCmap), ...
        double(sig) .* 0.7);

    % Per-session overlay data
    corrOverlays = cell(nSess, 1);
    actOverlays  = cell(nSess, 1);
    corrTitles   = cell(nSess, 1);
    actTitles    = cell(nSess, 1);
    for iSess = 1:nSess
        sp = sessPanels{iSess};
        cLim = max(abs(sp.crr(sp.actMsk)));
        if isempty(cLim) || cLim == 0, cLim = 1; end
        corrScaled = (sp.crr / cLim + 1) / 2;
        cIdx = max(1, min(nC, round(corrScaled * (nC-1)) + 1));
        corrOverlays{iSess}.RGB   = ind2rgb(cIdx, divCmap);
        corrOverlays{iSess}.Alpha = double(sp.actMsk) .* 0.7;
        corrTitles{iSess} = sprintf('S%d Correlation (r: %.2f to %.2f)', iSess, -cLim, cLim);

        aIdx = max(1, min(256, round(sp.actN * 255) + 1));
        actOverlays{iSess}.RGB   = ind2rgb(aIdx, actCmap);
        actOverlays{iSess}.Alpha = double(sp.sigM) .* 0.7;
        actTitles{iSess} = sprintf('S%d Activation (%d sig px)', iSess, sum(sp.sigM(:)));
    end

    % Combined overlay data
    cLim_c = max(abs(combCorr(combMask)));
    if isempty(cLim_c) || cLim_c == 0, cLim_c = 1; end
    corrScaled_c = (combCorr / cLim_c + 1) / 2;
    cIdx_c = max(1, min(nC, round(corrScaled_c * (nC-1)) + 1));
    combCorrOverlay.RGB   = ind2rgb(cIdx_c, divCmap);
    combCorrOverlay.Alpha = double(combMask) .* 0.7;
    combCorrTitle = sprintf('Combined Corr (%d sess, r: %.2f to %.2f)', nSess, -cLim_c, cLim_c);

    aIdx_c = max(1, min(256, round(combAct * 255) + 1));
    combActOverlay.RGB   = ind2rgb(aIdx_c, actCmap);
    combActOverlay.Alpha = double(combSig) .* 0.7;
    combActTitle = sprintf('Combined Activation (%d sess, %d sig px)', nSess, sum(combSig(:)));

    % ====================================================================
    %  Create figure with toggle button
    % ====================================================================
    nCols = nSess + 1;
    combFig = figure('Name', sprintf('[%s] Combined Sessions (%.1f deg)', ud.subject, ang), ...
        'NumberTitle', 'off', 'Color', 'k', ...
        'Position', [40 120 min(520 * nCols, 2400) 550]);

    % Toggle button at top
    btnToggle = uicontrol(combFig, 'Style', 'pushbutton', ...
        'String', 'Switch to: Activation Map', ...
        'Units', 'normalized', 'Position', [0.35 0.94 0.30 0.05], ...
        'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [0.9 0.9 0.9], ...
        'Callback', @(src, ~) toggleCombinedView(src));

    tl = tiledlayout(combFig, 1, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl.InnerPosition = [0.04 0.02 0.92 0.88];  % leave room for button

    % --- Draw panels (start in correlation mode) ---
    overlayHandles = gobjects(nCols, 1);
    titleHandles   = gobjects(nCols, 1);

    for iSess = 1:nSess
        sp = sessPanels{iSess};
        ax = nexttile(tl, iSess);
        imagesc(ax, sp.bm_x, sp.bm_z, sp.bmD);
        colormap(ax, gray); hold(ax, 'on');
        hOv = image(ax, sp.bm_x, sp.bm_z, corrOverlays{iSess}.RGB);
        hOv.AlphaData = corrOverlays{iSess}.Alpha;
        hold(ax, 'off'); axis(ax, 'image');
        set(ax, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');
        xlim(ax, xLim); ylim(ax, yLim);
        xlabel(ax, 'mm'); ylabel(ax, 'Depth (mm)');
        titleHandles(iSess) = title(ax, corrTitles{iSess}, 'Color', 'w');
        overlayHandles(iSess) = hOv;
    end

    % Combined panel
    axComb = nexttile(tl, nCols);
    imagesc(axComb, bm_x_ref, bm_z_ref, bmD_ref);
    colormap(axComb, gray); hold(axComb, 'on');
    hOvC = image(axComb, bm_x_ref, bm_z_ref, combCorrOverlay.RGB);
    hOvC.AlphaData = combCorrOverlay.Alpha;
    hold(axComb, 'off'); axis(axComb, 'image');
    set(axComb, 'XColor', 'w', 'YColor', 'w', 'Color', 'k');
    xlim(axComb, xLim); ylim(axComb, yLim);
    xlabel(axComb, 'mm'); ylabel(axComb, 'Depth (mm)');
    titleHandles(nCols) = title(axComb, combCorrTitle, 'Color', [0.2 1 0.2], 'FontWeight', 'bold');
    overlayHandles(nCols) = hOvC;

    % Store toggle state in button UserData
    toggleData = struct();
    toggleData.currentView     = 'correlation';
    toggleData.overlayHandles  = overlayHandles;
    toggleData.titleHandles    = titleHandles;
    toggleData.corrOverlays    = {corrOverlays};
    toggleData.actOverlays     = {actOverlays};
    toggleData.corrTitles      = {corrTitles};
    toggleData.actTitles       = {actTitles};
    toggleData.combCorrOverlay = combCorrOverlay;
    toggleData.combActOverlay  = combActOverlay;
    toggleData.combCorrTitle   = combCorrTitle;
    toggleData.combActTitle    = combActTitle;
    toggleData.nSess           = nSess;
    toggleData.nCols           = nCols;
    btnToggle.UserData = toggleData;

    % ====================================================================
    %  ROI drawing on combined panel
    % ====================================================================
    fprintf('\n  Zoom into the combined panel (rightmost), then press ENTER to draw ROI.\n');
    axes(axComb);
    zoom on;
    input('  Press ENTER here when ready to draw ROI...', 's');
    zoom off;

    fprintf('  Draw your ROI polygon on the combined panel...\n');
    axes(axComb);
    roi = drawpolygon('Color', 'g', 'LineWidth', 1.5);
    wait(roi);

    axChildren = get(axComb, 'Children');
    axChildren = axChildren(isvalid(axChildren));
    isImg = arrayfun(@(h) strcmp(get(h, 'Type'), 'image'), axChildren);
    imgHandles = axChildren(isImg);
    roiMaskDisplay = createMask(roi, imgHandles(end));

    [Xf, Zf] = meshgrid(ud.fusi_xaxis_mm, ud.fusi_zaxis_mm);
    [X_OUT_local, Z_OUT_local] = meshgrid(ud.x_out_mm, ud.z_out_mm);
    if ss1.bm_rot_init ~= 0
        cx = (ud.x_out_mm(1) + ud.x_out_mm(end)) / 2;
        cz = (ud.z_out_mm(1) + ud.z_out_mm(end)) / 2;
        theta = deg2rad(ss1.bm_rot_init);
        Xf_rot =  cos(theta)*(Xf - cx) + sin(theta)*(Zf - cz) + cx;
        Zf_rot = -sin(theta)*(Xf - cx) + cos(theta)*(Zf - cz) + cz;
    else
        Xf_rot = Xf; Zf_rot = Zf;
    end
    roiOnGrid = interp2(X_OUT_local, Z_OUT_local, double(roiMaskDisplay), Xf_rot, Zf_rot, 'nearest', 0);
    roiInNative = roiOnGrid > 0.5;

    nROIpx = sum(roiInNative(:));
    fprintf('  ROI: %d pixels in native fUSI space.\n', nROIpx);
    if nROIpx < 5, fprintf('  ROI too small. Aborting.\n'); return; end

    % ====================================================================
    %  Multi-session epoch extraction
    % ====================================================================
    expLen = st.restPeriod * 2;
    allDataStim = []; allDataCont = [];
    for iSess = 1:nSess
        DopSess = sessArr(iSess).DopAll;
        nFramesSess = size(DopSess, 3);
        if st.numAct * expLen > nFramesSess
            fprintf('  Session %d: not enough frames. Skipping.\n', iSess); continue;
        end
        roiDopSess = DopSess .* roiInNative;
        roiTC = squeeze(mean(mean(roiDopSess, 1), 2));
        gm = mean(roiTC); rBV = (roiTC - gm) ./ gm; rBV = movmean(rBV, 3);
        stimPN = sessArr(iSess).stimPatternNorm;
        [R, ~] = corrcoef(rBV(1:min(end,length(stimPN))), stimPN(1:min(end,length(rBV))));
        dataFirst = zeros(st.numAct, st.restPeriod);
        dataSecond = zeros(st.numAct, st.restPeriod);
        for ep = 1:st.numAct
            start = (ep-1)*expLen+1;
            dataFirst(ep,:) = rBV(start:start+st.restPeriod-1);
            secondStart = start+st.restPeriod;
            if secondStart+st.restPeriod-1 <= length(rBV)
                dataSecond(ep,:) = rBV(secondStart:secondStart+st.restPeriod-1);
            end
        end
        if st.stimFirst, sessStim=dataFirst; sessCont=dataSecond;
        else, sessCont=dataFirst; sessStim=dataSecond; end
        fprintf('  Session %d: %d epochs, r=%.3f\n', iSess, st.numAct, R(1,2));
        allDataStim = [allDataStim; sessStim]; allDataCont = [allDataCont; sessCont]; %#ok<AGROW>
    end
    nTotalEpochs = size(allDataStim, 1);
    if nTotalEpochs == 0, fprintf('  No valid epochs.\n'); return; end

    dataStim = allDataStim; dataCont = allDataCont;
    stimMean = mean(dataStim,1); stimSEM = std(dataStim,0,1)/sqrt(nTotalEpochs);
    contMean = mean(dataCont,1); contSEM = std(dataCont,0,1)/sqrt(nTotalEpochs);
    t_epoch = 0:st.restPeriod;
    stimMean = [contMean(end), stimMean]; stimSEM = [contSEM(end), stimSEM];
    contMean = [stimMean(end), contMean]; contSEM = [stimSEM(end), contSEM];
    [~, pTtest] = ttest(mean(dataCont,2), mean(dataStim,2));
    stimEpochMeans = mean(dataStim,2); contEpochMeans = mean(dataCont,2);
    pooledSD = sqrt((var(stimEpochMeans)+var(contEpochMeans))/2);
    if pooledSD > 0, cohensD = (mean(stimEpochMeans)-mean(contEpochMeans))/pooledSD;
    else, cohensD = NaN; end
    fprintf('  Pooled %d epochs across %d sessions | p=%.4f | d=%.3f\n', nTotalEpochs, nSess, pTtest, cohensD);

    % Store epoch data for cross-animal averaging
    epochData = struct();
    epochData.dataStim = dataStim; epochData.dataCont = dataCont;
    epochData.stimMean = stimMean; epochData.stimSEM = stimSEM;
    epochData.contMean = contMean; epochData.contSEM = contSEM;
    epochData.restPeriod = st.restPeriod; epochData.numAct = nTotalEpochs;
    epochData.numActPerSession = st.numAct; epochData.nSessions = nSess;
    epochData.subject = ud.subject; epochData.angle = ang;
    epochData.pTtest = pTtest; epochData.rCorr = NaN; epochData.pCorr = pTtest;
    ud.animalStore{ud.currentAnimal}.epochData = epochData;
    fig.UserData = ud;
    fprintf('  Epoch data stored for %s. Use Compare animals after analyzing both.\n', ud.subject);

    % ====================================================================
    %  Epoch figure
    % ====================================================================
    stimColor = [0.85 0.33 0.10]; contColor = [0.00 0.45 0.74];
    epochFig = figure('Name', sprintf('[%s] Combined Epoch (%d sess, %.1f deg)', ud.subject, nSess, ang), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1100 450]);
    tl2 = tiledlayout(epochFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl2, 1); hold(ax, 'on');
    fill(ax, [t_epoch, fliplr(t_epoch)], [stimMean+stimSEM, fliplr(stimMean-stimSEM)], stimColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill(ax, [t_epoch, fliplr(t_epoch)], [contMean+contSEM, fliplr(contMean-contSEM)], contColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, t_epoch, stimMean, '-', 'Color', stimColor, 'LineWidth', 2.5);
    plot(ax, t_epoch, contMean, '-', 'Color', contColor, 'LineWidth', 2.5);
    legend(ax, {'', '', 'Lights ON', 'Lights OFF'}, 'Location', 'best', 'FontSize', 11);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'rCBV');
    title(ax, sprintf('%s: %s (%d sessions, %d epochs)', ud.epochTitlePrefix, ud.subjectLabel, nSess, nTotalEpochs), ...
        'FontSize', 14, 'FontWeight', 'bold');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
    hold(ax, 'off'); grid(ax, 'on');

    ax = nexttile(tl2, 2);
    barData = [mean(stimEpochMeans), mean(contEpochMeans)];
    barSEM_ = [std(stimEpochMeans)/sqrt(length(stimEpochMeans)), ...
               std(contEpochMeans)/sqrt(length(contEpochMeans))];
    b = bar(ax, barData, 0.6); hold(ax, 'on');
    b.FaceColor = 'flat'; b.CData = [stimColor; contColor];
    errorbar(ax, 1:2, barData, barSEM_, 'k.', 'LineWidth', 1.5, 'CapSize', 12);
    set(ax, 'XTickLabel', {'Lights ON', 'Lights OFF'}, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Mean rCBV');
    title(ax, sprintf('p = %.4f, d = %.3f', pTtest, cohensD), 'FontSize', 14, 'FontWeight', 'bold');
    hold(ax, 'off');

    % Save
    baseDir = ud.saveDir;
    dateFolder = datestr(now, 'yyyy-mm-dd');
    outDir = fullfile(baseDir, dateFolder);
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    exportPdfSvg(combFig, fullfile(outDir, sprintf('%s_%.1fdeg_combined_%dsess_corr.pdf', ud.subject, ang, nSess)));
    exportPdfSvg(epochFig, fullfile(outDir, sprintf('%s_%.1fdeg_combined_%dsess_epoch.pdf', ud.subject, ang, nSess)));
    fprintf('=== Combined session figures saved ===\n\n');
end

%% ===================== TOGGLE COMBINED VIEW =====================
function toggleCombinedView(btn)
    td = btn.UserData;
    oh = td.overlayHandles;
    th = td.titleHandles;
    nS = td.nSess;
    nC = td.nCols;
    corrOv = td.corrOverlays{1};
    actOv  = td.actOverlays{1};
    corrT  = td.corrTitles{1};
    actT   = td.actTitles{1};

    if strcmp(td.currentView, 'correlation')
        % Switch to activation
        for iSess = 1:nS
            if isvalid(oh(iSess))
                oh(iSess).CData     = actOv{iSess}.RGB;
                oh(iSess).AlphaData = actOv{iSess}.Alpha;
            end
            if isvalid(th(iSess)), th(iSess).String = actT{iSess}; end
        end
        if isvalid(oh(nC))
            oh(nC).CData     = td.combActOverlay.RGB;
            oh(nC).AlphaData = td.combActOverlay.Alpha;
        end
        if isvalid(th(nC))
            th(nC).String = td.combActTitle;
        end
        btn.String = 'Switch to: Correlation Map';
        td.currentView = 'activation';
    else
        % Switch to correlation
        for iSess = 1:nS
            if isvalid(oh(iSess))
                oh(iSess).CData     = corrOv{iSess}.RGB;
                oh(iSess).AlphaData = corrOv{iSess}.Alpha;
            end
            if isvalid(th(iSess)), th(iSess).String = corrT{iSess}; end
        end
        if isvalid(oh(nC))
            oh(nC).CData     = td.combCorrOverlay.RGB;
            oh(nC).AlphaData = td.combCorrOverlay.Alpha;
        end
        if isvalid(th(nC))
            th(nC).String = td.combCorrTitle;
        end
        btn.String = 'Switch to: Activation Map';
        td.currentView = 'correlation';
    end
    btn.UserData = td;
end

%% ===================== CROSS-ANIMAL EPOCH AVERAGE (G key) =====================
function generateCrossAnimalEpochPlot(fig)
    ud = fig.UserData;
    nAnimals = length(ud.animalStore);

    fprintf('\n=== Generating cross-animal epoch average ===\n');

    allDataStim = [];
    allDataCont = [];
    subjectMeans = struct();
    restPeriod = [];

    for si = 1:nAnimals
        ed = ud.animalStore{si}.epochData;
        if isempty(ed)
            fprintf(2, '  Missing epoch data for %s.\n', ud.animalStore{si}.subject);
            fprintf(2, '  Switch to %s, run analysis, then use Compare animals again.\n', ...
                ud.animalStore{si}.subject);
            return;
        end
        if isfield(ed, 'nSessions') && ed.nSessions > 1
            fprintf('  %s: %d epochs (%d sessions) at %.1f deg (r=%.3f, p=%.4f, t-test p=%.4f)\n', ...
                ed.subject, ed.numAct, ed.nSessions, ed.angle, ed.rCorr, ed.pCorr, ed.pTtest);
        else
            fprintf('  %s: %d epochs at %.1f deg (r=%.3f, p=%.4f, t-test p=%.4f)\n', ...
                ed.subject, ed.numAct, ed.angle, ed.rCorr, ed.pCorr, ed.pTtest);
        end

        allDataStim = [allDataStim; ed.dataStim];  %#ok<AGROW>
        allDataCont = [allDataCont; ed.dataCont];  %#ok<AGROW>
        subjectMeans(si).stimMean = ed.stimMean;
        subjectMeans(si).contMean = ed.contMean;
        subjectMeans(si).subject  = ed.subject;
        restPeriod = ed.restPeriod;
    end

    nPooled = size(allDataStim, 1);
    fprintf('  Pooled %d total epochs across %d animals\n', nPooled, nAnimals);

    pooledStimMean = mean(allDataStim, 1);
    pooledStimSEM  = std(allDataStim, 0, 1) / sqrt(nPooled);
    pooledContMean = mean(allDataCont, 1);
    pooledContSEM  = std(allDataCont, 0, 1) / sqrt(nPooled);
    t_epoch = 0:restPeriod;
    % t=0 of each condition = t=10 of the opposite condition
    stimEndVal_p = pooledStimMean(end);  contEndVal_p = pooledContMean(end);
    stimEndSEM_p = pooledStimSEM(end);   contEndSEM_p = pooledContSEM(end);
    pooledStimMean = [contEndVal_p, pooledStimMean];
    pooledStimSEM  = [contEndSEM_p, pooledStimSEM];
    pooledContMean = [stimEndVal_p, pooledContMean];
    pooledContSEM  = [stimEndSEM_p, pooledContSEM];

    pooledContEpochMean = mean(allDataCont, 2);
    pooledStimEpochMean = mean(allDataStim, 2);
    [~, pPooled] = ttest(pooledContEpochMean, pooledStimEpochMean);

    fCross = figure('Name', 'Cross-Animal Epoch Average (Lights On vs Off)', ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [200 200 700 400]);
    ax = axes(fCross, 'Units', 'normalized', 'Position', [0.12 0.15 0.82 0.75]);
    hold(ax, 'on');

    t_stim = 0:restPeriod;
    t_cont = restPeriod:(2*restPeriod);
    fill(ax, [t_stim fliplr(t_stim)], ...
        100*[pooledStimMean+pooledStimSEM fliplr(pooledStimMean-pooledStimSEM)], ...
        'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hStim = plot(ax, t_stim, 100*pooledStimMean, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Lights On');
    fill(ax, [t_cont fliplr(t_cont)], ...
        100*[pooledContMean+pooledContSEM fliplr(pooledContMean-pooledContSEM)], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hCont = plot(ax, t_cont, 100*pooledContMean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Lights Off');

    hold(ax, 'off');
    xlim(ax, [0 2*restPeriod]);
    ylim(ax, [-5 5]);
    xlabel(ax, 'Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(ax, 'Change in LGN Blood Volume (%)', 'FontSize', 18, 'FontWeight', 'bold');
    title(ax, sprintf('Epoch Average Across Animals (N=%d)', nAnimals), ...
        'FontSize', 20, 'FontWeight', 'bold');
    legend(ax, [hStim, hCont], {'Lights On', 'Lights Off'}, ...
        'Location', 'northeast', 'FontSize', 14, 'FontWeight', 'bold');
    set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);

    fprintf('  Pooled t-test: p=%.4f\n', pPooled);

    % --- Effect size and Bonferroni correction ---
    pooledSD_xa = sqrt((var(pooledStimEpochMean) + var(pooledContEpochMean)) / 2);
    if pooledSD_xa > 0
        cohensD_xa = (mean(pooledStimEpochMean) - mean(pooledContEpochMean)) / pooledSD_xa;
    else
        cohensD_xa = NaN;
    end
    nComparisons = nAnimals + 1;  % per-subject tests + pooled
    pBonf_xa = min(pPooled * nComparisons, 1.0);
    fprintf('  Pooled Cohen''s d: %.3f\n', cohensD_xa);
    fprintf('  Bonferroni-corrected p (k=%d): %.4f\n', nComparisons, pBonf_xa);

    % Per-subject summary
    fprintf('\n  --- Per-Subject Effect Sizes ---\n');
    for si = 1:nAnimals
        ed = ud.animalStore{si}.epochData;
        sStim = mean(ed.dataStim, 2);
        sCont = mean(ed.dataCont, 2);
        sdPool = sqrt((var(sStim) + var(sCont)) / 2);
        if sdPool > 0
            dSubj = (mean(sStim) - mean(sCont)) / sdPool;
        else
            dSubj = NaN;
        end
        [~, pSubj] = ttest(sCont, sStim);
        pBonf_subj = min(pSubj * nComparisons, 1.0);
        fprintf('  %s: d=%.3f | p=%.4f (uncorrected) | p=%.4f (Bonferroni, k=%d)\n', ...
            ed.subject, dSubj, pSubj, pBonf_subj, nComparisons);
    end
    fprintf('  Pooled: d=%.3f | p=%.4f (uncorrected) | p=%.4f (Bonferroni, k=%d)\n', ...
        cohensD_xa, pPooled, pBonf_xa, nComparisons);

    % --- Generate combined bar plot ---
    barTitle = sprintf('Mean CBV Across Animals (N=%d)', nAnimals);
    generateEpochBarPlot(pooledStimEpochMean, pooledContEpochMean, barTitle);

    fprintf('=== Cross-animal figure generated ===\n\n');
end

%% ===================== EPOCH BAR PLOT =====================
function fBar = generateEpochBarPlot(stimEpochMeans, contEpochMeans, titleStr)
    % Bar plot comparing mean CBV for lights on vs lights off.
    % stimEpochMeans / contEpochMeans: per-epoch mean rBV vectors.

    stimGrandMean = mean(stimEpochMeans) * 100;
    contGrandMean = mean(contEpochMeans) * 100;
    nEpochs = length(stimEpochMeans);
    stimSEM = std(stimEpochMeans) / sqrt(nEpochs) * 100;
    contSEM = std(contEpochMeans) / sqrt(nEpochs) * 100;

    [~, pVal] = ttest(contEpochMeans, stimEpochMeans);

    fBar = figure('Name', titleStr, 'NumberTitle', 'off', 'Color', 'w', ...
        'Position', [250 250 500 500]);
    ax = axes(fBar);
    hold(ax, 'on');

    bar(ax, 1, stimGrandMean, 0.6, 'FaceColor', [0.85 0.2 0.2], ...
        'EdgeColor', 'k', 'LineWidth', 1.5);
    bar(ax, 2, contGrandMean, 0.6, 'FaceColor', [0.2 0.4 0.85], ...
        'EdgeColor', 'k', 'LineWidth', 1.5);

    errorbar(ax, 1, stimGrandMean, stimSEM, 'k', 'LineWidth', 1.5, 'CapSize', 12);
    errorbar(ax, 2, contGrandMean, contSEM, 'k', 'LineWidth', 1.5, 'CapSize', 12);

    % Significance bracket and annotation
    allVals = [stimGrandMean + stimSEM, stimGrandMean - stimSEM, ...
               contGrandMean + contSEM, contGrandMean - contSEM];
    yMax = max(allVals);
    yRange = max(allVals) - min(allVals);
    if yRange == 0, yRange = 1; end
    yBar  = yMax + 0.12 * yRange;
    yTick = yMax + 0.08 * yRange;
    yText = yBar + 0.02 * yRange;

    plot(ax, [1 2], [yBar yBar], 'k-', 'LineWidth', 1.5);
    plot(ax, [1 1], [yTick yBar], 'k-', 'LineWidth', 1.5);
    plot(ax, [2 2], [yTick yBar], 'k-', 'LineWidth', 1.5);

    if pVal < 0.001
        starStr = '***';
    elseif pVal < 0.01
        starStr = '**';
    elseif pVal < 0.05
        starStr = '*';
    else
        starStr = 'n.s.';
    end
    text(ax, 1.5, yText, starStr, ...
        'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');

    hold(ax, 'off');
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Lights On', 'Lights Off'});
    ylabel(ax, 'Change in LGN Blood Volume (%)', 'FontSize', 18, 'FontWeight', 'bold');
    title(ax, titleStr, 'FontSize', 20, 'FontWeight', 'bold');
    set(ax, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5);
    xlim(ax, [0.3 2.7]);
    ylim(ax, [-5 5]);

    fprintf('  Bar plot generated: %s (p=%.4f)\n', titleStr, pVal);
end

%% ===================== DUAL EXPORT: PDF + SVG =====================
function exportPdfSvg(figHandle, pdfPath, varargin)
    % Export as vector PDF
    exportgraphics(figHandle, pdfPath, 'ContentType', 'vector', varargin{:});
    fprintf('  Saved: %s\n', pdfPath);
    % Export as SVG
    svgPath = strrep(pdfPath, '.pdf', '.svg');
    print(figHandle, '-dsvg', svgPath);
    fprintf('  Saved: %s\n', svgPath);
end

%% ===================== BUILD Dop.bin FROM BEAMFORMED DATA =====================
function generateDopBin(fusiDir)
% Power Doppler from beamformed IQ, saved as Dop.bin in fusiDir.
% Same calculation as the first half of fUS_correlationZnorm.m.

    seqParamPath = dir(fullfile(fusiDir, 'SequenceParam*'));
    if isempty(seqParamPath)
        error('Cannot build Dop.bin. No SequenceParam file in %s.', fusiDir);
    end
    SP = load(fullfile(fusiDir, seqParamPath(1).name));

    fileName = SP.UserParam.savePrefix;
    svdVal   = [70, SP.UserParam.NumFrames];
    sizeZ    = SP.CParam.PData.Size(1);
    sizeX    = SP.CParam.PData.Size(2);
    sizeY    = SP.CParam.PData.Size(3);   % 1 for 2D; kept for Vantage NXT volumes
    nFrames  = SP.UserParam.NumFrames;
    nLoops   = SP.UserParam.NumLoops;

    Dop = zeros(sizeZ, sizeX, sizeY, nLoops + 1);

    fprintf('  Computing Dop.bin (%d loops)...\n', nLoops);
    for iAcq = 1:nLoops
        bfFile = fullfile(fusiDir, [fileName, '_BF_', num2str(iAcq, '%03.3d'), '.bin']);
        if ~isfile(bfFile)
            error('Cannot build Dop.bin. Missing beamformed file:\n  %s', bfFile);
        end
        fprintf('    Loop %d/%d\n', iAcq, nLoops);
        Img = readBinFile(bfFile, [sizeZ, sizeX, sizeY, nFrames], 'complex');
        ImgFilt = SVDfilter(Img, svdVal);
        tmpDop = sqrt(mean(abs(ImgFilt).^2, 4));
        Dop(:,:,iAcq) = (tmpDop - min(tmpDop(:))) ./ max(tmpDop(:));
    end

    dopFile = fullfile(fusiDir, 'Dop.bin');
    fid = fopen(dopFile, 'w');
    if fid < 0
        error('Could not open %s for writing.', dopFile);
    end
    fwrite(fid, Dop, 'double');
    fclose(fid);
    fprintf('  Saved %s\n', dopFile);
end


%% ===================== BUILD Dop.bin USING BUTTERWORTH FILTERING =====================
function generateDopBin2(fusiDir)
% Generate power Doppler images from beamformed complex IQ data.
%
% Processing:
%   1. Load SequenceParam metadata.
%   2. Read each beamformed IQ file.
%   3. Subtract the first IQ frame.
%   4. Apply a zero-phase Butterworth high-pass filter.
%   5. Sum squared IQ magnitudes across time to calculate power Doppler.
%   6. Normalize each power Doppler image.
%   7. Save the image series as Dop.bin.
%
% Output:
%   Dop.bin: double-precision array with dimensions
%            [sizeZ, sizeX, sizeY, nLoops+1].
%
% The final frame is left as zeros to preserve the existing output format.

    % ---------------------------------------------------------
    % Butterworth parameters
    % ---------------------------------------------------------
    IQ_RATE_HZ = 250;
    BUTTERWORTH_CUTOFF_HZ = 30;
    BUTTERWORTH_ORDER = 4;

    % ---------------------------------------------------------
    % Load acquisition metadata
    % ---------------------------------------------------------
    seqParamPath = dir(fullfile(fusiDir, 'SequenceParam*'));

    if isempty(seqParamPath)
        error('No SequenceParam file found in %s.', fusiDir);
    end

    SP = load(fullfile(fusiDir, seqParamPath(1).name));

    fileName = SP.UserParam.savePrefix;

    sizeZ = SP.CParam.PData.Size(1);
    sizeX = SP.CParam.PData.Size(2);
    sizeY = SP.CParam.PData.Size(3);

    nFrames = SP.UserParam.NumFrames;
    nLoops  = SP.UserParam.NumLoops;

    % ---------------------------------------------------------
    % Design Butterworth high-pass filter in SOS form
    % ---------------------------------------------------------
    Wn = BUTTERWORTH_CUTOFF_HZ / (IQ_RATE_HZ / 2);

    [z, p, k] = butter( ...
        BUTTERWORTH_ORDER, Wn, 'high');

    [sos, g] = zp2sos(z, p, k);

    fprintf(['Butterworth high-pass: order=%d, cutoff=%g Hz, ' ...
             'sampling=%g Hz\n'], ...
        BUTTERWORTH_ORDER, ...
        BUTTERWORTH_CUTOFF_HZ, ...
        IQ_RATE_HZ);

    % ---------------------------------------------------------
    % Allocate output array
    % ---------------------------------------------------------
    Dop = zeros(sizeZ, sizeX, sizeY, nLoops + 1);

    fprintf('Computing Dop.bin (%d loops)...\n', nLoops);

    % ---------------------------------------------------------
    % Process each acquisition loop
    % ---------------------------------------------------------
    for iAcq = 1:nLoops

        bfFile = fullfile(fusiDir, ...
            [fileName, '_BF_', num2str(iAcq, '%03.3d'), '.bin']);

        if ~isfile(bfFile)
            error('Missing beamformed file:\n  %s', bfFile);
        end

        fprintf('  Loop %d/%d\n', iAcq, nLoops);

        % Read complex beamformed IQ data:
        % [depth, lateral, plane, time]
        Img = readBinFile( ...
            bfFile, ...
            [sizeZ, sizeX, sizeY, nFrames], ...
            'complex');

        % -----------------------------------------------------
        % Subtract first IQ frame
        % -----------------------------------------------------
        firstFrame = Img(:,:,:,1);

        Img = bsxfun(@minus, Img, firstFrame);

        % -----------------------------------------------------
        % Reshape to [spatial pixels, time]
        % -----------------------------------------------------
        matrix = reshape(Img, [], nFrames);

        % -----------------------------------------------------
        % Butterworth high-pass filtering along time
        % -----------------------------------------------------
        % Transpose so time is the first dimension for filtfilt.
        % Filtering is performed independently for every pixel.
        filtered = filtfilt(sos, g, matrix.');

        filtered = filtered.';

        % Restore [depth, lateral, plane, time]
        ImgFilt = reshape( ...
            filtered, [sizeZ, sizeX, sizeY, nFrames]);

        % -----------------------------------------------------
        % Calculate power Doppler by summing over time
        % -----------------------------------------------------
        tmpDop = sum(abs(ImgFilt).^2, 4);

        % -----------------------------------------------------
        % Normalize each loop's power Doppler image
        % -----------------------------------------------------
        tmpMin = min(tmpDop(:));
        tmpMax = max(tmpDop(:));

        if tmpMax > tmpMin
            tmpDop = (tmpDop - tmpMin) ./ (tmpMax - tmpMin);
        else
            tmpDop = zeros(size(tmpDop));
        end

        % Store this loop's power Doppler image
        Dop(:,:,:,iAcq) = tmpDop;

    end

    % ---------------------------------------------------------
    % Write Dop.bin
    % ---------------------------------------------------------
    dopFile = fullfile(fusiDir, 'Dop.bin');

    fid = fopen(dopFile, 'w');

    if fid < 0
        error('Could not open %s for writing.', dopFile);
    end

    fwrite(fid, Dop, 'double');
    fclose(fid);

    fprintf('Saved %s\n', dopFile);

end


function IQf = SVDfilter(IQ, cutoff)
% SVD clutter filter. Temporal dimension must be the last dimension.
% Same implementation as fUS_correlationZnorm.m.

    initsize = size(IQ);
    if cutoff(end) > initsize(end)
        cutoff = cutoff(1):initsize(end);
    end
    if numel(cutoff) == 1
        cutoff = cutoff(1):initsize(end);
    elseif numel(cutoff) == 2
        cutoff = cutoff(1):cutoff(2);
    end

    if or(isequal(cutoff, 1:size(IQ, 3)), cutoff(1) < 2)
        IQf = IQ;
        return
    end

    X = reshape(IQ, prod(initsize(1:end-1)), initsize(end));
    [U, ~] = svd(X' * X);
    V = X * U;
    Reconst = V(:, cutoff) * U(:, cutoff)';
    IQf = reshape(Reconst, initsize);
end

%% ===================== NORMALIZE MRI SLICE (0-1, NaN -> 0) =====================
function mrN = normalizeMriSlice(mr)
    mrN = mr; vpMR = ~isnan(mrN);
    if any(vpMR(:))
        lo = min(mrN(vpMR)); hi = max(mrN(vpMR));
        if hi > lo, mrN = (mrN - lo) / (hi - lo); end
        mrN(isnan(mrN)) = 0;
    else
        mrN = zeros(size(mr));
    end
end