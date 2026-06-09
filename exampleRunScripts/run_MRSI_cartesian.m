%% rosette_pipeline.m
%   Folder layout — file_path_load resolves and creates anything missing:
%       <datasetFolder>/
%           met/        single .dat (water-suppressed)
%           mrs_ref/    single .dat (water reference)
%           img_ref/    T1 NIfTI — auto-built from /dicom/ if empty
%           dicom/      DICOM series (only used if /img_ref/ is empty)
%           outputs/    created; holds empty + filled NIfTIs
%               lcm_out/    LCModel raw + .ps outputs
%               maps/       per-metabolite NIfTIs, SNR, CRLB


%% === USER INPUTS ===========================================================
datasetFolder = 'F:\fida\divya\20260605_phantom_test\subject04';

% CARTESIAN run: leave kFile EMPTY.  No trajectory, no density compensation,
% plain FFT recon (op_CSIRecon auto-detects empty kFile / isCartesian flag).
% For a rosette/concentric run, set kFile to the trajectory .txt instead.
kFile         = '';

% Reconstruction choices  (dcfMethod / ftMethod are IGNORED for Cartesian)
dcfMethod     = 'pipe_menon';        % 'nn' | 'voronoi' | 'pipe_menon' | 'none'
ftMethod      = 'nufft';     % 'nufft' | 'dft' | 'tikhonov'   (tikhonov auto-skips DCF)
smoothFwhm    = 20;          % spatial Gaussian FWHM (mm)

lcmbin        = '/home/divya/.lcmodel/bin/lcmodel';   % set '' to skip LCModel
basisFile     = 'C:\Users\divya\Downloads\fida codes\fid_a\basis_phantom\Phantom_LacAceCrCho_TE2p3_SW4000_2048.basis';  % Lac/Ace/Cr/Cho phantom basis
lcmOwner      = 'Divyasri Krishnkaumar, Sunnybrook Research Institute';
hzpppm        = 123.25;      % 1H @ ~2.89 T -> 123.25 MHz; auto-checked below
deltat        = 0.00025;     % placeholder; auto-derived below (1/spectralWidth)
nunfil        = 2048;        % placeholder; auto-derived below (# spectral points)
echot         = 2.3;         % placeholder; auto-derived below (ms)
ppmst         = 4.25;
ppmend        = 0.2;

%% === RESOLVE PATHS ========================================================
% kFile required only for non-Cartesian.  Empty kFile => Cartesian (FFT, no DCF).
isCartesian = isempty(kFile);
if ~isCartesian
    assert(exist(kFile,'file')==2, 'k-file missing: %s', kFile);
else
    fprintf('No k-file supplied -> CARTESIAN mode (FFT recon, no DCF).\n');
end
paths = file_path_load(datasetFolder);

%% === 1. LOAD ===============================================================
fprintf('--- 1. Loading TWIX pair ---\n');
[timeCombined_rs, timeCombined_rs_w] = io_CSIload_twix_pair( ...
    paths.metFile, paths.refFile, kFile);

%% === 2. SPATIAL RECONSTRUCTION (DCF + FT) ==================================
fprintf('\n--- 2. Spatial reconstruction (dcf=%s, ft=%s) ---\n', dcfMethod, ftMethod);
ft   = op_CSIRecon(timeCombined_rs,   kFile, dcfMethod, ftMethod);
ft_w = op_CSIRecon(timeCombined_rs_w, kFile, dcfMethod, ftMethod);

%% === 3. COIL COMBINATION ==================================================
fprintf('\n--- 3. Coil combination ---\n');
[coilCombined_w, phase, weights] = op_CSICombineCoils1(ft_w);
coilCombined                     = op_CSICombineCoils1(ft, 1, phase, weights);

%% === 4. AVERAGE + WATER MASK ==============================================
fprintf('\n--- 4. Averaging + water mask ---\n');
ccav   = op_CSIAverage(coilCombined);
ccav_w = op_CSIAverage(coilCombined_w);
ccav_w = op_CSISegment_simple(ccav_w);
%ccav_w.mask.brainmasks=ccav_w.mask.brainmask;
ccav.mask = ccav_w.mask;

mask   = ccav_w.mask.brainmasks;

%% === 5. SPECTRAL FT =======================================================
fprintf('\n--- 5. Spectral FT ---\n');
ftSpec   = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);

%% === 6. LIPID + WATER REMOVAL + B0 CORRECTION =============================
fprintf('\n--- 6. Lipid / water removal + B0 correction ---\n');
ftSpec_rmlip = op_CSIssp(ftSpec, 0.8, 1.88);
ftSpec_rmw   = op_CSIRemoveLipids(ftSpec_rmlip, ...
                  'lipidPPMRange',  [4.5 5.0], ...
                  'linewidthRange', [1 10]);
[ftSpec_B0corr, ftSpec_B0corr_w, freqMap, R2Map] = ...
    op_CSIB0Correction_v2(ftSpec_rmw, ftSpec_w);

ftSpec_masked = op_CSIapplymask(ftSpec_B0corr);

%% === 7. SPATIAL SMOOTHING =================================================
fprintf('\n--- 7. Spatial smoothing (Gaussian, FWHM=%g mm) ---\n', smoothFwhm);
ftSpec_smooth   = op_CSIApodize(ftSpec_masked,  ...
                       'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
ftSpec_smooth_w = op_CSIApodize(ftSpec_B0corr_w, ...
                       'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
op_CSIPlot(ftSpec_smooth);

%% === 8. WRITE 4D NIFTI + OPEN VIEWER ======================================
fprintf('\n--- 8. Writing 4D MRSI NIfTI ---\n');
mrsiout = mrsi_on_t1_map(ccav_w, ftSpec_smooth, paths.mrsiOut, paths.emptyNiiPath);
nii_viewer(paths.t1Path, mrsiout.mrsi4D_time);

specDim = ftSpec_smooth.dims.f;
if specDim == 0, specDim = ftSpec_smooth.dims.t; end   % pre/post spectral-FT safety

deltat = 1 / ftSpec_smooth.spectralWidth;              % spectral dwell time (s)
nunfil = ftSpec_smooth.sz(specDim);                    % # spectral points
echot  = ftSpec_smooth.te;                             % echo time (ms)
hzpppm = ftSpec_smooth.txfrq / 1e6;                    % MHz per ppm

fprintf('\n--- 8.6. LCModel timing (auto-derived) ---\n');
fprintf('  DELTAT = %.8f s  (SW = %.1f Hz)\n', deltat, ftSpec_smooth.spectralWidth);
fprintf('  NUNFIL = %d points\n', nunfil);
fprintf('  ECHOT  = %g ms\n', echot);
fprintf('  HZPPPM = %.4f MHz/ppm\n', hzpppm);
fprintf('  BASIS  = %s\n', basisFile);


%% === 9. LCMODEL + METABOLITE MAPS =========================================
if isempty(lcmbin)
    fprintf('\n--- 9. LCModel: SKIPPED (lcmbin empty) ---\n');
    fprintf('\n=== Pipeline complete.  Outputs in %s ===\n', paths.outputDir);
    return
end

fprintf('\n--- 9. LCModel + maps ---\n');
run_lcm_rosette_portable(paths.metFile, ftSpec_smooth, ftSpec_smooth_w, mask, ...
    'basisFile',  basisFile, ...
    'lcmodelBin', lcmbin, ...
    'ownerStr',   lcmOwner, ...
    'rawFolder',  paths.rawFolder, ...
    'hzpppm',     hzpppm, ...
    'deltat',     deltat, ...
    'nunfil',     nunfil, ...
    'echot',      echot, ...
    'ppmst',      ppmst, ...
    'ppmend',     ppmend);

Nx = numel(getCoordinates(ftSpec, 'x'));
Ny = numel(getCoordinates(ftSpec, 'y'));
[map, crlb, LW, SNR] = op_CSILCModelMaps(Nx, Ny, paths.rawFolder, ...
    'figure_folder_name', 'maps');
save('map.mat', 'map', '-v7.3');
save('crlb.mat', 'crlb', '-v7.3');
save('LW.mat', 'LW', '-v7.3');
save('SNR.mat', 'SNR', '-v7.3');
save('map.mat', 'map', '-v7.3');
save('crlb.mat', 'crlb', '-v7.3');
save('LW.mat', 'LW', '-v7.3');
save('SNR.mat', 'SNR', '-v7.3');
save('SNR.mat', 'SNR', '-v7.3');
out = create_separate_metabolite_niftis_v2(ftSpec_smooth, map, crlb, LW, SNR, ...
        paths.t1Path, paths.emptyNiiPath, paths.mapsDir);
%nii_viewer(paths.t1Path, paths.mapsDir);
nii_viewer(paths.t1Path, out.metabolite_4d_files);
allMaps = [out.metabolite_4d_files, out.crlb_3d_files,{out.shared_lw_3d_file, out.shared_snr_3d_file}];
nii_viewer(paths.t1Path, allMaps);
fprintf('\n=== Pipeline complete.  Outputs in %s ===\n', paths.outputDir);
