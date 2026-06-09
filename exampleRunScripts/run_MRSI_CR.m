%% rosette_r14.m  —  CONCENTRIC-RINGS reconstruction pipeline
%   Dataset: F:\fida\divya\20260605_phantom_test\subject03  (CR_40x40 phantom)
%   Trajectory: conc_ring_anticlockwise_posY_FOV240_grid40_20rings.txt
%
%   The scan was played with the XA60_RosetteSpinEcho multi-trajectory
%   sequence but using a CONCENTRIC-RINGS gradient, so the loader labels it
%   "rosette".  That is fine
%   Folder layout (file_path_load resolves / creates anything missing):
%       subject03/
%           met/        CR_40x40_isoctr.dat        (water-suppressed)
%           mrs_ref/    CR_40x40_isoctr_w.dat      (water reference)
%           dicom/      MPRAGE etc -> T1 auto-built into img_ref/
%           outputs/    empty + filled NIfTIs

%% === USER INPUTS ===========================================================
datasetFolder = 'F:\fida\divya\20260605_phantom_test\subject03';
kFile         = 'C:\Users\divya\Downloads\fida codes\fid_a\processingTools\MRSI\kFiles\conc_ring_anticlockwise_posY_FOV240_grid40_20rings.txt';

% Reconstruction choices
dcfMethod     = 'pipe_menon';   % 'nn' | 'voronoi' | 'pipe_menon' | 'none'
ftMethod      = 'nufft';        % 'nufft' | 'dft' | 'tikhonov'
smoothFwhm    = 20;             % spatial Gaussian FWHM (mm)

% LCModel — left EMPTY so the pipeline stops after recon + overlay.
% To fit later you need a basis matched to THIS acquisition's DELTAT/NUNFIL/TE
% (auto-derived in step 8.6); set lcmbin + basisFile then.

lcmOwner      = 'Divyasri Krishnkaumar, Sunnybrook Research Institute';
hzpppm        = 123.25;
deltat        = 0.00063;   % placeholders; auto-derived in step 8.6
nunfil        = 576;
echot         = 15;
ppmst         = 4.25;
ppmend        = 0.2;

%% === RESOLVE PATHS ========================================================
assert(exist(kFile,'file')==2, 'k-file missing: %s', kFile);
paths = file_path_load(datasetFolder);   % builds T1 from dicom/ if img_ref/ empty

%% === 1. LOAD ===============================================================
fprintf('--- 1. Loading TWIX pair ---\n');
[timeCombined_rs, timeCombined_rs_w] = io_CSIload_twix_pair( ...
    paths.metFile, paths.refFile, kFile);
fprintf('   data size: %s   seq: %s\n', mat2str(timeCombined_rs.sz), timeCombined_rs.seq);

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

%% === 8.6. AUTO-DERIVE LCMODEL TIMING FROM THE DATA =======================
specDim = ftSpec_smooth.dims.f;
if specDim == 0, specDim = ftSpec_smooth.dims.t; end
deltat = 1 / ftSpec_smooth.spectralWidth;
nunfil = ftSpec_smooth.sz(specDim);
echot  = ftSpec_smooth.te;
hzpppm = ftSpec_smooth.txfrq / 1e6; %this is taken somewhere

lcmbin        = '/home/divya/.lcmodel/bin/lcmodel';   % set '' to skip LCModel
basisFile     = 'C:\Users\divya\Downloads\fida codes\fid_a\basis_r14_phantom\LacAceCrCho_TE15_SW4000_NP4096.basis';

%% === 9. LCMODEL + METABOLITE MAPS =========================================
if isempty(lcmbin)
    fprintf('\n--- 9. LCModel: SKIPPED (lcmbin empty) ---\n');
    fprintf('To fit: build a basis matched to DELTAT=%.5g / NUNFIL=%d / TE=%g\n', ...
            deltat, nunfil, echot);
    fprintf('\n=== Reconstruction complete.  Outputs in %s ===\n', paths.outputDir);
    return
end

fprintf('\n--- 9. LCModel + maps ---\n');
load('SNR.mat')
load('map.mat')
load('crlb.mat')
load('LW.mat')

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

save('SNR.mat', 'SNR', '-v7.3');
save('map.mat', 'map', '-v7.3');
save('crlb.mat', 'crlb', '-v7.3');
save('LW.mat', 'LW', '-v7.3');


out = create_separate_metabolite_niftis_v2(ftSpec_smooth, map, crlb, LW, SNR, ...
        paths.t1Path, paths.emptyNiiPath, paths.mapsDir);
nii_viewer(paths.t1Path, out.metabolite_4d_files);
fprintf('\n=== Pipeline complete.  Outputs in %s ===\n', paths.outputDir);
