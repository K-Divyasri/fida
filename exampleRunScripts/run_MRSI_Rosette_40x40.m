datasetFolder = 'F:\fida\divya\newPhantom_20251119 - Copy\subject07';
kFile         = 'C:\Users\divya\Downloads\fida codes\fid_a\processingTools\MRSI\kFiles\Rosette_traj_40x40.txt';

% Reconstruction choices
dcfMethod     = 'pipe_menon';        % 'nn' | 'voronoi' | 'pipe_menon' | 'none'
ftMethod      = 'nufft';     % 'nufft' | 'dft' | 'tikhonov'   (tikhonov auto-skips DCF)

% Off-resonance MFI recon (corrects intra-readout dephasing that
% underestimates metabolites far from the water carrier on the rosette
% trajectory). Only used with ftMethod='nufft'. Validate PER DATASET:
% run once false (baseline), once true, compare to GT. Flip offresSign if
% the far-from-water peaks get WORSE.
useOffres  = false;   % baseline first for this 2-bottle set; flip to test MFI
nMFI       = 5;
offresSign = +1;

smoothFwhm    = 0;            % spatial Gaussian FWHM (mm). 0 = NO smoothing.
                             % TWO-bottle phantom: 20 mm blends inner<->outer
                             % (partial volume) and wrecks per-compartment mM.
                             % MUST be 0 for quantification here.

% LCModel — leave lcmbin empty to skip LCModel + map generation

lcmOwner      = 'Divyasri Krishnkaumar, Sunnybrook Research Institute';
hzpppm        = 123.25;      % PLACEHOLDER ONLY - overwritten from ftSpec_smooth.txfrq before step 9
deltat        = 0.00063;     % matches rosette data SW 1587 Hz
nunfil        = 576;         % matches rosette data point count
echot         = 15;          % matches rosette data TE
ppmst         = 4.25;
ppmend        = 0.2;

%% === RESOLVE PATHS ========================================================
assert(exist(kFile,'file')==2, 'k-file missing: %s', kFile);
paths = file_path_load(datasetFolder);

%% === 1. LOAD ===============================================================
fprintf('--- 1. Loading TWIX pair ---\n');
[timeCombined_rs, timeCombined_rs_w] = io_CSIload_twix_pair( ...
    paths.metFile, paths.refFile, kFile);

%% === 2. SPATIAL RECONSTRUCTION (DCF + FT) ==================================
fprintf('\n--- 2. Spatial reconstruction (dcf=%s, ft=%s) ---\n', dcfMethod, ftMethod);
ft   = op_CSIRecon(timeCombined_rs,   kFile, dcfMethod, ftMethod, ...
    'offres', useOffres, 'nMFI', nMFI, 'offresSign', offresSign);
ft_w = op_CSIRecon(timeCombined_rs_w, kFile, dcfMethod, ftMethod, ...
    'offres', useOffres, 'nMFI', nMFI, 'offresSign', offresSign);

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

%% === 4b. LEFT-SHIFT (time domain, before spectral FT) =====================
% Left-shift removes the FID first-point 1st-order phase (ADC/echo delay).
% NOTE: zero-fill moved to AFTER B0 correction (step 6) so lipid removal and
% B0 run on the short native FID (fast, no zero-padded region to distort);
% only the final spectrum is interpolated. See op_CSIspecZeroFill.
fprintf('\n--- 4b. Left-shift ---\n');
ccav   = op_CSIleftshift(ccav);            % uses ccav.pointsToLeftshift
ccav_w = op_CSIleftshift(ccav_w);

%% === 5. SPECTRAL FT =======================================================
fprintf('\n--- 5. Spectral FT ---\n');
ftSpec   = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);
save('ftSpec.mat', 'ftSpec', '-v7.3');
save('ftSpec_w.mat', 'ftSpec_w', '-v7.3');

%% === 6. LIPID + WATER REMOVAL + B0 CORRECTION =============================
fprintf('\n--- 6. Lipid / water removal + B0 correction ---\n');
%ftSpec_rmlip = op_CSIssp(ftSpec, 0.8, 1.88);
ftSpec_rmw   = op_CSIRemoveLipids(ftSpec, ...
                  'lipidPPMRange',  [4.5 5.0], ...
                  'linewidthRange', [1 10]);
% tmp        = op_CSIRemoveWater(ftSpec,[4.4 5.0],50);
% ftSpec_rmw = op_CSIRemoveLipids(tmp,'lipidPPMRange',[4.4 5.0],'linewidthRange',[5 40],'beta',5e-4);
[ftSpec_B0corr, ftSpec_B0corr_w, freqMap, R2Map] = ...
    op_CSIB0Correction_v2(ftSpec_rmw, ftSpec_w);

% Zero-fill NOW (after B0), interpolating the processed spectrum 576 -> 4096
% for finer ppm sampling / better LCModel peak separation. Padding the native
fprintf('  Zero-fill (spectral) 576 -> 4096 ...\n');
ftSpec_B0corr   = op_CSIspecZeroFill(ftSpec_B0corr,   4096);
ftSpec_B0corr_w = op_CSIspecZeroFill(ftSpec_B0corr_w, 4096);

ftSpec_masked = op_CSIapplymask(ftSpec_B0corr);

%% === 7. SPATIAL SMOOTHING =================================================
% op_CSIApodize divides by sigma, so FWHM=0 would blow up. Skip it entirely
% when smoothFwhm<=0 and pass the un-smoothed spectra straight through (the
% variables keep the _smooth name so steps 8-9 need no change).
if smoothFwhm > 0
    fprintf('\n--- 7. Spatial smoothing (Gaussian, FWHM=%g mm) ---\n', smoothFwhm);
    ftSpec_smooth   = op_CSIApodize(ftSpec_masked,  ...
                           'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
    ftSpec_smooth_w = op_CSIApodize(ftSpec_B0corr_w, ...
                           'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
else
    fprintf('\n--- 7. Spatial smoothing: SKIPPED (smoothFwhm=%g) ---\n', smoothFwhm);
    ftSpec_smooth   = ftSpec_masked;
    ftSpec_smooth_w = ftSpec_B0corr_w;
end
op_CSIPlot(ftSpec_smooth);

%% === 8. WRITE 4D NIFTI + OPEN VIEWER ======================================
fprintf('\n--- 8. Writing 4D MRSI NIfTI ---\n');
mrsiout = mrsi_on_t1_map(ccav_w, ftSpec_smooth, paths.mrsiOut, paths.emptyNiiPath);
nii_viewer(paths.t1Path, mrsiout.mrsi4D_time);

lcmbin        = '/home/divya/.lcmodel/bin/lcmodel';   % e.g. '/home/divya/.lcmodel/bin/lcmodel'
% Basis rebuilt at the data dwell (BADELT = DELTAT = 6.3e-4 s, SW 1587 Hz)
% so LCModel no longer resamples it (kills the MYBASI 8 warning).
basisFile     = 'C:\Users\divya\Downloads\fida codes\fid_a\basis_rosette_phantom\RosettePhantom_LacAceCrCho_TE15_SW1587.basis';



%% === 9. LCMODEL + METABOLITE MAPS =========================================
if isempty(lcmbin)
    fprintf('\n--- 9. LCModel: SKIPPED (lcmbin empty) ---\n');
    fprintf('\n=== Pipeline complete.  Outputs in %s ===\n', paths.outputDir);
    return
end

% Auto-derive LCModel timing from the data so the control matches the basis.
% The hard-coded hzpppm=123.25 above was WRONG (real value 123.19877) and was
% the sole remaining cause of the MYBASI 8 warning once BADELT was fixed.
specDim = ftSpec_smooth.dims.f;
if specDim == 0, specDim = ftSpec_smooth.dims.t; end
hzpppm = ftSpec_smooth.txfrq / 1e6;
deltat = 1 / ftSpec_smooth.spectralWidth;
nunfil = ftSpec_smooth.sz(specDim);
echot  = ftSpec_smooth.te;
fprintf('  LCModel timing (from data): HZPPPM=%.5f  DELTAT=%.8f  NUNFIL=%d  ECHOT=%g\n', ...
        hzpppm, deltat, nunfil, echot);

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
    'ppmend',     ppmend, ...
    'doecc',      false);   % LCModel eddy-current OFF: op_CSIB0Correction_v2
                            % already did the water-phase (Klose) correction
                            % in-pipeline, and DOECC=T would divide the metab
                            % FID by the zero-filled water tail (0/0 -> NaN ->
                            % SOLVE FATAL). DOWS/water scaling stays ON.

Nx = numel(getCoordinates(ftSpec, 'x'));
Ny = numel(getCoordinates(ftSpec, 'y'));
[map, crlb, LW, SNR] = op_CSILCModelMaps(Nx, Ny, paths.rawFolder, ...
    'figure_folder_name', 'maps', ...
    'metabList', {'Lac' 'Act' 'Cr' 'Cho'});   % this phantom's basis labels
  % was apply_water_scaling_all
save('map.mat', 'map', '-v7.3');
save('crlb.mat', 'crlb', '-v7.3');
save('LW.mat', 'LW', '-v7.3');
save('SNR.mat', 'SNR', '-v7.3');
apply_water_scaling_all; 

out = create_separate_metabolite_niftis_v2(ftSpec_smooth, map, crlb, LW, SNR, ...
        paths.t1Path, paths.emptyNiiPath, paths.mapsDir);
%nii_viewer(paths.t1Path, paths.mapsDir);
nii_viewer(paths.t1Path, out.metabolite_4d_files);
allMaps = [out.metabolite_4d_files, out.crlb_3d_files,{out.shared_lw_3d_file, out.shared_snr_3d_file}];
nii_viewer(paths.t1Path, allMaps);
fprintf('\n=== Pipeline complete.  Outputs in %s ===\n', paths.outputDir);
