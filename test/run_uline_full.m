%% run_MRSI_Rosette_40x40.m  —  rosette phantom, straight run, NO MFI.
% load -> recon -> combine -> average -> [2 Hz filter -> zero-pad 4x] -> spectral
% FT -> [1st-order phase 1.2e-4 s] -> water/lipid removal + B0 -> LCModel -> maps.
% Prof preprocessing (step 4c + phase after FT): apodize-then-zeropad avoids the
% truncation ripple that a bare spectral zero-fill caused.

datasetFolder = 'F:\fida\divya\20260831_ulinePhantom_test\subject04';
kFile         = 'C:\Users\divya\Downloads\fida codes\fid_a\processingTools\MRSI\kFiles\Rosette_traj_40x40.txt';

dcfMethod = 'pipe_menon';    % 'nn' | 'voronoi' | 'pipe_menon' | 'none'
ftMethod  = 'nufft';         % 'nufft' | 'dft' | 'tikhonov'

lcmOwner  = 'Divyasri Krishnkaumar, Sunnybrook Research Institute';
ppmst     = 4.25;
ppmend    = 0;

%% === RESOLVE PATHS ========================================================
assert(exist(kFile,'file')==2, 'k-file missing: %s', kFile);
paths = file_path_load(datasetFolder);

%% === 1. LOAD ==============================================================
fprintf('--- 1. Loading TWIX pair ---\n');
[timeCombined_rs, timeCombined_rs_w] = io_CSIload_twix_pair( ...
    paths.metFile, paths.refFile, kFile);

%% === 2. SPATIAL RECONSTRUCTION (DCF + NUFFT, no MFI) ======================
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

%% === 4b. LEFT-SHIFT (removes FID first-point 1st-order phase) =============
fprintf('\n--- 4b. Left-shift ---\n');
ccav   = op_CSIleftshift(ccav);            % uses ccav.pointsToLeftshift
ccav_w = op_CSIleftshift(ccav_w);

%% === 4c. PROF PREPROCESSING: 2 Hz filter -> zero-pad 4x (time domain) =====
% Apodize FIRST so the FID decays, THEN zero-pad 576 -> 2304 -> finer ppm grid
% with NO truncation ripple. First-order phase is applied after the spectral FT.
lb       = 2;             % 2 Hz filter HURT (broadens narrow phantom lines) -> off.
zfTarget = 576*4;% * 4;       % zero-pad 576 -> 2304  (this is the part that helps)
fprintf('\n--- 4c. Time filter (%g Hz) + zero-pad -> %d ---\n', lb, zfTarget);
ccav   = op_CSIFilter(ccav,   lb);         % 2 Hz exponential LB (op_filter analogue)
ccav_w = op_CSIFilter(ccav_w, lb);
ccav   = op_CSItimeZeroFill(ccav,   zfTarget);
ccav_w = op_CSItimeZeroFill(ccav_w, zfTarget);

%% === 5. SPECTRAL FT =======================================================
fprintf('\n--- 5. Spectral FT ---\n');
ftSpec   = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);

% First-order phase: 1.2e-4 s HURT (data already left-shifted + LCModel phases
% it) -> off. Applied only if foPhase ~= 0. (sweep_ge_preproc: every phase made
% it worse; zero-pad-only was best, MAE 0.18.)
foPhase = 0;              % seconds; set nonzero only if a future dataset needs it
if foPhase ~= 0
    ftSpec   = op_CSIaddPhase(ftSpec,   0, foPhase, 4.65);
    ftSpec_w = op_CSIaddPhase(ftSpec_w, 0, foPhase, 4.65);
end

save('ftSpec.mat',   'ftSpec',   '-v7.3');
save('ftSpec_w.mat', 'ftSpec_w', '-v7.3');

%% === 6. LIPID + WATER REMOVAL + B0 CORRECTION =============================
fprintf('\n--- 6. Lipid / water removal + B0 correction ---\n');
% Water removal (op_CSIRemoveLipids basis FIXED to the data's ifft convention;
% before, it removed ~0% of water). Band narrowed to protect Cho/Ins 4.05, Cr 3.91.
ftSpec_rmw = op_CSIRemoveLipids(ftSpec, ...
                'lipidPPMRange',  [4.55 4.95], ...
                'linewidthRange', [1 10]);
[ftSpec_B0corr, ftSpec_B0corr_w, freqMap, R2Map] = ...
    op_CSIB0Correction_v2(ftSpec_rmw, ftSpec_w);

% Spatial smoothing: 0 for this 3-compartment phantom (any blur mixes the
% Small/Large/Outer concentrations). If you do smooth, it must be BOTH scans
% UNMASKED then mask BOTH (mask-first biased metab/water low at the edges).
smoothFwhm = 0;                                    % mm; 0 = none (recommended here)
if smoothFwhm > 0
    ftSpec_sm   = op_CSIApodize(ftSpec_B0corr,  'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
    ftSpec_sm_w = op_CSIApodize(ftSpec_B0corr_w,'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
else
    ftSpec_sm   = ftSpec_B0corr;
    ftSpec_sm_w = ftSpec_B0corr_w;
end
ftSpec_sm_w.mask = ftSpec_B0corr.mask;
ftSpec_smooth    = op_CSIapplymask(ftSpec_sm);
ftSpec_smooth_w  = op_CSIapplymask(ftSpec_sm_w);


op_CSIPlot(ftSpec_smooth);

%% === 7. WRITE 4D NIFTI + VIEWER ===========================================
fprintf('\n--- 7. Writing 4D MRSI NIfTI ---\n');
mrsiout = mrsi_on_t1_map(ccav_w, ftSpec_smooth, paths.mrsiOut, paths.emptyNiiPath);
nii_viewer(paths.t1Path, mrsiout.mrsi4D_time);

lcmbin    = '/home/divya/.lcmodel/bin/lcmodel';
basisFile = 'C:\Users\divya\Downloads\fida codes\fid_a\basis_rosette_phantom\RosettePhantom_LacAceCrCho_TE15_SW1587.basis';

%% === 8. LCMODEL + METABOLITE MAPS =========================================
% LCModel timing derived from the data so the control matches the basis.
specDim = ftSpec_smooth.dims.f;
if specDim == 0, specDim = ftSpec_smooth.dims.t; end
hzpppm = ftSpec_smooth.txfrq / 1e6;
deltat = 1 / ftSpec_smooth.spectralWidth;
nunfil = ftSpec_smooth.sz(specDim);
echot  = ftSpec_smooth.te;
fprintf('  LCModel timing (from data): HZPPPM=%.5f  DELTAT=%.8f  NUNFIL=%d  ECHOT=%g\n', ...
        hzpppm, deltat, nunfil, echot);
if abs(echot - 15) > 1
    warning('Data TE=%.1f ms but basis is TE=15 ms. Rebuild the basis at the data TE.', echot);
end

fprintf('\n--- 8. LCModel + maps ---\n');
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
    'doecc',      false);   % B0 correction already did the water-phase in-pipeline

Nx = numel(getCoordinates(ftSpec, 'x'));
Ny = numel(getCoordinates(ftSpec, 'y'));
[map, crlb, LW, SNR] = op_CSILCModelMaps(Nx, Ny, paths.rawFolder, ...
    'figure_folder_name', 'maps', ...
    'metabList', {'Lac' 'Act' 'Cr' 'Cho'});
save('map.mat',  'map',  '-v7.3');
save('crlb.mat', 'crlb', '-v7.3');
save('LW.mat',   'LW',   '-v7.3');
save('SNR.mat',  'SNR',  '-v7.3');

apply_water_scaling_all;   % 3-compartment GT + raw mM (relaxation off)

out = create_separate_metabolite_niftis_v2(ftSpec_smooth, map, crlb, LW, SNR, ...
        paths.t1Path, paths.emptyNiiPath, paths.mapsDir);
%nii_viewer(paths.t1Path, out.metabolite_4d_files);
allMaps = [out.metabolite_4d_files, out.crlb_3d_files, {out.shared_lw_3d_file, out.shared_snr_3d_file}];
%nii_viewer(paths.t1Path, allMaps);
fprintf('\n=== Pipeline complete.  Outputs in %s ===\n', paths.outputDir);
