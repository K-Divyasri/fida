%% run_Rosette_40x40_Cr.m
%  Rosette MRSI pipeline with CREATINE referencing (not water scaling).
%  Outputs metabolite maps in mM anchored to Cr, segments the inner bottle by
%  lactate, and prints measured vs recipe ground truth.
%
%  RUN_FULL_PIPELINE = true  -> reconstruct + LCModel-fit from .dat (slow).
%                    = false -> load the already-computed map.mat and just do
%                               Cr referencing + segmentation + stats (fast).


addpath(genpath('C:\Users\divya\Downloads\fida codes\fida27May\fida-main'));
addpath('C:\Users\divya\Downloads\fida codes\fid_a');

%% ================= USER INPUTS =================
datasetFolder = 'F:\fida\divya\28thJULYpHANTOM40X40\subject01';
kFile         = 'C:\Users\divya\Downloads\fida codes\fid_a\processingTools\MRSI\kFiles\Rosette_traj_40x40.txt';
basisFile     = 'C:\Users\divya\Downloads\fida codes\fid_a\basis_rosette_phantom\RosettePhantom_LacAceCrCho_TE15_SW1587.basis';
lcmbin        = '/home/divya/.lcmodel/bin/lcmodel';
lcmOwner      = 'Divyasri Krishnkaumar, Sunnybrook Research Institute';
smoothFwhm    = 20;          % spatial Gaussian FWHM (mm); 0 = none
dcfMethod     = 'pipe_menon';
ftMethod      = 'nufft';

CrTrue        = 10;          % mM -- known creatine (both bottles, recipe)
innerLacVal   = 26;          % measured inner lactate (mM) for the seg threshold

% recipe ground truth (mM) -- phantom2
truth.inner = struct('Lac',8,'Act',6.4,'Cr',10,'Cho',9);
truth.outer = struct('Lac',0,'Act',32, 'Cr',10,'Cho',3);

RUN_FULL_PIPELINE = false;   % false = load existing map.mat (fast)
%% ==============================================

paths  = file_path_load(datasetFolder);
lcmDir = paths.rawFolder;

%% ---- 1-9. FULL PIPELINE (recon -> LCModel -> maps) --- slow, optional ----
if RUN_FULL_PIPELINE
    fprintf('--- 1. Load TWIX pair ---\n');
    [tc, tc_w] = io_CSIload_twix_pair(paths.metFile, paths.refFile, kFile);

    fprintf('--- 2. Spatial recon ---\n');
    ft   = op_CSIRecon(tc,   kFile, dcfMethod, ftMethod);
    ft_w = op_CSIRecon(tc_w, kFile, dcfMethod, ftMethod);

    fprintf('--- 3. Coil combine ---\n');
    [cc_w, phase, weights] = op_CSICombineCoils1(ft_w);
    cc                     = op_CSICombineCoils1(ft, 1, phase, weights);

    fprintf('--- 4. Average + mask ---\n');
    ccav   = op_CSIAverage(cc);
    ccav_w = op_CSIAverage(cc_w);
    ccav_w = op_CSISegment_simple(ccav_w);
    ccav.mask = ccav_w.mask;
    mask   = ccav_w.mask.brainmasks;

    fprintf('--- 4b. Left-shift (remove FID first-point 1st-order phase) ---\n');
    ccav   = op_CSIleftshift(ccav);      % uses ccav.pointsToLeftshift
    ccav_w = op_CSIleftshift(ccav_w);

    fprintf('--- 5. Spectral FT ---\n');
    ftSpec   = op_CSIFourierTransform(ccav);
    ftSpec_w = op_CSIFourierTransform(ccav_w);

    fprintf('--- 6. Lipid/water removal + B0 ---\n');
    ftSpec_rmw = op_CSIRemoveLipids(ftSpec, 'lipidPPMRange',[4.5 5.0], 'linewidthRange',[1 10]);
    [ftSpec_B0, ftSpec_B0_w] = op_CSIB0Correction_v2(ftSpec_rmw, ftSpec_w);
    ftSpec_masked = op_CSIapplymask(ftSpec_B0);

    fprintf('--- 7. Smoothing ---\n');
    if smoothFwhm > 0
        ftSpec_smooth   = op_CSIApodize(ftSpec_masked, 'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
        ftSpec_smooth_w = op_CSIApodize(ftSpec_B0_w,   'functionType','gaussian','fullWidthHalfMax',smoothFwhm);
    else
        ftSpec_smooth = ftSpec_masked;  ftSpec_smooth_w = ftSpec_B0_w;
    end

    fprintf('--- 8. Overlay NIfTI + viewer ---\n');
    mrsiout = mrsi_on_t1_map(ccav_w, ftSpec_smooth, paths.mrsiOut, paths.emptyNiiPath);
    nii_viewer(paths.t1Path, mrsiout.mrsi4D_time);

    fprintf('--- 9. LCModel + maps ---\n');
    hzpppm = ftSpec_smooth.txfrq/1e6;  deltat = 1/ftSpec_smooth.spectralWidth;
    sd = ftSpec_smooth.dims.f; if sd==0, sd = ftSpec_smooth.dims.t; end
    nunfil = ftSpec_smooth.sz(sd);  echot = ftSpec_smooth.te;
    run_lcm_rosette_portable(paths.metFile, ftSpec_smooth, ftSpec_smooth_w, mask, ...
        'basisFile',basisFile, 'lcmodelBin',lcmbin, 'ownerStr',lcmOwner, ...
        'rawFolder',lcmDir, 'hzpppm',hzpppm, 'deltat',deltat, ...
        'nunfil',nunfil, 'echot',echot, 'ppmst',4.25, 'ppmend',0.2, ...
        'doecc', true);   % eddy-current correction ON (homogenises Cr)
    Nx = numel(getCoordinates(ftSpec,'x'));  Ny = numel(getCoordinates(ftSpec,'y'));
    [map, crlb, LW, SNR] = op_CSILCModelMaps(Nx, Ny, lcmDir, 'figure_folder_name','maps');
else
    fprintf('Loading existing LCModel maps from:\n  %s\n', lcmDir);
    map  = load(fullfile(lcmDir,'map.mat')).map;
    crlb = load(fullfile(lcmDir,'crlb.mat')).crlb;
    LW   = load(fullfile(lcmDir,'LW.mat')).LW;
    SNR  = load(fullfile(lcmDir,'SNR.mat')).SNR;
end

%% ---- 10. SEGMENT inner bottle by lactate CRLB (CRLB-driven, scale-free) ---
% Inner bottle = voxels with a RELIABLE lactate fit (Lac CRLB < 10%). CRLB is a
% percentage, so this is independent of Cr referencing -> safe on raw or scaled
% map. Outer voxels have no lactate (CRLB ~999%) and drop out.
fprintf('\n--- 10. Segment inner/outer (lactate CRLB < 10%%) ---\n');
good        = isfinite(map.Cr) & (map.Cr > 0) & (crlb.Cr > 0) & (crlb.Cr <= 20);
phantomMask = good;                            % phantom = voxels with a good Cr fit
roi = segment_inner_phantom(map, crlb, phantomMask, 10);

%% ---- 11. CREATINE REFERENCING (per voxel: every voxel's Cr = CrTrue) -----
fprintf('\n--- 11. Creatine referencing (Cr = %g mM) ---\n', CrTrue);
scale = nan(size(map.Cr));
scale(good) = CrTrue ./ map.Cr(good);          % per-voxel Cr ratio
fn = fieldnames(map);
for i = 1:numel(fn), map.(fn{i}) = map.(fn{i}) .* scale; end
fprintf('  %d voxels Cr-referenced.\n', nnz(good));

%% ---- 12. STATS: measured vs recipe --------------------------------------
fprintf('\n--- 12. Measured vs recipe ground truth ---\n');
mets = {'Lac','Act','Cr','Cho'};
T = table();
for r = {'inner', roi.innerMask; 'outer', roi.outerMask}'
    rn = r{1}; rm = r{2};
    for i = 1:numel(mets)
        m = mets{i};
        v = map.(m)(rm & isfinite(map.(m)) & map.(m) >= 0);
        T = [T; table(string(rn), string(m), numel(v), median(v), mean(v), std(v), ...
                truth.(rn).(m), 'VariableNames', ...
                {'ROI','Metab','N','Median','Mean','SD','Recipe'})]; %#ok<AGROW>
    end
end
disp(T);

% histograms with ground-truth lines
plot_roi_hist(map, crlb, roi);

fprintf('\n=== Done. map is Cr-referenced; compare Median vs Recipe above. ===\n');
