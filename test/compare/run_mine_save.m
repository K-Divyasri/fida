%% run_mine_save.m — USER's pipeline + intermediate saves for comparison.
% Runs the user's recon on the same uline phantom and saves the key stages to
%   F:\fida\divya\20260831_ulinePhantom_test\subject04\compare\mine\
% RUN THIS WITH THE USER'S TOOLBOX ON THE PATH (fida27May\fida-main + fid_a).
% Run rosette4divya_save.m separately (Jason's toolbox), then compare_recons.m.

datasetFolder = 'F:\fida\divya\20260831_ulinePhantom_test\subject04';
kFile         = 'C:\Users\divya\Downloads\fida codes\fid_a\processingTools\MRSI\kFiles\Rosette_traj_40x40.txt';
dcfMethod = 'nn';   ftMethod = 'dft';

outDir = 'F:\fida\divya\20260831_ulinePhantom_test\subject04\compare\mine';
if ~exist(outDir,'dir'), mkdir(outDir); end

paths = file_path_load(datasetFolder);

%% RECON (user's chain) ====================================================
[timeCombined_rs, timeCombined_rs_w] = io_CSIload_twix_pair(paths.metFile, paths.refFile, kFile);

ft   = op_CSIRecon(timeCombined_rs,   kFile, dcfMethod, ftMethod);
ft_w = op_CSIRecon(timeCombined_rs_w, kFile, dcfMethod, ftMethod);

[coilCombined_w, phase, weights] = op_CSICombineCoils1(ft_w);
coilCombined                     = op_CSICombineCoils1(ft, 1, phase, weights);

ccav   = op_CSIAverage(coilCombined);
ccav_w = op_CSIAverage(coilCombined_w);
ccav_w = op_CSISegment_simple(ccav_w);
ccav.mask = ccav_w.mask;

% save the RECON OUTPUT here (before left-shift) so it lines up with Jason's ccav
snap(outDir,'ccav',   ccav);
snap(outDir,'ccav_w', ccav_w);

%% downstream (user's) so the final spectrum is comparable too =============
ccav   = op_CSIleftshift(ccav);
ccav_w = op_CSIleftshift(ccav_w);
ftSpec   = op_CSIFourierTransform(ccav);
ftSpec_w = op_CSIFourierTransform(ccav_w);
ftSpec_rmw = op_CSIRemoveLipids(ftSpec, 'lipidPPMRange',[4.5 5.0], 'linewidthRange',[1 10]);
[ftSpec_B0corr, ftSpec_B0corr_w, freqMap, R2Map] = op_CSIB0Correction_v2(ftSpec_rmw, ftSpec_w);
ftSpec_smooth = op_CSIapplymask(ftSpec_B0corr);      % user's final (no smoothing)

snap(outDir,'ftSpec',        ftSpec);
snap(outDir,'ftSpec_smooth', ftSpec_smooth);
fprintf('\nMINE intermediates saved to %s\n', outDir);

%% ---- local: save reduced struct (data+dims+sz+ppm) ----------------------
function snap(outDir, name, S)
    snapStruct = struct();
    snapStruct.data = S.data;
    snapStruct.dims = S.dims;
    snapStruct.sz   = S.sz;
    if isfield(S,'ppm'), snapStruct.ppm = S.ppm; else, snapStruct.ppm = []; end
    save(fullfile(outDir,[name '.mat']), 'snapStruct', '-v7.3');
    fprintf('  saved %s  (data %s)\n', name, mat2str(size(S.data)));
end
