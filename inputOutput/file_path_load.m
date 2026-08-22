function paths = file_path_load(datasetFolder)
%FILE_PATH_LOAD  Resolve every input/output path for an MRSI dataset folder.
%
%   paths = file_path_load(datasetFolder)
%
%   Expected folder layout — anything missing is created or generated:
%       <datasetFolder>/
%           met/         single Siemens TWIX .dat   (water-suppressed)
%           mrs_ref/     single Siemens TWIX .dat   (water reference)
%           img_ref/     T1 NIfTI  (auto-generated from /dicom/ if empty)
%           dicom/       raw DICOM series           (only needed if /img_ref/ empty)
%           outputs/     created if missing
%               lcm_out/  created
%               maps/     created
%
%   The empty NIfTI template (used by mrsi_on_t1_map) is produced by
%   running twix2nifti_v2 on the water-reference .dat if it isn't already
%   present in /outputs/.
%
%   Returns a struct with fields:
%       datasetFolder, metFile, refFile, t1Path,
%       emptyNiiPath, mrsiOut,
%       outputDir, rawFolder, mapsDir.

arguments
    datasetFolder (1,:) char {mustBeFolder}
end

% ---- Subfolder layout --------------------------------------------------
metDir    = fullfile(datasetFolder, 'met');
mrsRefDir = fullfile(datasetFolder, 'mrs_ref');
imgRefDir = fullfile(datasetFolder, 'img_ref');
dicomDir  = fullfile(datasetFolder, 'dicom');
outputDir = fullfile(datasetFolder, 'outputs');
rawFolder = fullfile(outputDir,    'lcm_out');
mapsDir   = fullfile(outputDir,    'maps');

% Create write-side dirs if missing
for d = {imgRefDir, outputDir, rawFolder, mapsDir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

% ---- Single .dat in /met/ and /mrs_ref/ --------------------------------
metFile = pickSingleDat(metDir,    'metabolite');
refFile = pickSingleDat(mrsRefDir, 'water-reference');

% ---- T1 NIfTI in /img_ref/ ---------------------------------------------
t1Path = pickT1(imgRefDir);
if isempty(t1Path)
    fprintf('[file_path_load] /img_ref/ is empty — running dicm2nii on /dicom/...\n');
    assert(exist(dicomDir,'dir')==7, ...
        'No /img_ref/ NIfTI and no /dicom/ folder found in %s', datasetFolder);
    dicm2nii(dicomDir, imgRefDir, 'nii.gz');
    t1Path = pickT1(imgRefDir);
    assert(~isempty(t1Path), ...
        'dicm2nii produced no T1-like NIfTI in %s', imgRefDir);
end

% ---- Empty NIfTI template in /outputs/ ---------------------------------
[~, refBase, ~] = fileparts(refFile);
emptyNiiPath = fullfile(outputDir, [refBase '_empty.nii.gz']);
if ~exist(emptyNiiPath, 'file')
    altNii = fullfile(outputDir, [refBase '_empty.nii']);
    if exist(altNii, 'file')
        gzip(altNii); delete(altNii);
    else
        fprintf('[file_path_load] Empty NIfTI not found — running twix2nifti_v2...\n');
        emptyOutBase = fullfile(outputDir, [refBase '_empty']);   % v2 appends .nii
        twix2nifti_v2(refFile, emptyOutBase);
        gzip([emptyOutBase '.nii']);
        delete([emptyOutBase '.nii']);
    end
end
assert(exist(emptyNiiPath,'file')==2, ...
    'Failed to find or generate empty NIfTI: %s', emptyNiiPath);

% ---- Filled MRSI 4D NIfTI output path ----------------------------------
mrsiOut = fullfile(outputDir, [refBase '_filled.nii.gz']);

% ---- Pack --------------------------------------------------------------
paths = struct( ...
    'datasetFolder', datasetFolder, ...
    'metFile',       metFile, ...
    'refFile',       refFile, ...
    't1Path',        t1Path, ...
    'emptyNiiPath',  emptyNiiPath, ...
    'mrsiOut',       mrsiOut, ...
    'outputDir',     outputDir, ...
    'rawFolder',     rawFolder, ...
    'mapsDir',       mapsDir);

fprintf('\n=== file_path_load: resolved paths ===\n');
fprintf('  met          : %s\n', paths.metFile);
fprintf('  ref          : %s\n', paths.refFile);
fprintf('  t1           : %s\n', paths.t1Path);
fprintf('  emptyNii     : %s\n', paths.emptyNiiPath);
fprintf('  mrsiOut      : %s\n', paths.mrsiOut);
fprintf('  outputDir    : %s\n', paths.outputDir);
fprintf('  rawFolder    : %s\n', paths.rawFolder);
fprintf('  mapsDir      : %s\n', paths.mapsDir);
fprintf('======================================\n\n');
end


% ========================================================================
function f = pickSingleDat(dirPath, label)
    assert(exist(dirPath,'dir')==7, '%s folder missing: %s', label, dirPath);
    d = dir(fullfile(dirPath, '*.dat'));
    assert(~isempty(d), 'No .dat in %s/  — expected one %s file', dirPath, label);
    if numel(d) > 1
        warning('file_path_load:multipleDat', ...
            '%d .dat files in %s/ — using the first: %s', numel(d), dirPath, d(1).name);
    end
    f = fullfile(d(1).folder, d(1).name);
end


function t1 = pickT1(imgRefDir)
    if ~exist(imgRefDir,'dir'), t1 = ''; return; end
    nii = [dir(fullfile(imgRefDir,'*.nii.gz')); dir(fullfile(imgRefDir,'*.nii'))];
    if isempty(nii), t1 = ''; return; end

    names = {nii.name};
    % Deterministic ranking. Alphabetical order picked whatever sorted first,
    % which differs between exports (a bare t1_..._0.nii.gz sorts before the
    % _MPR_cor reformat because '.' < '_') and can flip/shift the overlay.
    % Rank instead: mprage/t1 keyword first, prefer a coronal reformat, and
    % PENALISE the non-distortion-corrected (_ND) copies.
    score = zeros(1, numel(names));
    for i = 1:numel(names)
        n = lower(names{i}); s = 0;
        if contains(n,'mprage')||contains(n,'t1w')||contains(n,'t1_'), s = s + 100; end
        if contains(n,'mpr_cor'), s = s + 40;                % same reformat as before
        elseif contains(n,'mpr_tra'), s = s + 30;
        elseif contains(n,'mpr_sag'), s = s + 20; end
        if contains(n,'_nd'), s = s - 50; end                % avoid _ND variants
        score(i) = s;
    end
    [~, order] = sort(score, 'descend');                     % stable: ties keep dir() order
    chosen = nii(order(1));
    t1 = fullfile(chosen.folder, chosen.name);

    % Some exports write a 4D anatomical (e.g. 2 volumes). nii_viewer/overlay
    % expect 3D; take volume 1 and write a 3D sidecar so geometry is unambiguous.
    t1 = ensure3D(t1);
end

function p = ensure3D(p)
    try
        h = nii_tool('hdr', p);
        nvol = double(h.dim(5));
        if h.dim(1) >= 4 && nvol > 1
            [d,n,e] = fileparts(p);
            if strcmpi(e,'.gz'), [~,n2,e2] = fileparts(n); n = n2; e = [e2 e]; end
            out = fullfile(d, [n '_vol1' e]);
            if exist(out,'file') ~= 2
                v = nii_tool('load', p);
                v.img = v.img(:,:,:,1);
                v.hdr.dim(1) = 3; v.hdr.dim(5) = 1;
                nii_tool('save', v, out);
                fprintf('[pickT1] 4D T1 (%d vols) -> wrote 3D volume-1: %s\n', nvol, out);
            end
            p = out;
        end
    catch ME
        warning('pickT1:ensure3D','Could not 3D-ify %s (%s) -- using as-is.', p, ME.message);
    end
end
