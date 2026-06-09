function MRSIStruct = op_CSIPSFCorrection_v(MRSIStruct, k_space_file, NameValueArgs)
%OP_CSIPSFCORRECTION_V  Density compensation — Voronoi cell area method.
%
%   Weight = Voronoi cell area proportional to 1/density.
%
%   Key design choices:
%     - Boundary (infinite) cells get p90 of finite areas.
%     - mapBackShared divides weight equally among duplicates
%       mapping to the same unique cell.
%     - Smooth taper beyond Kmax instead of hard zero.
%     - sum(w) = Nk so op_NUFFTSpatial1 recovers correct amplitudes.
%
%   USAGE:
%     dComp = op_CSIPSFCorrection_v(MRSIStruct, kFile)
%     dComp = op_CSIPSFCorrection_v(MRSIStruct, kFile, 'modelType', 'Gaussian')
%     dComp = op_CSIPSFCorrection_v(MRSIStruct, kFile, 'isPlotWeights', true)

arguments
    MRSIStruct    (1,1) struct
    k_space_file  (1,:) char {mustBeFile}
    NameValueArgs.modelType     (1,:) char {mustBeMember(NameValueArgs.modelType, ...
        {'Uniform','Gaussian','FlatEdge'})} = 'Uniform'
    NameValueArgs.sigma         (1,1) double  = sqrt(-(0.1^2)/(2*log(0.01)))
    NameValueArgs.steep         (1,1) double  = -100
    NameValueArgs.isPlotWeights (1,1) logical = false
end

fprintf('\n=== DENSITY COMPENSATION (voronoi) ===\n');

%% --- Guard ---
if isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'spatialft') ...
        && MRSIStruct.flags.spatialft
    error('Density compensation must be performed prior to spatial FT.');
end

%% --- Read k-space trajectory ---
[kx, ky] = readKTrajectory(k_space_file);
Nk = numel(kx);
fprintf('[DEBUG] K-space: %d pts  kx=[%.4f, %.4f]  ky=[%.4f, %.4f]\n', ...
    Nk, min(kx), max(kx), min(ky), max(ky));

%% --- Data dimensions ---
[Nt, Ncoils, Nave, Nkpts, Nkshot, hasAvg] = getDataDims(MRSIStruct);
fprintf('[DEBUG] Dims: Nt=%d  Ncoils=%d  Nave=%d  Nkpts=%d  Nkshot=%d\n', ...
    Nt, Ncoils, Nave, Nkpts, Nkshot);
assert(Nk == Nkpts * Nkshot, ...
    'K-file has %d points but data expects %d.', Nk, Nkpts*Nkshot);

%% --- Compute raw weights ---
w_raw = dcf_voronoi(kx, ky, MRSIStruct, ...
    NameValueArgs.modelType, NameValueArgs.sigma, NameValueArgs.steep);

%% --- Normalise so sum(w) = Nk ---
w = w_raw * (Nk / sum(w_raw));
fprintf('[DEBUG] Normalised weights: min=%.4e  max=%.4e  sum=%.1f  mean=%.4f\n', ...
    min(w), max(w), sum(w), mean(w));

%% --- Radial profile diagnostic ---
printRadialProfile(kx, ky, w);

%% --- Optional weight plot ---
if NameValueArgs.isPlotWeights
    figure;
    scatter(kx, ky, 8, w, 'filled'); axis equal; colorbar;
    title(sprintf('DCF weights (voronoi, %s)', NameValueArgs.modelType));
    xlabel('kx'); ylabel('ky');
end

%% --- Reshape and apply ---
W_mat = reshape(w, [Nkpts, Nkshot]);
if hasAvg && ndims(MRSIStruct.data) == 5
    fullW = repmat(reshape(W_mat, [1 1 1 Nkpts Nkshot]), [Nt Ncoils Nave 1 1]);
else
    fullW = repmat(reshape(W_mat, [1 1 Nkpts Nkshot]),   [Nt Ncoils 1 1]);
end
fullW = cast(fullW, 'like', MRSIStruct.data);
assert(isequal(size(MRSIStruct.data), size(fullW)), ...
    'Size mismatch: data=%s  weights=%s', ...
    mat2str(size(MRSIStruct.data)), mat2str(size(fullW)));

MRSIStruct.data = MRSIStruct.data .* fullW;

%% --- Store metadata ---
MRSIStruct.densityComp = struct( ...
    'method',       'voronoi', ...
    'modelType',    NameValueArgs.modelType, ...
    'weights',      w(:), ...
    'weights_raw',  w_raw(:), ...
    'kx',           kx(:), ...
    'ky',           ky(:), ...
    'weight_matrix',W_mat);

fprintf('=== DCF COMPLETE (voronoi) ===\n\n');
end


%% =================================================================
%%  CORE DENSITY FUNCTION
%% =================================================================
function w = dcf_voronoi(kx, ky, S, modelType, sigma, steep)
    fprintf('[VOR] Computing Voronoi areas ...\n');

    % Remove duplicate origin points
    [xc, yc, nRem] = removeDuplicateOrigin(kx, ky);
    fprintf('[VOR] Unique pts: %d  (removed %d origin duplicates)\n', numel(xc), nRem);

    [v, c] = voronoin([xc, yc]);

    [FOV, Nx, Kmax] = getGridParams(S);
    fprintf('[VOR] FOV=%.1f mm  Nx=%d  Kmax=%.4f 1/mm\n', FOV, Nx, Kmax);

    % Compute clipped cell areas
    nCells  = numel(c);
    areas   = nan(nCells, 1);
    isBndry = false(nCells, 1);

    for i = 1:nCells
        idx = c{i};
        if any(idx == 1) || numel(idx) < 3
            isBndry(i) = true;
            continue
        end
        vx = v(idx, 1);  vy = v(idx, 2);
        r  = hypot(vx, vy);
        out = r > Kmax;
        if any(out)
            vx(out) = Kmax * vx(out) ./ r(out);
            vy(out) = Kmax * vy(out) ./ r(out);
        end
        a = polyarea(vx, vy);
        if a > 0, areas(i) = a; end
    end

    % Boundary cells get p90 of finite areas
    finMask  = ~isnan(areas) & ~isBndry;
    needArea = isnan(areas) | isBndry;
    bndArea  = prctile(areas(finMask), 90);
    areas(needArea) = bndArea;

    % Safety fill
    medArea = median(areas(areas > 0 & ~isnan(areas)));
    areas(isnan(areas) | areas <= 0) = medArea;

    fprintf('[VOR] Finite: %d  Boundary: %d  bndArea=%.4e  medArea=%.4e\n', ...
        sum(finMask), sum(needArea), bndArea, medArea);

    % Density shaping
    w_c = areas;
    R   = hypot(xc, yc);
    switch lower(modelType)
        case 'uniform'
            far = R > Kmax;
            if any(far)
                taper    = max(1 - (R(far) - Kmax) / (0.1 * Kmax), 0);
                w_c(far) = w_c(far) .* taper;
            end
        case 'gaussian'
            w_c = w_c .* exp(-(xc.^2 + yc.^2) / (2 * sigma^2));
        case 'flatedge'
            w_c = w_c ./ (1 + exp(steep * (R - 0.64 * Kmax)));
    end

    % Map back with sharing among duplicates
    w = mapBackShared(kx, ky, xc, yc, w_c);
    fprintf('[VOR] Final weights: min=%.4e  max=%.4e\n', min(w), max(w));
end


%% =================================================================
%%  SHARED UTILITIES
%% =================================================================
function [kx, ky] = readKTrajectory(f)
    T    = readtable(f);
    cols = lower(string(T.Properties.VariableNames));
    ix   = find(cols == "kx", 1);
    iy   = find(cols == "ky", 1);
    assert(~isempty(ix) && ~isempty(iy), ...
        'Kx/Ky columns not found. Available: %s', ...
        strjoin(T.Properties.VariableNames, ', '));
    kx = double(T{:, ix});
    ky = double(T{:, iy});
end

function [xc, yc, nRemoved] = removeDuplicateOrigin(kx, ky)
    atOrigin = (kx == 0) & (ky == 0);
    nDup     = sum(atOrigin);
    keep     = true(size(kx));
    if nDup > 1
        keep(atOrigin)          = false;
        keep(find(atOrigin, 1)) = true;
    end
    xc       = kx(keep);
    yc       = ky(keep);
    nRemoved = nDup - min(nDup, 1);
end

function w = mapBackShared(kx, ky, xc, yc, wc)
    N  = numel(kx);
    M  = numel(xc);
    nn = zeros(N, 1);
    for i = 1:N
        [~, nn(i)] = min(hypot(xc - kx(i), yc - ky(i)));
    end
    cnt = accumarray(nn, 1, [M, 1]);
    w   = wc(nn) ./ cnt(nn);
end

function [Nt, Ncoils, Nave, Nkpts, Nkshot, hasAvg] = getDataDims(S)
    Nt     = S.sz(S.dims.t);
    Ncoils = S.sz(S.dims.coils);
    Nkpts  = S.sz(S.dims.kpts);
    Nkshot = S.sz(S.dims.kshot);
    hasAvg = isfield(S.dims, 'averages') && S.dims.averages > 0;
    Nave   = 1;
    if hasAvg, Nave = S.sz(S.dims.averages); end
end

function [FOV, Nx, Kmax] = getGridParams(S)
    if isfield(S,'fov') && isfield(S.fov,'x') && S.fov.x > 0
        FOV = S.fov.x;
    else
        FOV = 240;
        fprintf('[WARN] FOV not found, defaulting to 240 mm\n');
    end
    if isfield(S,'voxelSize') && isfield(S.voxelSize,'x') && S.voxelSize.x > 0
        Nx = round(FOV / S.voxelSize.x);
    else
        Nx = 48;
        fprintf('[WARN] voxelSize not found, defaulting to Nx=48\n');
    end
    Kmax = Nx / (2 * FOV);
end

function printRadialProfile(kx, ky, w)
    R     = hypot(kx, ky);
    Rmax  = max(R);
    edges = linspace(0, Rmax, 11);
    fprintf('[DEBUG] --- Radial profile (voronoi) ---\n');
    fprintf('[DEBUG]   %8s %6s %8s %8s %8s\n', ...
        'R_mid','Npts','w_mean','w_med','w_std');
    for b = 1:10
        m = R >= edges(b) & R < edges(b+1);
        if any(m)
            fprintf('[DEBUG]   %8.4f %6d %8.4f %8.4f %8.4f\n', ...
                (edges(b)+edges(b+1))/2, sum(m), ...
                mean(w(m)), median(w(m)), std(w(m)));
        end
    end
    orig = R < 1e-6;
    fprintf('[DEBUG] Origin: N=%d  w_each=%.6f  w_total=%.4f\n', ...
        sum(orig), mean(w(orig)), sum(w(orig)));
    p = prctile(w, [1 5 25 50 75 95 99]);
    fprintf('[DEBUG] Percentiles [1 5 25 50 75 95 99]: %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n', p);
end