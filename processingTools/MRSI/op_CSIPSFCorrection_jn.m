function MRSIStruct = op_CSIPSFCorrection_jn(MRSIStruct, k_space_file, NameValueArgs)
%OP_CSIPSFCORRECTION_JN  Density compensation for non-Cartesian MRSI.
%
%   Methods:
%     'nearest'    – k-nearest-neighbour density estimate
%     'voronoi'    – Voronoi cell area
%     'pipe_menon' – Pipe-Menon iterative (requires Fessler IRT)
%
%   The normalised weights satisfy  sum(w) = Nk  so that the adjoint
%   NUFFT in op_NUFFTSpatial1 (which divides by Nk) recovers correct
%   amplitudes.

arguments
    MRSIStruct   (1,1) struct
    k_space_file (1,:) char {mustBeFile}
    NameValueArgs.method       (1,:) char {mustBeMember(NameValueArgs.method, ...
        {'nearest','voronoi','pipe_menon'})} = 'nearest'
    NameValueArgs.modelType    (1,:) char {mustBeMember(NameValueArgs.modelType, ...
        {'Uniform','Gaussian','FlatEdge'})} = 'Uniform'
    NameValueArgs.sigma        (1,1) double = sqrt(-(0.1^2)/(2*log(0.01)))
    NameValueArgs.steep        (1,1) double = -100
    NameValueArgs.numNeighbors (1,1) double = 10
    NameValueArgs.numIterations(1,1) double = 25
    NameValueArgs.isPlotWeights(1,1) logical = false
end

method = lower(NameValueArgs.method);
fprintf('\n=== DENSITY COMPENSATION (%s) ===\n', method);

% --- Guard: must run before spatial FT ---
if isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'spatialft') ...
        && MRSIStruct.flags.spatialft
    error('Density compensation must be performed prior to spatial FT.');
end

% --- Read k-space trajectory ---
[kx, ky] = readKTrajectory(k_space_file);
Nk = numel(kx);
fprintf('[DEBUG] K-space: %d pts  kx=[%.4f, %.4f]  ky=[%.4f, %.4f]\n', ...
    Nk, min(kx), max(kx), min(ky), max(ky));

% --- Data dimensions ---
[Nt, Ncoils, Nave, Nkpts, Nkshot, hasAvg] = getDataDims(MRSIStruct);
fprintf('[DEBUG] Dims: Nt=%d  Ncoils=%d  Nave=%d  Nkpts=%d  Nkshot=%d\n', ...
    Nt, Ncoils, Nave, Nkpts, Nkshot);
assert(Nk == Nkpts * Nkshot, ...
    'K-file has %d points but data expects %d.', Nk, Nkpts*Nkshot);

% --- Compute raw weights (method-specific) ---
switch method
    case 'nearest'
        w_raw = dcf_nearest(kx, ky, NameValueArgs.numNeighbors);
    case 'voronoi'
        w_raw = dcf_voronoi(kx, ky, MRSIStruct, ...
            NameValueArgs.modelType, NameValueArgs.sigma, NameValueArgs.steep);
    case 'pipe_menon'
        w_raw = dcf_pipe_menon(kx, ky, MRSIStruct, ...
            NameValueArgs.numIterations);
end

% --- Normalise so sum(w) = Nk ---
w = w_raw * (Nk / sum(w_raw));
fprintf('[DEBUG] Normalised weights: min=%.4e  max=%.4e  sum=%.1f  mean=%.4f\n', ...
    min(w), max(w), sum(w), mean(w));

% --- Radial profile diagnostic ---
printRadialProfile(kx, ky, w, method);

% --- Optional weight plot ---
if NameValueArgs.isPlotWeights
    figure;
    scatter(kx, ky, 8, w, 'filled'); axis equal; colorbar;
    title(sprintf('DCF weights (%s)', method)); xlabel('kx'); ylabel('ky');
end

% --- Reshape weights to match data and apply ---
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

% --- Store metadata ---
MRSIStruct.densityComp = struct( ...
    'method', method, 'weights', w(:), 'weights_raw', w_raw(:), ...
    'kx', kx(:), 'ky', ky(:), 'weight_matrix', W_mat);

fprintf('=== DCF COMPLETE ===\n\n');
end


function w = dcf_nearest(kx, ky, K)
    N = numel(kx);
    fprintf('[NN] Computing k-NN density  (N=%d, K=%d) ...\n', N, K);

    rho  = zeros(N, 1);       % geometric density per point
    kRep = zeros(N, 1);       % number of co-located duplicates

    for i = 1:N
        d       = hypot(kx - kx(i), ky - ky(i));
        kRep(i) = sum(d == 0);                        %includes self
        sd      = sort(d);
        kdist   = sd(min(K + kRep(i), N));             %skip duplicates justincase
        rho(i)  = K / (pi * max(kdist, eps)^2);         %same as jamie's/lubna's code
    end

    %geometric weight (ignoring duplicates)
    w_geo = 1 ./ max(rho, eps);
    fprintf('[NN] Geometric weights:  min=%.4e  max=%.4e  ratio=%.1f\n', ...
        min(w_geo), max(w_geo), max(w_geo)/min(w_geo));

    %radial smoothing (geometry only, kRep not yet applied)
    window = max(round(N / 50), 5);
    w_geo  = radialSmooth(kx, ky, w_geo, window); %%this is technically cheating so talk to mark and jamie
    fprintf('[NN] After radial smooth (win=%d):  min=%.4e  max=%.4e\n', ...
        window, min(w_geo), max(w_geo));

    %divide by kRep so each duplicate gets its fair share
    w = w_geo ./ kRep;


    fprintf('[NN] After kRep correction:  min=%.4e  max=%.4e  (max kRep=%d, %d dup pts)\n', ...
        min(w), max(w), max(kRep), sum(kRep > 1));

    %clip the max 2% -cap outliers
    cap = prctile(w, 98);
    w   = min(w, cap);
    fprintf('[NN] After 98th-pctl cap (%.4e):  min=%.4e  max=%.4e\n', cap, min(w), max(w));
end


%% =================================================================
%%  METHOD 2 — VORONOI
%%  Weight = Voronoi cell area ∝ 1/density.
%%
%%  Key design choices:
%%    - Boundary (infinite) cells get p90 of finite areas.
%%    - mapBackShared divides weight by number of duplicates
%%      mapping to the same unique cell.
%%    - Smooth taper beyond Kmax instead of hard zero.
%% =================================================================
function w = dcf_voronoi(kx, ky, S, modelType, sigma, steep)
    fprintf('[VOR] Computing Voronoi areas ...\n');

    % Remove duplicate origin points for Voronoi tessellation
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
        if any(out)      % clip to k-space disk
            vx(out) = Kmax * vx(out) ./ r(out);
            vy(out) = Kmax * vy(out) ./ r(out);
        end
        a = polyarea(vx, vy);
        if a > 0, areas(i) = a; end
    end

    % Assign boundary cells the 90th-percentile finite area
    finMask  = ~isnan(areas) & ~isBndry;
    needArea = isnan(areas) | isBndry;
    bndArea  = prctile(areas(finMask), 90);
    areas(needArea) = bndArea;

    % Safety fill for any remaining NaN/zero
    medArea = median(areas(areas > 0 & ~isnan(areas)));
    areas(isnan(areas) | areas <= 0) = medArea;

    fprintf('[VOR] Finite: %d  Boundary: %d  bndArea=%.4e  medArea=%.4e\n', ...
        sum(finMask), sum(needArea), bndArea, medArea);

    % Density shaping (Uniform / Gaussian / FlatEdge)
    w_c = areas;
    R   = hypot(xc, yc);
    switch lower(modelType)
        case 'uniform'
            far = R > Kmax;
            if any(far)
                taper = max(1 - (R(far) - Kmax) / (0.1 * Kmax), 0);
                w_c(far) = w_c(far) .* taper;
            end
        case 'gaussian'
            w_c = w_c .* exp(-(xc.^2 + yc.^2) / (2 * sigma^2));
        case 'flatedge'
            w_c = w_c ./ (1 + exp(steep * (R - 0.64 * Kmax)));
    end

    % Map unique-point weights back, sharing among duplicates
    w = mapBackShared(kx, ky, xc, yc, w_c);
    fprintf('[VOR] Final weights: min=%.4e  max=%.4e\n', min(w), max(w));
end


%% =================================================================
%%  METHOD 3 — PIPE-MENON ITERATIVE
%%  w converges to 1/diag(E*E^H) via fixed-point iteration using
%%  the full NUFFT forward/adjoint operators (not just st.p).
%%
%%  Key design choices:
%%    - abs(Gw) + relative floor to avoid negative/zero denominators.
%%    - Per-iteration mean-normalisation to prevent magnitude drift.
%%    - Symmetric percentile clipping + radial smoothing post-iteration.
%%    - Falls back to st.p if full NUFFT diverges.
%% =================================================================
function w = dcf_pipe_menon(kx, ky, S, iters)
    fprintf('[PM] Pipe-Menon (%d iterations) ...\n', iters);

    if exist('nufft_init', 'file') ~= 2
        fprintf('[PM] Fessler IRT not found — falling back to NN\n');
        w = dcf_nearest(kx, ky, 10);
        return
    end

    [FOV, Nx] = getGridParams(S);
    dx = FOV / Nx;
    Nk = numel(kx);

    % NUFFT setup — must match op_NUFFTSpatial1
    om     = [2*pi*ky(:)*dx, 2*pi*kx(:)*dx];   % [ky, kx] ordering
    Nd     = [Nx, Nx];
    Jd     = [6, 6];
    Kd     = 2 * Nd;
    nshift = Nd / 2;
    st     = nufft_init(om, Nd, Jd, Kd, nshift);
    fprintf('[PM] NUFFT: Nd=%s  Kd=%s  dx=%.2f mm\n', mat2str(Nd), mat2str(Kd), dx);

    % Fixed-point iteration:  w_{n+1} = w_n / |G * w_n|
    w = ones(Nk, 1);
    converged = true;

    for ii = 1:iters
        Gw  = nufft(nufft_adj(w, st), st);        % Gram matrix * w
        den = max(abs(Gw), 1e-4 * max(abs(Gw)));  % relative floor
        w   = w ./ den;
        w   = w / mean(w);                         % prevent drift

        if any(~isfinite(w))
            fprintf('[PM] *** Diverged at iteration %d ***\n', ii);
            converged = false;
            break
        end
        if ii <= 3 || mod(ii, 5) == 0 || ii == iters
            fprintf('[PM]   iter %2d:  min=%.3e  max=%.3e  ratio=%.1f\n', ...
                ii, min(w), max(w), max(w)/max(min(w), eps));
        end
    end

    % Fallback: sparse interpolation matrix (always stable)
    if ~converged
        fprintf('[PM] Falling back to st.p iteration ...\n');
        w = ones(Nk, 1);
        for ii = 1:iters
            den = real(st.p * (st.p' * w));
            w   = w ./ max(den, 1e-10);
        end
    end
    w = abs(w);

    % Symmetric percentile clipping
    wFloor = prctile(w, 2);
    wCap   = prctile(w, 98);
    w = max(min(w, wCap), wFloor);
    fprintf('[PM] Clipped [p2, p98] = [%.4e, %.4e]:  min=%.4e  max=%.4e\n', ...
        wFloor, wCap, min(w), max(w));

    % Radial smoothing to remove iteration noise
    window = max(round(Nk / 80), 5);
    w = radialSmooth(kx, ky, w, window);
    fprintf('[PM] Radial smooth (win=%d):  min=%.4e  max=%.4e\n', ...
        window, min(w), max(w));
end


%% =================================================================
%%  SHARED UTILITIES
%% =================================================================

function ws = radialSmooth(kx, ky, w, window)
%RADIALSMOOTH  Moving-average of weights sorted by k-space radius.
    [~, si] = sort(hypot(kx, ky));
    tmp     = movmean(w(si), window);
    ws      = w;
    ws(si)  = tmp;
end

function [kx, ky] = readKTrajectory(f)
%READKTRAJECTORY  Read Kx, Ky columns from a k-space CSV/table file.
    T    = readtable(f);
    cols = lower(string(T.Properties.VariableNames));
    ix   = find(cols == "kx", 1);
    iy   = find(cols == "ky", 1);
    assert(~isempty(ix) && ~isempty(iy), ...
        'Kx/Ky columns not found. Available: %s', strjoin(T.Properties.VariableNames, ', '));
    kx = double(T{:, ix});
    ky = double(T{:, iy});
end

function [xc, yc, nRemoved] = removeDuplicateOrigin(kx, ky)
%REMOVEDUPLICATEORIGIN  Keep one copy of (0,0); remove the rest.
    atOrigin = (kx == 0) & (ky == 0);
    nDup     = sum(atOrigin);
    keep     = true(size(kx));
    if nDup > 1
        keep(atOrigin)           = false;
        keep(find(atOrigin, 1))  = true;
    end
    xc       = kx(keep);
    yc       = ky(keep);
    nRemoved = nDup - min(nDup, 1);
end

function w = mapBackShared(kx, ky, xc, yc, wc)
%MAPBACKSHARED  Map unique-point weights to the full k-space list,
%   dividing equally among duplicates that share the same cell.
    N  = numel(kx);
    M  = numel(xc);
    nn = zeros(N, 1);          % nearest-unique index for each point
    for i = 1:N
        [~, nn(i)] = min(hypot(xc - kx(i), yc - ky(i)));
    end
    cnt = accumarray(nn, 1, [M, 1]);
    w   = wc(nn) ./ cnt(nn);  % vectorised — no second loop needed
end

function [Nt, Ncoils, Nave, Nkpts, Nkshot, hasAvg] = getDataDims(S)
%GETDATADIMS  Extract dimension sizes from an MRSI struct.
    Nt     = S.sz(S.dims.t);
    Ncoils = S.sz(S.dims.coils);
    Nkpts  = S.sz(S.dims.kpts);
    Nkshot = S.sz(S.dims.kshot);
    hasAvg = isfield(S.dims, 'averages') && S.dims.averages > 0;
    Nave   = 1;
    if hasAvg, Nave = S.sz(S.dims.averages); end
end

function [FOV, Nx, Kmax] = getGridParams(S)
%GETGRIDPARAMS  FOV (mm), matrix size, and Kmax (1/mm).
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
    Kmax = Nx / (2 * FOV);    % = 1 / (2 * dx)
end

function printRadialProfile(kx, ky, w, method)
%PRINTRADIALPROFILE  10-bin radial weight summary + origin analysis.
    R    = hypot(kx, ky);
    Rmax = max(R);
    edges = linspace(0, Rmax, 11);

    fprintf('[DEBUG] --- Radial profile (%s) ---\n', method);
    fprintf('[DEBUG]   %8s %6s %8s %8s %8s\n', ...
        'R_mid', 'Npts', 'w_mean', 'w_med', 'w_std');
    for b = 1:10
        m = R >= edges(b) & R < edges(b+1);
        if any(m)
            fprintf('[DEBUG]   %8.4f %6d %8.4f %8.4f %8.4f\n', ...
                (edges(b)+edges(b+1))/2, sum(m), mean(w(m)), median(w(m)), std(w(m)));
        end
    end

    orig = R < 1e-6;
    fprintf('[DEBUG] Origin (R~0): N=%d  w_each=%.6f  w_total=%.4f\n', ...
        sum(orig), mean(w(orig)), sum(w(orig)));
    p = prctile(w, [1 5 25 50 75 95 99]);
    fprintf('[DEBUG] Percentiles [1 5 25 50 75 95 99]: %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n', p);
end
