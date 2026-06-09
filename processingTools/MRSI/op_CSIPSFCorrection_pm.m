function MRSIStruct = op_CSIPSFCorrection_pm(MRSIStruct, k_space_file, NameValueArgs)
%OP_CSIPSFCORRECTION_PM  Density compensation — Pipe-Menon iterative method.
%
%   USAGE:
%     dComp = op_CSIPSFCorrection_pm(MRSIStruct, kFile)
%     dComp = op_CSIPSFCorrection_pm(MRSIStruct, kFile, 'numIterations', 30)

arguments
    MRSIStruct    (1,1) struct
    k_space_file  (1,:) char {mustBeFile}
    NameValueArgs.numIterations (1,1) double  = 25
    NameValueArgs.isPlotWeights (1,1) logical = false
end

fprintf('\n=== DENSITY COMPENSATION (pipe-menon) ===\n');

%% Guard
if isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'spatialft') ...
        && MRSIStruct.flags.spatialft
    error('Density compensation must be performed prior to spatial FT.');
end

%% Read k-space
[kx, ky] = readKTrajectory(k_space_file);
Nk = numel(kx);
fprintf('K-space: %d points\n', Nk);

%% Data dimensions
[Nt, Ncoils, Nave, Nkpts, Nkshot, hasAvg] = getDataDims(MRSIStruct);
fprintf('Dims: Nt=%d  Ncoils=%d  Nave=%d  Nkpts=%d  Nkshot=%d\n', ...
    Nt, Ncoils, Nave, Nkpts, Nkshot);
assert(Nk == Nkpts * Nkshot, ...
    'K-file has %d points but data expects %d.', Nk, Nkpts*Nkshot);

%% Compute weights
w_raw = pipe_menon(kx, ky, MRSIStruct, NameValueArgs.numIterations);

%% Normalise so sum(w) = Nk
w = w_raw * (Nk / sum(w_raw));
fprintf('Weights: min=%.4e  max=%.4e  mean=%.4f\n', min(w), max(w), mean(w));

%% Optional plot
if NameValueArgs.isPlotWeights
    figure;
    scatter(kx, ky, 8, w, 'filled');
    axis equal; colorbar;
    title('Pipe-Menon DCF weights');
    xlabel('kx'); ylabel('ky');
end

%% Reshape and apply
W_mat = reshape(w, [Nkpts, Nkshot]);
if hasAvg && ndims(MRSIStruct.data) == 5
    fullW = repmat(reshape(W_mat, [1 1 1 Nkpts Nkshot]), [Nt Ncoils Nave 1 1]);
else
    fullW = repmat(reshape(W_mat, [1 1 Nkpts Nkshot]),   [Nt Ncoils 1 1]);
end
fullW = cast(fullW, 'like', MRSIStruct.data);
MRSIStruct.data = MRSIStruct.data .* fullW;

%% Store metadata
MRSIStruct.densityComp = struct( ...
    'method',        'pipe_menon', ...
    'numIterations', NameValueArgs.numIterations, ...
    'weights',       w(:), ...
    'kx',            kx(:), ...
    'ky',            ky(:));

fprintf('=== DCF COMPLETE (pipe-menon) ===\n\n');
end


%% =================================================================
%%  PIPE-MENON CORE
%% =================================================================
function w = pipe_menon(kx, ky, S, iters)

    %% Check for Fessler IRT
    if exist('nufft_init', 'file') ~= 2
        fprintf('Fessler IRT not found — using nearest-neighbour fallback\n');
        w = nn_fallback(kx, ky);
        return
    end

    Nk  = numel(kx);
    FOV = S.fov.x;
    Nx  = round(FOV / S.voxelSize.x);
    dx  = FOV / Nx;

    %% NUFFT setup — must match op_NUFFTSpatial1
    om     = [2*pi*ky(:)*dx,  2*pi*kx(:)*dx];  % [ky, kx] ordering
    Nd     = [Nx, Nx];
    Jd     = [6,  6];
    Kd     = 2 * Nd;
    nshift = Nd / 2;
    st     = nufft_init(om, Nd, Jd, Kd, nshift);
    fprintf('Pipe-Menon: %d iters  Nx=%d  dx=%.1f mm\n', iters, Nx, dx);

    %% Fixed-point iteration
    % w converges to 1/diag(E*E^H) — the optimal DCF weights
    w = ones(Nk, 1);

    for ii = 1:iters
        % G*w = E * (E^H * w) — Gram matrix applied to w
        Gw = nufft(nufft_adj(w, st), st);

        % Dense k-space (large Gw) → small weight
        % Sparse k-space (small Gw) → large weight
        den = abs(Gw);
        den = max(den, 1e-4 * max(den));  % prevent divide by zero
        w   = w ./ den;

        % Keep mean = 1 to prevent numerical drift
        w = w / mean(w);

        fprintf('  iter %2d:  min=%.3e  max=%.3e\n', ii, min(w), max(w));
    end

    w = abs(w);

    %% clippingggg
    % The 76 origin copies converge to near-zero (DC is heavily oversampled)
    % The floor raises them slightly so DC is not completely suppressed
    % The cap removes runaway large weights from sparse trajectory gaps
    lo = prctile(w, 2);
    hi = prctile(w, 98);
    w  = max(min(w, hi), lo);
    fprintf('Clipped to [%.4e, %.4e]\n', lo, hi);

    window = max(round(Nk / 80), 5);   % ~120 points for Nk=9576
    w = radialSmooth(kx, ky, w, window);
    fprintf('Radial smooth (win=%d):  min=%.4e  max=%.4e\n', ...
        window, min(w), max(w));
end


%% =================================================================
%%  NEAREST-NEIGHBOUR FALLBACK
%% =================================================================
function w = nn_fallback(kx, ky)
    K = 10;
    N = numel(kx);
    rho  = zeros(N, 1);
    kRep = zeros(N, 1);

    for i = 1:N
        d       = hypot(kx - kx(i), ky - ky(i));
        kRep(i) = sum(d == 0);
        sd      = sort(d);
        kdist   = sd(min(K + kRep(i), N));
        rho(i)  = K / (pi * max(kdist, eps)^2);
    end

    w_geo = 1 ./ max(rho, eps);

    % Smooth before kRep correction
    window = max(round(N / 50), 5);
    w_geo  = radialSmooth(kx, ky, w_geo, window);

    w   = w_geo ./ kRep;
    cap = prctile(w, 98);
    w   = min(w, cap);
    fprintf('NN fallback: min=%.4e  max=%.4e\n', min(w), max(w));
end


%% =================================================================
%%  UTILITIES
%% =================================================================
function ws = radialSmooth(kx, ky, w, window)
% Sort points by radius, apply moving average, map back to original order
    [~, si] = sort(hypot(kx, ky));
    tmp     = movmean(w(si), window);
    ws      = w;
    ws(si)  = tmp;
end

function [kx, ky] = readKTrajectory(f)
    T    = readtable(f);
    cols = lower(string(T.Properties.VariableNames));
    ix   = find(cols == "kx", 1);
    iy   = find(cols == "ky", 1);
    assert(~isempty(ix) && ~isempty(iy), 'Kx/Ky columns not found.');
    kx = double(T{:, ix});
    ky = double(T{:, iy});
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