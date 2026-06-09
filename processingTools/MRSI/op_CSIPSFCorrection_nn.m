function MRSIStruct = op_CSIPSFCorrection_nn(MRSIStruct, k_space_file, NameValueArgs)
%OP_CSIPSFCORRECTION_JN  Density compensation — k-nearest-neighbour method.
%
%   Density estimated as K / (pi * d_K^2) where d_K is the distance
%   to the K-th distinct neighbour.


arguments
    MRSIStruct    (1,1) struct
    k_space_file  (1,:) char {mustBeFile}
    NameValueArgs.numNeighbors  (1,1) double  = 10
    NameValueArgs.isPlotWeights (1,1) logical = false
end

fprintf('\n=== DENSITY COMPENSATION (nearest neighbour) ===\n');

%% --- Guard: must run before spatial FT ---
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
w_raw = dcf_nearest(kx, ky, NameValueArgs.numNeighbors);

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
    title('DCF weights (nearest neighbour)'); xlabel('kx'); ylabel('ky');
end

%% --- Reshape weights to match data and apply ---
W_mat = reshape(w, [Nkpts, Nkshot]);
if hasAvg && ndims(MRSIStruct.data) == 5
    fullW = repmat(reshape(W_mat, [1 1 1 Nkpts Nkshot]), [Nt Ncoils Nave 1 1]);
else
    fullW = repmat(reshape(W_mat, [1 1 Nkpts Nkshot]),   [Nt Ncoils 1 1]);
end
fullW = cast(fullW, 'like', MRSIStruct.data);

MRSIStruct.data = MRSIStruct.data .* fullW;

%% --- Store metadata ---
MRSIStruct.densityComp = struct( ...
    'method',       'nearest', ...
    'numNeighbors', NameValueArgs.numNeighbors, ...
    'weights',      w(:), ...
    'weights_raw',  w_raw(:), ...
    'kx',           kx(:), ...
    'ky',           ky(:), ...
    'weight_matrix',W_mat);

fprintf('=== DCF COMPLETE (nearest neighbour) ===\n\n');
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
        kdist   = sd(min(K + kRep(i), N));             %skip duplicates justincase outofbound
        rho(i)  = K / (pi * max(kdist, eps)^2);         %same as jamie's/lubna's code
    end

    %geometric weight (ignoring duplicates)
    w_geo = 1 ./ max(rho, eps);

    %radial smoothing (geometry only, kRep not yet applied)
    window = max( ...
        round(N / 50), 5);
    w_geo  = radialSmooth(kx, ky, w_geo, window); %%this is technically cheating so talk to mark and jamie

    %divide by kRep so each duplicate gets its fair share
    w = w_geo ./ kRep;

    %clip the max 2% -cap outliers
    cap = prctile(w, 98);
    w   = min(w, cap);
end



function ws = radialSmooth(kx, ky, w, window)
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
    assert(~isempty(ix) && ~isempty(iy), ...
        'Kx/Ky columns not found. Available: %s', ...
        strjoin(T.Properties.VariableNames, ', '));
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

function printRadialProfile(kx, ky, w)
    R     = hypot(kx, ky);
    Rmax  = max(R);
    edges = linspace(0, Rmax, 11);
    fprintf('[DEBUG] --- Radial profile (nearest) ---\n');
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