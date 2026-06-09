function ftSpatial = op_NUFFTSpatial1(dComp, kFile_path)
% op_NUFFTSpatial1 
% Accepts your native data structure directly:
%   dComp:   [576×16×6×126×76]   dims: t=1,coils=2,avgs=3,kpts=4,kshot=5
%   dComp_w: [576×16×126×76]     dims: t=1,coils=2,kpts=3,kshot=4


    fprintf('\n=== START: op_NUFFTSpatial1 ===\n');
    fprintf('Input: %s\n', mat2str(size(dComp.data)));

    % ── Step 1: Reshape 
    % [576,16,6,126,76] → [72576,16,6,76]  dims: t=1,coils=2,avgs=3,ky=4
    % [576,16,126,76]   → [72576,16,76]    dims: t=1,coils=2,ky=3

    if ndims(dComp.data) == 5
        Nt           = dComp.sz(dComp.dims.t);           % 576
        kPtsPerCycle = dComp.sz(dComp.dims.kpts);        % 126
        nCoil        = dComp.sz(dComp.dims.coils);       % 16
        nAvg         = dComp.sz(dComp.dims.averages);    % 6
        nKy          = dComp.sz(dComp.dims.kshot);       % 76

        % permute [kpts, t, coils, avgs, shot] then merge kpts*t
        A = permute(dComp.data, [4, 1, 2, 3, 5]);        % [126,576,16,6,76]
        A = reshape(A, [kPtsPerCycle*Nt, nCoil, nAvg, nKy]); % [72576,16,6,76]

        dComp.data         = A;
        dComp.sz           = double([size(A,1), size(A,2), size(A,3), size(A,4)]);
        dComp.dims.t       = 1;
        dComp.dims.coils   = 2;
        dComp.dims.averages= 3;
        dComp.dims.ky      = 4;
        dComp.dims.kpts    = 0;
        dComp.dims.kshot   = 0;

    elseif ndims(dComp.data) == 4
        Nt           = dComp.sz(dComp.dims.t);           % 576
        kPtsPerCycle = dComp.sz(dComp.dims.kpts);        % 126
        nCoil        = dComp.sz(dComp.dims.coils);       % 16
        nKy          = dComp.sz(dComp.dims.kshot);       % 76

        % permute [kpts, t, coils, shot] then merge kpts*t
        A = permute(dComp.data, [3, 1, 2, 4]);           % [126,576,16,76]
        A = reshape(A, [kPtsPerCycle*Nt, nCoil, nKy]);   % [72576,16,76]

        dComp.data      = A;
        dComp.sz        = double([size(A,1), size(A,2), size(A,3)]);
        dComp.dims.t    = 1;
        dComp.dims.coils= 2;
        dComp.dims.ky   = 3;
        dComp.dims.kpts = 0;
        dComp.dims.kshot= 0;
    else
        error('op_NUFFTSpatial1: unsupported ndims=%d', ndims(dComp.data));
    end

    fprintf('After reshape: %s\n', mat2str(size(dComp.data)));
    fprintf('dims: t=%d coils=%d avgs=%d ky=%d\n', ...
        dComp.dims.t, dComp.dims.coils, dComp.dims.averages, dComp.dims.ky);

    % ── Step 2: Read k-file (1/mm) ───────────────────────────────────────
    [~, kXY] = readKFile_simple(kFile_path);   % [Nk×2] = [Kx, Ky] in 1/mm
    kXY = double(kXY(:, 1:2));
    Nk  = size(kXY, 1);

    % ── Step 3: Sizes ─────────────────────────────────────────────────────
    sz         = dComp.sz;
    dims       = dComp.dims;
    Nt_total   = sz(dims.t);            % 72576
    nKy        = sz(dims.ky);           % 76
    kPtsPerCycle = Nk / nKy;           % 9576/76 = 126
    NPtemporal = Nt_total / kPtsPerCycle; % 72576/126 = 576

    xCoords = getCoordinates(dComp, 'x');   % [48] mm
    yCoords = getCoordinates(dComp, 'y');   % [48] mm
    Nx = numel(xCoords);
    Ny = numel(yCoords);

    fprintf('Nk=%d  kPtsPerCycle=%d  NPtemporal=%d  Nx=%d  Ny=%d\n', ...
        Nk, kPtsPerCycle, NPtemporal, Nx, Ny);

    % ── Step 4: reshapeDimensions — identical to DFT path ─────────────────
    % [72576,16,6,76] → [72576, 76, 96]   Nextra = coils*avgs
    % [72576,16,76]   → [72576, 76, 16]   Nextra = coils
    [dComp, prevPermute, prevSize] = reshapeDimensions(dComp, {'t', 'ky'});
    X      = getData(dComp);               % [72576, 76, Nextra]
    Nextra = size(X, 3);
    fprintf('After reshapeDims: %s  Nextra=%d\n', mat2str(size(X)), Nextra);

    % ── Step 5: NUFFT init ────────────────────────────────────────────────
    % Use [Ky, Kx] axis convention — confirmed correct by op_NUFFTSpatial calibration.
    % Omega in radians/index: om = 2π * k(1/mm) * Δ(mm)
    dx = median(abs(diff(xCoords(:))));
    dy = median(abs(diff(yCoords(:))));
    x_shift = median((0:Nx-1).' - xCoords(:)/dx);
    y_shift = median((0:Ny-1).' - yCoords(:)/dy);
    n_shift = [y_shift, x_shift];

    om = [2*pi*kXY(:,2)*dy, 2*pi*kXY(:,1)*dx];   % [Ky→y, Kx→x]

    Nd = [Ny, Nx];
    Jd = [6, 6];
    Kd = 2 * Nd;

    fprintf('n_shift=[%.3f, %.3f]  Nd=%s  Kd=%s\n', ...
        n_shift(1), n_shift(2), mat2str(Nd), mat2str(Kd));
    fprintf('Initializing NUFFT...\n');
    st = nufft_init(om, Nd, Jd, Kd, n_shift);
    fprintf('NUFFT initialized.\n');

    % Normalization: calibration confirmed alpha ≈ 1/Nk
    normFactor = Nk;   % = 9576
    fprintf('normFactor = %d (= Nk)\n', normFactor);

    % ── Step 6: Reconstruction loop ───────────────────────────────────────
    img = zeros(NPtemporal, Ny, Nx, Nextra, 'like', double(X));

    for it = 1:NPtemporal
        i0 = (it-1)*kPtsPerCycle + 1;
        i1 =  it   *kPtsPerCycle;

        kSlice = X(i0:i1, :, :);
        Y      = double(reshape(kSlice, [], Nextra));         % [Nk, Nextra]

        Z  = nufft_adj(Y, st) / normFactor;                  % [Ny, Nx, Nextra]
        img(it, :, :, :) = reshape(Z, [Ny, Nx, Nextra]);

    end
    fprintf('Reconstruction done: %s\n', mat2str(size(img)));

    % ── Step 7: reshapeBack — identical to DFT path ───────────────────────
    dComp = setData(dComp, double(img));

    kyDim       = getDimension(dComp, 'ky');
    prevPermute = removeDimPrevPermute(prevPermute, kyDim);
    prevPermute = addDimPrevPermute(prevPermute, 'y', kyDim);
    prevPermute = addDimPrevPermute(prevPermute, 'x', kyDim + 1);
    prevSize(1) = NPtemporal;
    prevSize(2) = Ny;
    prevSize    = [prevSize(1:2), Nx, prevSize(3:end)];
    dComp       = reshapeBack(dComp, prevPermute, prevSize);

    fprintf('After reshapeBack: %s\n', mat2str(size(dComp.data)));

    % ── Step 8: Spectral values ───────────────────────────────────────────
    adcDT  = getAdcDwellTime(dComp);
    specDT = adcDT * kPtsPerCycle;
    dComp  = setSpectralWidth(dComp,     1/specDT);
    dComp  = setSpectralDwellTime(dComp, specDT);
    dComp  = setSpectralTime(dComp,      0:specDT:specDT*(NPtemporal-1));
    dComp  = setDimension(dComp, 'kx', 0);
    dComp  = setDimension(dComp, 'ky', 0);
    dComp  = setFlags(dComp, 'spatialFT', true);

    fprintf('spectralDwellTime=%.4e s  spectralWidth=%.2f Hz\n', specDT, 1/specDT);
    fprintf('=== END: op_NUFFTSpatial1 ===\n\n');

    ftSpatial = dComp;
end


% ── Helper: readKFile_simple ──────────────────────────────────────────────
function [kTable, kArray] = readKFile_simple(kFile)
    kTable = []; kArray = [];
    if isempty(kFile) || ~isfile(kFile), return; end
    try
        T = readtable(kFile);
        kTable = T;
        vn = lower(string(T.Properties.VariableNames));
        ix = find(vn == "kx", 1);
        iy = find(vn == "ky", 1);
        if ~isempty(ix) && ~isempty(iy)
            kArray = [T{:,ix}, T{:,iy}];
        elseif width(T) >= 3
            kArray = [T{:,2}, T{:,3}];
        end
    catch
        A = readmatrix(kFile);
        if size(A,2) >= 3
            kArray = A(:,2:3);
            kTable = A;
        end
    end
end