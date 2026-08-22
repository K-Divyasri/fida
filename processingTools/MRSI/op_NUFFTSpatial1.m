function ftSpatial = op_NUFFTSpatial1(dComp, kFile_path, varargin)
% op_NUFFTSpatial1
% Accepts your native data structure directly:
%   dComp:   [576×16×6×126×76]   dims: t=1,coils=2,avgs=3,kpts=4,kshot=5
%   dComp_w: [576×16×126×76]     dims: t=1,coils=2,kpts=3,kshot=4
%
% NAME-VALUE OPTIONS (all optional):
%   'offres'     (false)  EXPERIMENTAL off-resonance-corrected recon (MFI).
%                         Within one temporal point the 126 k-points span
%                         specDT (~0.63 ms); off-carrier spins dephase across
%                         them -> spatial blur / signal loss growing with
%                         |chemical-shift offset|. MFI reconstructs at several
%                         demodulation frequencies and interpolates per
%                         spectral bin to undo it. Default OFF = original path.
%   'nMFI'       (5)      number of demodulation frequencies across the SW.
%   'offresSign' (-1)     sign of the demod phase exp(1i*sign*2*pi*f*tau).
%                         If the ppm response gets WORSE, flip to +1.

    prm = struct('offres', false, 'nMFI', 5, 'offresSign', -1);
    for a = 1:2:numel(varargin)
        prm.(varargin{a}) = varargin{a+1};
    end

    fprintf('\n=== START: op_NUFFTSpatial1 ===\n');
    fprintf('Input: %s\n', mat2str(size(dComp.data)));

    % ── Step 1: Reshape
    if ndims(dComp.data) == 5
        Nt           = dComp.sz(dComp.dims.t);           % 576
        kPtsPerCycle = dComp.sz(dComp.dims.kpts);        % 126
        nCoil        = dComp.sz(dComp.dims.coils);       % 16
        nAvg         = dComp.sz(dComp.dims.averages);    % 6
        nKy          = dComp.sz(dComp.dims.kshot);       % 76

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

    % ── Step 2: Read k-file (1/mm)
    [~, kXY] = readKFile_simple(kFile_path);
    kXY = double(kXY(:, 1:2));
    Nk  = size(kXY, 1);

    % ── Step 3: Sizes
    sz         = dComp.sz;
    dims       = dComp.dims;
    Nt_total   = sz(dims.t);            % 72576
    nKy        = sz(dims.ky);           % 76
    kPtsPerCycle = Nk / nKy;           % 9576/76 = 126
    NPtemporal = Nt_total / kPtsPerCycle; % 72576/126 = 576

    xCoords = getCoordinates(dComp, 'x');
    yCoords = getCoordinates(dComp, 'y');
    Nx = numel(xCoords);
    Ny = numel(yCoords);

    fprintf('Nk=%d  kPtsPerCycle=%d  NPtemporal=%d  Nx=%d  Ny=%d\n', ...
        Nk, kPtsPerCycle, NPtemporal, Nx, Ny);

    % ── Step 4: reshapeDimensions
    [dComp, prevPermute, prevSize] = reshapeDimensions(dComp, {'t', 'ky'});
    X      = getData(dComp);               % [72576, 76, Nextra]
    Nextra = size(X, 3);
    fprintf('After reshapeDims: %s  Nextra=%d\n', mat2str(size(X)), Nextra);

    % ── Step 5: NUFFT init
    dx = median(abs(diff(xCoords(:))));
    dy = median(abs(diff(yCoords(:))));
    x_shift = median((0:Nx-1).' - xCoords(:)/dx);
    y_shift = median((0:Ny-1).' - yCoords(:)/dy);
    n_shift = [y_shift, x_shift];

    om = [2*pi*kXY(:,2)*dy, 2*pi*kXY(:,1)*dx];

    Nd = [Ny, Nx];
    Jd = [6, 6];
    Kd = 2 * Nd;

    fprintf('n_shift=[%.3f, %.3f]  Nd=%s  Kd=%s\n', ...
        n_shift(1), n_shift(2), mat2str(Nd), mat2str(Kd));
    fprintf('Initializing NUFFT...\n');
    st = nufft_init(om, Nd, Jd, Kd, n_shift);
    fprintf('NUFFT initialized.\n');

    normFactor = Nk;
    fprintf('normFactor = %d (= Nk)\n', normFactor);

    % Intra-cycle sample times + spectral params (needed for MFI / spectral hdr)
    adcDT  = getAdcDwellTime(dComp);
    specDT = adcDT * kPtsPerCycle;
    SW     = 1/specDT;

    % ── Step 6: Reconstruction ────────────────────────────────────────────
    if ~prm.offres
        % ---- Original single-frequency recon (unchanged) ----
        img = reconLoopDemod(X, st, kPtsPerCycle, NPtemporal, Ny, Nx, Nextra, ...
                             normFactor, []);
        fprintf('Reconstruction done: %s\n', mat2str(size(img)));
    else
        % ---- EXPERIMENTAL off-resonance MFI recon ----
        tau  = (0:kPtsPerCycle-1).' * adcDT;                 % intra-cycle time [kpts×1]
        fvec = linspace(-SW/2, SW/2, prm.nMFI);              % demod frequencies (Hz)
        Nf   = NPtemporal;
        % spectral frequency axis matching the downstream fftshift(ifft(.)) FT
        fax  = ((0:Nf-1).' - floor(Nf/2)) * (SW/Nf);         % [Nf×1] Hz
        fprintf('MFI: nMFI=%d  SW=%.1f Hz  demod=%s Hz  sign=%+d\n', ...
            prm.nMFI, SW, mat2str(round(fvec)), prm.offresSign);

        Sfinal = complex(zeros(Nf, Ny, Nx, Nextra));
        for m = 1:prm.nMFI
            demod = exp(1i * prm.offresSign * 2*pi * fvec(m) * tau);   % [kpts×1]
            img_m = reconLoopDemod(X, st, kPtsPerCycle, NPtemporal, Ny, Nx, Nextra, ...
                                   normFactor, demod);
            S_m = fftshift(ifft(img_m, [], 1), 1);            % -> spectrum (match FT)
            W   = mfiWeight(fax, fvec, m);                    % [Nf×1] interp weight
            Sfinal = Sfinal + W .* S_m;                       % broadcast over y,x,extra
            fprintf('  MFI node %d/%d (f=%+.0f Hz) done\n', m, prm.nMFI, fvec(m));
        end
        % back to time domain so the downstream spectral FT reproduces Sfinal
        img = fft(ifftshift(Sfinal, 1), [], 1);
        fprintf('MFI reconstruction done: %s\n', mat2str(size(img)));
    end

    % ── Step 7: reshapeBack
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

    % ── Step 8: Spectral values
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


% ── Helper: per-temporal-point NUFFT-adjoint loop, with optional demod ─────
function img = reconLoopDemod(X, st, kPtsPerCycle, NPtemporal, Ny, Nx, Nextra, normFactor, demod)
    img = complex(zeros(NPtemporal, Ny, Nx, Nextra));
    useDemod = ~isempty(demod);
    for it = 1:NPtemporal
        i0 = (it-1)*kPtsPerCycle + 1;
        i1 =  it   *kPtsPerCycle;

        kSlice = X(i0:i1, :, :);                 % [kpts, nKy, Nextra]
        if useDemod
            kSlice = kSlice .* demod;            % demod is [kpts×1], broadcasts
        end
        Y = double(reshape(kSlice, [], Nextra)); % [Nk, Nextra]

        Z = nufft_adj(Y, st) / normFactor;       % [Ny, Nx, Nextra]
        img(it, :, :, :) = reshape(Z, [Ny, Nx, Nextra]);
    end
end


% ── Helper: MFI linear-interpolation weight (partition of unity) ───────────
function W = mfiWeight(fax, fvec, m)
    L = numel(fvec);
    W = zeros(size(fax));
    f = fvec(m);
    if m > 1, lo = fvec(m-1); else, lo = -inf; end
    if m < L, hi = fvec(m+1); else, hi =  inf; end

    % rising edge (lo -> f)
    if m > 1
        idx = fax >= lo & fax <= f;
        W(idx) = (fax(idx) - lo) / (f - lo);
    else
        W(fax <= f) = 1;                 % clamp below range to first node
    end
    % falling edge (f -> hi)
    if m < L
        idx = fax > f & fax <= hi;
        W(idx) = (hi - fax(idx)) / (hi - f);
    else
        W(fax > f) = 1;                  % clamp above range to last node
    end
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
