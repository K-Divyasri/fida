function MRSIStruct = op_CSIPSFCorrection_jn2(MRSIStruct, k_space_file, NameValueArgs)
% Density compensation for CSI/MRSI with support for concentric rings layouts.
% Now supports:
%   (A) ky-mode:        dims.ky > 0  -> expects weights on [t x ky]
%   (B) kgrid-mode:     dims.kpts>0 & dims.kshot>0 -> weights on [kpts x kshot]
%
% In kgrid-mode (your timeCombined_rs / timeCombined_rs_w), weights are computed
% on the [kpts x kshot] grid from the k-file and broadcast over all other dims.
%
% USAGE:
%   out = op_CSIPSFCorrection_jn2(in, k_file, 'dcfMethod','kNN', 'modelType','Uniform', ...)
%
% Key NameValueArgs:
%   dcfMethod : 'kNN' (default) | 'PipeMenon'
%   modelType : 'Uniform' (default) | 'Gaussian' | 'FlatEdge'
%   pmIters, pmKernel, pmSigma, pmWidth, pmBeta  (Pipe–Menon)
%   kNN, edgeAlpha, capCenter, capPct            (kNN)
%   isPlotRosette, isPlotWeights                 (plots)
%
% OUTPUT:
%   MRSIStruct : input struct with .data multiplied by DC weights

    arguments
        MRSIStruct (1,1) struct
        k_space_file (1,:) char {mustBeFile}

        NameValueArgs.dcfMethod (1,:) {mustBeMember(NameValueArgs.dcfMethod,{'kNN','PipeMenon'})} = 'kNN'
        NameValueArgs.modelType (1,:) {mustBeMember(NameValueArgs.modelType,{'Uniform','Gaussian','FlatEdge'})} = 'Uniform'
        NameValueArgs.sigma (1,1) double = sqrt(-(0.10^2)/(2*log(0.01)))   % ~0.466
        NameValueArgs.steep (1,1) double = -60

        NameValueArgs.kNN (1,1) double {mustBeInteger, mustBePositive} = 10
        NameValueArgs.edgeAlpha (1,1) double {mustBeGreaterThanOrEqual(NameValueArgs.edgeAlpha,0), mustBeLessThanOrEqual(NameValueArgs.edgeAlpha,1)} = 0.7
        NameValueArgs.capCenter (1,1) logical = true
        NameValueArgs.capPct (1,1) double {mustBeGreaterThan(NameValueArgs.capPct,50), mustBeLessThanOrEqual(NameValueArgs.capPct,100)} = 99

        NameValueArgs.pmIters (1,1) double {mustBeInteger, mustBePositive} = 20
        NameValueArgs.pmKernel (1,:) {mustBeMember(NameValueArgs.pmKernel,{'gaussian','kaiser'})} = 'gaussian'
        NameValueArgs.pmSigma (1,1) double = NaN
        NameValueArgs.pmWidth (1,1) double = NaN
        NameValueArgs.pmBeta (1,1) double = 8.6

        NameValueArgs.isPlotRosette (1,1) logical = false
        NameValueArgs.isPlotWeights (1,1) logical = false
    end

    fprintf('\n=== op_CSIPSFCorrection_jn2 (kgrid-aware) ===\n');
    fprintf('Method: %s   Target: %s\n', NameValueArgs.dcfMethod, NameValueArgs.modelType);

    assert(isfield(MRSIStruct,'dims') && isfield(MRSIStruct,'sz'), 'Input must have .dims and .sz');

    % ---- Dimension detection
    tdim   = safe_dim(MRSIStruct.dims, 't');
    kydim  = safe_dim(MRSIStruct.dims, 'ky');      % ky-mode if >0
    kptsd  = safe_dim(MRSIStruct.dims, 'kpts');    % kgrid-mode needs both >0
    kshotd = safe_dim(MRSIStruct.dims, 'kshot');

    assert(tdim>0, 'dims.t must be set.');

    mode_ky     = (kydim  > 0);
    mode_kgrid  = (kptsd  > 0) && (kshotd > 0);

    assert(mode_ky || mode_kgrid, ...
        'Expected either dims.ky > 0 (ky-mode) OR dims.kpts & dims.kshot > 0 (kgrid-mode).');

    % ---- Read k-file
    [kx, ky] = readKfile(k_space_file);
    assert(numel(kx)==numel(ky), 'Kx/Ky length mismatch');
    Ntraj = numel(kx);
    fprintf('K-file samples: %d\n', Ntraj);

    % ---- Unique compaction
    XY = [kx(:), ky(:)];
    [XYu, ia, ic] = unique(XY, 'rows', 'stable');  %#ok<ASGLU>
    Nu = size(XYu,1);
    dupCounts = accumarray(ic, 1, [Nu, 1]);
    numDups   = sum(dupCounts > 1);
    if numDups > 0
        fprintf('Note: %d unique sites have duplicates (max dup count=%d)\n', numDups, max(dupCounts));
    end

    % ---- Compute DCF at unique sites
    switch NameValueArgs.dcfMethod
        case 'kNN'
            w_u = dcf_knn_unique(XYu, NameValueArgs.kNN, NameValueArgs.edgeAlpha, ...
                                 NameValueArgs.capCenter, NameValueArgs.capPct);
        case 'PipeMenon'
            w_u = dcf_pipemenon_unique(XYu, dupCounts, ...
                                       NameValueArgs.pmIters, NameValueArgs.pmKernel, ...
                                       NameValueArgs.pmSigma, NameValueArgs.pmWidth, NameValueArgs.pmBeta);
    end

    % ---- Optional radial profile shaping
    r_all  = hypot(kx, ky); rMaxAll = max(r_all); rNorm = r_all / max(rMaxAll, eps);
    switch lower(NameValueArgs.modelType)
        case 'uniform'
            desired = ones(size(rNorm));
        case 'gaussian'
            desired = exp( -(rNorm.^2) / (2*NameValueArgs.sigma^2) );
        case 'flatedge'
            s  = NameValueArgs.steep; r0 = 0.90;
            desired = 1 ./ (1 + exp(s*(rNorm - r0)));
    end
    desired = desired / median(desired);
    w_all = (w_u(ic)) .* desired;

    % --------------------------- APPLY WEIGHTS ---------------------------
    if mode_kgrid
        % ===== kgrid-mode: expect [kpts x kshot] from k-file =====
        Nkpts  = MRSIStruct.sz(kptsd);
        Nkshot = MRSIStruct.sz(kshotd);
        assert(Ntraj == Nkpts*Nkshot, ...
            'k-file length (%d) must equal kpts*kshot (%d*%d=%d).', Ntraj, Nkpts, Nkshot, Nkpts*Nkshot);

        % Reshape in acquisition order: fastest kpts, then kshot
        W_kgrid = reshape(w_all, [Nkpts, Nkshot]);

        medW = median(W_kgrid(:)); if medW>0, W_kgrid = W_kgrid/medW; end
        fprintf('Weights [kpts x kshot]=[%d x %d]  median=1.0  min=%.3g  max=%.3g\n', ...
            size(W_kgrid,1), size(W_kgrid,2), min(W_kgrid(:)), max(W_kgrid(:)));

        if NameValueArgs.isPlotWeights
            figure('Color','w'); imagesc(W_kgrid); axis image;
            title('[kpts x kshot] weights'); colorbar; xlabel('kshot'); ylabel('kpts');
        end

        % Broadcast to full nd array (dim-aware)
        MRSIStruct = applyW_dimaware_kgrid(MRSIStruct, W_kgrid, kptsd, kshotd);

    else
        % ===== ky-mode: fall back to legacy [t x ky] approach =====
        Nt  = MRSIStruct.sz(tdim);
        Nky = MRSIStruct.sz(kydim);

        assert(mod(Ntraj, Nky)==0, 'k-file length must be divisible by ky count.');
        ADCptsPerKy = Ntraj / Nky;         % samples along readout per ky index

        % Build [t x ky] by repeating the readout pattern to length Nt
        base = reshape(w_all, [ADCptsPerKy, Nky]);        % [kpts-per-ky x ky]
        reps = ceil(Nt / ADCptsPerKy);
        W_tky = repmat(base, [reps, 1]);                   % >= [Nt x ky]
        W_tky = W_tky(1:Nt, :);                            % crop to [Nt x ky]

        medW = median(W_tky(:)); if medW>0, W_tky = W_tky/medW; end
        fprintf('Weights [t x ky]=[%d x %d]  median=1.0  min=%.3g  max=%.3g\n', ...
            size(W_tky,1), size(W_tky,2), min(W_tky(:)), max(W_tky(:)));

        if NameValueArgs.isPlotWeights
            figure('Color','w');
            subplot(1,2,1); plot(mean(W_tky,2)); xlabel('time index'); ylabel('mean over ky'); title('Mean weight vs time');
            subplot(1,2,2); imagesc(W_tky); axis image; title('[t x ky] weights'); colorbar;
        end

        MRSIStruct = applyW_dimaware_tky(MRSIStruct, W_tky, tdim, kydim);
    end
    % --------------------------------------------------------------------

    fprintf('Applied weights (dim-aware). Done.\n');
end

% ====================== Helpers ======================
function d = safe_dim(dims, field)
    if isfield(dims, field), d = dims.(field); else, d = 0; end
end

function [kx, ky] = readKfile(fname)
    T = readtable(fname);
    if ismember('Kx', T.Properties.VariableNames), kx = T.Kx; 
    elseif ismember('kx', T.Properties.VariableNames), kx = T.kx;
    else, error('Kx column not found in k-file'); end
    if ismember('Ky', T.Properties.VariableNames), ky = T.Ky; 
    elseif ismember('ky', T.Properties.VariableNames), ky = T.ky;
    else, error('Ky column not found in k-file'); end
    kx = kx(:); ky = ky(:);
end

% ---------- DCF (kNN) on unique sites ----------
function w_u = dcf_knn_unique(XYu, k, edgeAlpha, capCenter, capPct)
    fprintf('DCF(kNN): k=%d  edgeAlpha=%.2f  capCenter=%d\n', k, edgeAlpha, capCenter);
    Nu = size(XYu,1);
    D  = pdist2(XYu, XYu, 'euclidean');

    r_u  = hypot(XYu(:,1), XYu(:,2));
    rMax = max(r_u);

    dens_u = zeros(Nu,1);
    for j = 1:Nu
        dj = sort(D(j,:));
        dj = dj(dj > 0);
        idx = min(k, numel(dj));
        if idx==0, djk = rMax*1e-6; else, djk = dj(idx); end

        dist2edge = rMax - r_u(j);
        if djk > dist2edge
            baseCorr = 1 + (1 - dist2edge/max(djk, eps));
            corr = 1 + edgeAlpha*(baseCorr - 1);
        else
            corr = 1;
        end
        dens_u(j) = corr * (k / (pi * max(djk, eps)^2));
    end

    if capCenter
        capVal = prctile(dens_u, capPct);
        dens_u = min(dens_u, capVal);
        fprintf('  kNN: capped dens at %g-th pct = %.3g\n', capPct, capVal);
    end

    dens_u = max(dens_u, 1e-6);
    dens_u = dens_u / min(dens_u);
    fprintf('  kNN: density stats (unique): min=%.3g  med=%.3g  max=%.3g\n', min(dens_u), median(dens_u), max(dens_u));

    w_u = 1 ./ dens_u;
    w_u = w_u / median(w_u);
    fprintf('  kNN: weights (unique): min=%.3g  med=1.0  max=%.3g\n', min(w_u), max(w_u));
end

% ---------- DCF (Pipe–Menon) on unique sites ----------
function w_u = dcf_pipemenon_unique(XYu, dupCounts, pmIters, pmKernel, pmSigma, pmWidth, pmBeta)
    Nu = size(XYu,1);
    fprintf('DCF(PipeMenon): iters=%d  kernel=%s\n', pmIters, pmKernel);

    D = pdist2(XYu, XYu, 'euclidean');

    sortedD = sort(D + diag(inf(Nu,1)), 2);
    med4 = median(sortedD(:, min(4, size(sortedD,2))));
    if isnan(pmSigma) || pmSigma<=0, pmSigmaEff = med4 / 1.5; else, pmSigmaEff = pmSigma; end
    if isnan(pmWidth) || pmWidth<=0, pmWidthEff = 2.0 * med4; else, pmWidthEff = pmWidth; end
    fprintf('  kernel params: sigma=%.4g  width=%.4g  beta=%.3g\n', pmSigmaEff, pmWidthEff, pmBeta);

    switch lower(pmKernel)
        case 'gaussian'
            K = exp( - (D.^2) / (2*pmSigmaEff^2) );
        case 'kaiser'
            R = pmWidthEff;
            K = zeros(size(D));
            idx = D < R;
            x = D(idx)/R;
            I0beta = besseli(0, pmBeta);
            K(idx) = besseli(0, pmBeta*sqrt(1 - x.^2)) / max(I0beta, eps);
        otherwise
            error('Unknown pmKernel.');
    end

    w = ones(Nu,1); counts = dupCounts(:);

    for it = 1:pmIters
        rho = K * (w .* counts);
        rho = max(rho, 1e-12);
        w_new = w ./ rho;
        wm = median(w_new); if wm>0, w_new = w_new / wm; end
        w = w_new;
        if mod(it,5)==0 || it==1 || it==pmIters
            fprintf('  iter %2d: w(min/med/max) = [%.3g / %.3g / %.3g]\n', it, min(w), median(w), max(w));
        end
    end

    w_u = w / median(w);
    fprintf('  PipeMenon: final weights (unique): min=%.3g  med=1.0  max=%.3g\n', min(w_u), max(w_u));
end

% ---------- Broadcasting: [kpts x kshot] ----------
function out = applyW_dimaware_kgrid(in, W_kgrid, kptsd, kshotd)
    out  = in;
    data = in.data;
    szD  = size(data);

    Nkpts  = in.sz(kptsd);
    Nkshot = in.sz(kshotd);
    assert(all(size(W_kgrid) == [Nkpts, Nkshot]), ...
        'W_kgrid must be [%d x %d].', Nkpts, Nkshot);

    wSize             = ones(1, ndims(data));
    wSize(kptsd)      = Nkpts;
    wSize(kshotd)     = Nkshot;
    Wnd               = reshape(W_kgrid, wSize);

    reps = szD ./ wSize;
    if any(abs(reps - round(reps)) > 1e-9)
        error('Weights cannot tile to match data size. Check dims/sizes.');
    end
    Wnd = repmat(Wnd, round(reps));

    out.data = data .* Wnd;
end

% ---------- Broadcasting: [t x ky] (legacy path) ----------
function out = applyW_dimaware_tky(in, W_tky, tdim, kydim)
    out  = in;
    data = in.data;
    szD  = size(data);

    Nt  = in.sz(tdim);
    Nky = in.sz(kydim);
    assert(all(size(W_tky) == [Nt, Nky]), ...
        'W_tky must be [%d x %d].', Nt, Nky);

    wSize         = ones(1, ndims(data));
    wSize(tdim)   = Nt;
    wSize(kydim)  = Nky;
    Wnd           = reshape(W_tky, wSize);

    reps = szD ./ wSize;
    if any(abs(reps - round(reps)) > 1e-9)
        error('Weights cannot tile to match data size. Check dims/sizes.');
    end
    Wnd = repmat(Wnd, round(reps));

    out.data = data .* Wnd;
end
