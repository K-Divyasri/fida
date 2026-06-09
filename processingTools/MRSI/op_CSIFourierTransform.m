
function MRSIStruct = op_CSIFourierTransform(MRSIStruct, k_file, fourierTransform)
    arguments
        MRSIStruct (1,1) struct
        k_file (1,:) char {mustBeFileorDefault} = ""
        fourierTransform.spatial (1,1) logical {mustHaveSpatial(fourierTransform.spatial, MRSIStruct)}
        fourierTransform.spectral (1,1) logical {mustHaveSpectral(fourierTransform.spectral, MRSIStruct)}
    end

    fourierTransform = setDefaultFlags(fourierTransform, MRSIStruct);

    % --- Spatial Transform (if flagged) ---F
    if fourierTransform.spatial
        disp('Calculating spatial dimension');
        if (k_file == "")
            % Cartesian => fast Fourier Transform
            MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct);
            Nt = MRSIStruct.sz(MRSIStruct.dims.t);
            MRSIStruct.spectralDwellTime = MRSIStruct.adcDwellTime;      % 1 time sample per ADC tick
            MRSIStruct.spectralWidth     = 1 / MRSIStruct.spectralDwellTime;
            MRSIStruct.spectralTime      = (0:Nt-1) * MRSIStruct.spectralDwellTime;
            MRSIStruct.adcTime      = (0:Nt-1) * MRSIStruct.adcDwellTime;
        else
            dComp=MRSIStruct;
            if ndims(dComp.data) == 5
                % --- 5D case ---
                % Original: [576 16 4 126 63]
                % Want: bring 126 to the front → permute to [126 576 16 4 63]
                % Then reshape 126*576 = 72576 → [72576 16 4 63]
                dComp.data = permute(dComp.data, [4 1 2 3 5]);
                dComp.data = reshape(dComp.data, [MRSIStruct.sz(4)*MRSIStruct.sz(1), MRSIStruct.sz(2), MRSIStruct.sz(3), MRSIStruct.sz(5)]);
                dComp.sz = size(dComp.data);
            
                % Update dims
                dComp.dims.t = 1;
                dComp.dims.coils = 2;
                dComp.dims.averages = 3;
                dComp.dims.ky = 4;
            
                % Zero out everything else
                fn = fieldnames(dComp.dims);
                for i = 1:numel(fn)
                    if ~ismember(fn{i}, {'t','coils','averages','ky'})
                        dComp.dims.(fn{i}) = 0;
                    end
                end
            
            elseif ndims(dComp.data) == 4
                % --- 4D case ---
                % Original: [576 16 126 63]
                % Want: bring 126 to the front → permute to [126 576 16 63]
                % Then reshape 126*576 = 72576 → [72576 16 63]
                dComp.data = permute(dComp.data, [3 1 2 4]);
                dComp.data = reshape(dComp.data, [MRSIStruct.sz(3)*MRSIStruct.sz(1), MRSIStruct.sz(2), MRSIStruct.sz(4)]);
                dComp.sz = size(dComp.data);
            
                % Update dims
                dComp.dims.t = 1;
                dComp.dims.coils = 2;
                dComp.dims.ky = 3;
            
                % Zero out everything else
                fn = fieldnames(dComp.dims);
                for i = 1:numel(fn)
                    if ~ismember(fn{i}, {'t','coils','ky'})
                        dComp.dims.(fn{i}) = 0;
                    end
                end
            
            else
                warning('Unexpected data dimensionality: %dD', ndims(dComp.data));
            end
            MRSIStruct=dComp;
            % Non-Cartesian => slow transform
            [kTable, kArray]  = readKFile(k_file);
            kPtsPerCycle      = getKPtsPerCycle(kTable);
            NPtemporal        = getTemporalPts(kTable, MRSIStruct);
            
            MRSIStruct = slowFourierTransform(MRSIStruct, kArray, kPtsPerCycle, NPtemporal);
            MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal);
            Nt = MRSIStruct.sz(MRSIStruct.dims.t);
            MRSIStruct.adcTime      = (0:Nt-1) * MRSIStruct.adcDwellTime;
        end
        MRSIStruct = setFlags(MRSIStruct, 'spatialFT', true);
    end

    % --- Spectral Transform (if flagged) ---
    if fourierTransform.spectral
        disp('Calculating spectral dimension');
        MRSIStruct = fastFourierTransformTime(MRSIStruct); 
    end
end



%% -------------------------------------------------------------------------
% sft2_Operator: kinda, but vectorized
%% -------------------------------------------------------------------------
function sft2_Oper = sft2_Operator(InTraj, OutTraj, Ift_flag)
    % In the original code:
    %   if (~Ift_flag) => exponent = -2*pi*1i
    %   else exponent =  2*pi*1i
    %   if Ift_flag => divide by Nx as well

    if ~Ift_flag
        Expy = -2*pi*1i;    % forward transform exponent
    else
        Expy =  2*pi*1i;    % inverse transform exponent
    end

    NOut   = size(OutTraj,1);
    NIn    = size(InTraj,1);

    % Vectorized approach
    xTerm   = OutTraj(:,1)*InTraj(:,1)';  % size [NOut x NIn]
    yTerm   = OutTraj(:,2)*InTraj(:,2)';  % size [NOut x NIn]
    sft2_Oper = exp(Expy*(xTerm + yTerm));

    if Ift_flag
        sft2_Oper = sft2_Oper / NIn;  % divide by number of input k-points
    end
end


%% -------------------------------------------------------------------------
% mustHaveSpatial / mustHaveSpectral
%% -------------------------------------------------------------------------
function mustHaveSpatial(a, in)
    if isfield(in, 'spatialFT') && (a == true && in.spatialFT == 1)
        error('Spatial Fourier Transform already done!');
    end
end

function mustHaveSpectral(a, in)
    if isfield(in, 'spectral') && (a == true && in.spectral == 1)
        error('Spectral Fourier Transform already done!');
    end
end

function mustBeFileorDefault(file)
    if ~isfile(file) && ~strcmp(file, "")
        error('Invalid k_file, must be an existing file or empty.');
    end
end


%% -------------------------------------------------------------------------
% setDefaultFlags: same stuff
%% -------------------------------------------------------------------------
function fourierTransform = setDefaultFlags(fourierTransform, in)
    if ~isfield(fourierTransform, 'spatial')
        if in.flags.spatialFT
            fourierTransform.spatial = 0;
        else
            fourierTransform.spatial = 1;
        end
    end
    if ~isfield(fourierTransform, 'spectral')
        if in.flags.spectralFT
            fourierTransform.spectral = 0;
        else
            fourierTransform.spectral = 1;
        end
    end
end


%% -------------------------------------------------------------------------
% applyFastFourierTransformSpatial: same stuff
%% -------------------------------------------------------------------------
function MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct)
    disp('Applying fast fourier transform (Cartesian)');

    % 1) half-pixel shift
    MRSIStruct = halfPixelShift(MRSIStruct);

    % 2) FFT along x
    data = getData(MRSIStruct);
    xDim = getDimension(MRSIStruct, 'kx');
    if mod(getSizeFromDimensions(MRSIStruct, {'kx'}), 2) == 1
        data = circshift(data, 1, xDim);
    end
    data = fftshift( fft( fftshift(data, xDim), [], xDim ), xDim);

    % 3) FFT along y
    yDim = getDimension(MRSIStruct, 'ky');
    if mod(getSizeFromDimensions(MRSIStruct, {'ky'}), 2) == 1
        data = circshift(data, 1, yDim);
    end
    data = fftshift( fft( fftshift(data, yDim), [], yDim ), yDim);

    MRSIStruct = setData(MRSIStruct, data);

    % 4) re-label dims
    MRSIStruct = setDimension(MRSIStruct, 'x',  getDimension(MRSIStruct, 'kx'));
    MRSIStruct = setDimension(MRSIStruct, 'y',  getDimension(MRSIStruct, 'ky'));
    MRSIStruct = setDimension(MRSIStruct, 'z',  getDimension(MRSIStruct, 'kz'));
    MRSIStruct = setDimension(MRSIStruct, 'kx', 0);
    MRSIStruct = setDimension(MRSIStruct, 'ky', 0);
    MRSIStruct = setDimension(MRSIStruct, 'kz', 0);
end


function MRSIStruct = fastFourierTransformTime(MRSIStruct)
    data = getData(MRSIStruct);
    timeDimension = getDimension(MRSIStruct, 't');

    % Fourier transform in the spectral domain (time -> frequency)
    data = fftshift(ifft(data, [], timeDimension), timeDimension);
    MRSIStruct = setData(MRSIStruct, data);

    % Compute PPM axis
    ppm = calculatePPM(MRSIStruct);
    if strcmp(MRSIStruct.nucleus,'1H')
        ppm = ppm + 4.65;
    end
    MRSIStruct = setPPM(MRSIStruct, ppm);

    % Update flags
    MRSIStruct = setFlags(MRSIStruct, 'spectralFT', true);

    % ---- NEW: Rename dimension label from 't' to 'f' ----
    fDim = getDimension(MRSIStruct, 't');  % Get dimension index for 't'
    MRSIStruct = setDimension(MRSIStruct, 'f', fDim);  % Set same index for 'f'
    MRSIStruct = setDimension(MRSIStruct, 't', 0);      % Clear old 't' label
end

%% -------------------------------------------------------------------------
% slowFourierTransform: same chunk approach as Brenden, preserving dims
%% -------------------------------------------------------------------------
function MRSIStruct = slowFourierTransform(MRSIStruct, kTrajectory, kPtsPerCycle, NPtemporal)
    [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct);

    % Build slow transform operator
    sftOperator = sft2_Operator(kTrajectory, imageTrajectory, 1);
    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t','ky'});

    data = getData(MRSIStruct);

    % Apply chunked slow transform
    image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle);
    MRSIStruct = setData(MRSIStruct, image);

    % Re-label dims
    kyDimension = getDimension(MRSIStruct, 'ky');
    prevPermute = removeDimPrevPermute(prevPermute, kyDimension);
    prevPermute = addDimPrevPermute(prevPermute, 'y', kyDimension);
    prevPermute = addDimPrevPermute(prevPermute, 'x', kyDimension + 1);

    % Match Brenden's logic for final size
    prevSize(1) = NPtemporal;
    prevSize(2) = length(yCoordinates);
    prevSize    = [prevSize(1:2), length(xCoordinates), prevSize(3:end)];

    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end




function image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle)
    yLength = length(getCoordinates(MRSIStruct, 'y'));
    xLength = length(getCoordinates(MRSIStruct, 'x'));

    % sftOperator = (1/Nk) * E^H    [Np x Nk]   from sft2_Operator(Ift_flag=1)
    % Rescale to get actual E^H and E:
    Nk = size(sftOperator, 2);    % = 9576
    EH = Nk * sftOperator;        % [Np x Nk]
    E  = EH';                     % [Nk x Np]
    Np = size(EH, 1);             % = 2304

    % Build normal equation matrix
    EHE = EH * E;                 % [Np x Np]

    %% --- Find optimal lambda via L-curve convexity test ---
    %lambda = find_lambda_lcurve(data, EH, E, EHE, Np, Nk, kPtsPerCycle);
    lambda = 4e-3;
    %% --- Build Tikhonov reconstruction matrix ---
    % B = (E^H*E + lambda*I)^{-1} * E^H    [Np x Nk]
    B = (EHE + lambda * eye(Np)) \ EH;
    fprintf('Tikhonov: Nk=%d  Np=%d  lambda=%.6e\n', Nk, Np, lambda);

    %% --- Reconstruct all time points ---
    imageDimensions = [NPtemporal, yLength, xLength, ...
                       getSizeFromDimensions(MRSIStruct, {'extras'})];
    image = zeros(imageDimensions);

    for iPoint = 1:NPtemporal
        startingPoint   = (iPoint - 1)*kPtsPerCycle + 1;
        endingPoint     =  iPoint * kPtsPerCycle;
        kSpaceSlice     = data(startingPoint:endingPoint, :, :);
        vectorizedSlice = reshape(kSpaceSlice, [], size(kSpaceSlice, 3));

        % Apply Tikhonov operator
        ftVectorizedSlice = B * vectorizedSlice;

        imageSlice = reshape(ftVectorizedSlice, [yLength, xLength, size(ftVectorizedSlice, 2)]);
        image(iPoint, :, :, :) = imageSlice;
    end
end


%% -------------------------------------------------------------------------
%  find_lambda_lcurve: L-curve convexity test using data already in workspace
%  Uses the first temporal frame (iPoint=1) for the sweep
%% -------------------------------------------------------------------------
function lambda_opt = find_lambda_lcurve(data, EH, E, EHE, Np, Nk, kPtsPerCycle)

    fprintf('  [lambda] Starting L-curve sweep...\n');

    %% --- Eigendecompose E^H*E ---
    % Done once here — reused for all 300 lambda values in the sweep
    % Much faster than solving a linear system 300 times
    fprintf('  [lambda] Eigendecomposing [%d x %d]...\n', Np, Np);
    [V, D]   = eig(EHE);
    dv       = real(diag(D));
    [dv, si] = sort(dv, 'descend');
    V        = V(:, si);

    fprintf('  [lambda] Eigenvalues: max=%.4e  min=%.4e  null(eig<1)=%d/%d\n', ...
        dv(1), dv(end), sum(dv < 1), Np);

    %% --- Extract first temporal frame as the k-space vector for sweep ---
    % data is [Nk_total x nCoils x nExtras]
    % First frame = rows 1 to kPtsPerCycle, all petals are in that block
    y_frame = data(1:kPtsPerCycle, :, :);    % [kPtsPerCycle x nCoils x nExtras]
    y       = y_frame(:, 1, 1);              % [kPtsPerCycle x 1] — coil 1 only

    % If data has more shots (petals), we need the full Nk vector
    % data is arranged as [kPtsPerCycle * NPtemporal x nCoils x nExtras]
    % The first temporal frame spans rows 1:kPtsPerCycle but that is only
    % one petal. We need all petals for one time point.
    % Re-extract: all kPtsPerCycle rows across all shots for iPoint=1
    % Actually data layout: rows go [shot1_t1, shot1_t2, ..., shot2_t1, ...]
    % so for iPoint=1: rows 1:kPtsPerCycle is one shot, not all Nk points.
    % Use the full first block across all shots:
    nShots = Nk / kPtsPerCycle;
    y_all  = zeros(Nk, 1);
    for shot = 1:nShots
        row_start           = (shot-1)*kPtsPerCycle + 1;   % in the data matrix
        row_end             = shot * kPtsPerCycle;
        kk_start            = (shot-1)*kPtsPerCycle + 1;   % in y_all
        kk_end              = shot * kPtsPerCycle;
        y_all(kk_start:kk_end) = data(row_start:row_end, 1, 1);
    end
    y = y_all;    % [Nk x 1]

    fprintf('  [lambda] y for sweep: [%d x 1]  norm=%.4e\n', numel(y), norm(y));

    %% --- Project data onto eigenbasis ---
    % c(i) = how much of spatial mode i is in the data
    c = V' * (EH * y);    % [Np x 1]

    %% --- Lambda sweep: 300 values ---
    lambdas     = logspace(-12, 4, 300) * dv(1);
    sol_norm_sq = zeros(1, numel(lambdas));
    res_norm_sq = zeros(1, numel(lambdas));

    fprintf('  [lambda] Sweeping %d values...\n', numel(lambdas));

    for i = 1:numel(lambdas)
        % Tikhonov solution using eigendecomposition shortcut:
        % x = V * (c ./ (d + lambda))  — no matrix inversion needed
        x_lam          = V * (c ./ (dv + lambdas(i)));
        sol_norm_sq(i) = norm(x_lam)^2;        % ||x||^2 — decreases with lambda
        res_norm_sq(i) = norm(E*x_lam - y)^2;  % ||Ex-y||^2 — increases with lambda
    end

    %% --- Convexity test on log-log L-curve ---
    % Log space is used for numerical stability —
    % the L-curve is nearly a right angle on linear scale,
    % making gradient() unreliable there
    lr = log10(res_norm_sq);   % log residual norm squared
    ls = log10(sol_norm_sq);   % log solution norm squared

    % Arc length parameterisation
    % Points are not evenly spaced along the curve in (lr,ls) space.
    % Parameterising by arc length s gives reliable derivatives regardless
    % of how fast we move along the curve at each lambda value.
    dx = diff(lr);
    dy = diff(ls);
    ds = sqrt(dx.^2 + dy.^2);
    s  = [0, cumsum(ds)];

    % First derivatives d/ds
    dlr  = gradient(lr, s);
    dls  = gradient(ls, s);

    % Second derivatives d^2/ds^2
    d2lr = gradient(dlr, s);
    d2ls = gradient(dls, s);

    % Signed curvature of parametric curve (lr(s), ls(s)):
    %   kappa = (lr' * ls'' - ls' * lr'') / (lr'^2 + ls'^2)^(3/2)
    %
    % Sign tells us the bending direction:
    %   kappa > 0 : convex bend  = genuine L-corner  → keep
    %   kappa < 0 : concave bend = spurious corner    → reject
    kappa_signed = (dlr .* d2ls - dls .* d2lr) ./ ...
                   (dlr.^2 + dls.^2).^1.5;

    % Mask 1: endpoints — gradient() less accurate at boundaries
    convex_mask            = kappa_signed > 0;
    convex_mask(1:5)       = false;
    convex_mask(end-4:end) = false;

    % Mask 2: over-regularized vertical arm
    % When lambda is huge, solution is crushed to near zero —
    % these points are not the real corner
    zero_mask              = sol_norm_sq < 0.001 * max(sol_norm_sq);
    convex_mask(zero_mask) = false;

    fprintf('  [lambda] Convex points: %d / %d\n', sum(convex_mask), numel(lambdas));

    % Find maximum curvature among convex points only
    kappa_for_max               = kappa_signed;
    kappa_for_max(~convex_mask) = -inf;

    [kappa_max, ci] = max(kappa_for_max);
    lambda_opt      = lambdas(ci);

    %% --- Print result ---
    fprintf('\n');
    fprintf('  ╔══════════════════════════════════════════╗\n');
    fprintf('  ║  OPTIMAL LAMBDA (L-curve convexity test) ║\n');
    fprintf('  ╠══════════════════════════════════════════╣\n');
    fprintf('  ║  lambda_opt  = %-12.6e             ║\n', lambda_opt);
    fprintf('  ║  log10(lam)  = %-8.4f                 ║\n', log10(lambda_opt));
    fprintf('  ║  index       = %d / %d                  ║\n', ci, numel(lambdas));
    fprintf('  ║  kappa_max   = %-12.6e             ║\n', kappa_max);
    fprintf('  ║  ||x||^2     = %-12.4e             ║\n', sol_norm_sq(ci));
    fprintf('  ║  ||Ex-y||^2  = %-12.4e             ║\n', res_norm_sq(ci));
    fprintf('  ╚══════════════════════════════════════════╝\n\n');
end


%% -------------------------------------------------------------------------
% getImageTrajectory
%% -------------------------------------------------------------------------
function [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct)
    xCoordinates = getCoordinates(MRSIStruct, 'x');
    yCoordinates = getCoordinates(MRSIStruct, 'y');

    [xx, yy]      = meshgrid(xCoordinates, yCoordinates);
    imageTrajectory = [xx(:), yy(:)];
end


%% -------------------------------------------------------------------------
% fastFourierTransformTime: same stuff
%% -------------------------------------------------------------------------



%% -------------------------------------------------------------------------
% calculatePPM: same stuff
%% -------------------------------------------------------------------------
function ppmVals = calculatePPM(MRSIStruct)
    gammaVal      = MRSIStruct.gamma;
    spectralWidth = getSpectralWidth(MRSIStruct);
    timeSize      = getSizeFromDimensions(MRSIStruct, {'t'});

    step       = spectralWidth / timeSize;
    lowerBound = -spectralWidth/2 + step/2;
    upperBound =  spectralWidth/2 - step/2;

    frequencyArray = lowerBound : step : upperBound;
    ppmVals        = -frequencyArray / (MRSIStruct.Bo * gammaVal);
end



function MRSIStruct = calculateSpectralValues(MRSIStruct, kPtsPerCycle, NPtemporal)
    fprintf('=== calculateSpectralValues (FIXED) ===\n');
    
    % CRITICAL FIX: Use NPtemporal for spectral calculations, not kPtsPerCycle
    % The spectral dwell time should be based on temporal points, not k-space points
    adcDwellTime = MRSIStruct.adcDwellTime;
    
    % FIXED: Spectral dwell time calculation
    % For Rosette data: spectral_dt = adc_dt * k_points_per_temporal_cycle
    spectralDwellTime = adcDwellTime * kPtsPerCycle;
    
    % Calculate spectral width and time vector
    spectralWidth = 1/spectralDwellTime;
    spectralTime = (0:(NPtemporal-1)) * spectralDwellTime;
    
    fprintf('Debug (FIXED): adcDwellTime = %.6e s\n', adcDwellTime);
    fprintf('Debug (FIXED): kPtsPerCycle = %d\n', kPtsPerCycle);
    fprintf('Debug (FIXED): NPtemporal = %d\n', NPtemporal);
    fprintf('Debug (FIXED): spectralDwellTime = %.6e s\n', spectralDwellTime);
    fprintf('Debug (FIXED): spectralWidth = %.2f Hz\n', spectralWidth);
    fprintf('Debug (FIXED): spectralTime length = %d\n', length(spectralTime));
    
    % Set values in structure
    MRSIStruct.spectralWidth = spectralWidth;
    MRSIStruct.spectralDwellTime = spectralDwellTime;
    MRSIStruct.spectralTime = spectralTime;
    
    fprintf('=== calculateSpectralValues (FIXED) END ===\n');
end





%% -------------------------------------------------------------------------
% halfPixelShift: same stuff
%% -------------------------------------------------------------------------
function MRSIStruct = halfPixelShift(MRSIStruct)
    kx = getCoordinates(MRSIStruct, 'kx');
    ky = getCoordinates(MRSIStruct, 'ky');

    halfPixelX = getVoxSize(MRSIStruct, 'x')/2;
    halfPixelY = getVoxSize(MRSIStruct, 'y')/2;

    kShift = kx*halfPixelX + ky'*halfPixelY;

    [MRSIStruct, prevPerm, prevSz] = reshapeDimensions(MRSIStruct, {'ky','kx'});
    data = getData(MRSIStruct);

    data = data .* exp(-1i*2*pi*kShift);

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPerm, prevSz);
end






