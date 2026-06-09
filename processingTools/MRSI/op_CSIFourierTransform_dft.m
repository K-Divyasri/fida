
function MRSIStruct = op_CSIFourierTransform(MRSIStruct, k_file, fourierTransform)
    arguments
        MRSIStruct (1,1) struct
        k_file (1,:) char {mustBeFileorDefault} = ""
        fourierTransform.spatial (1,1) logical {mustHaveSpatial(fourierTransform.spatial, MRSIStruct)}
        fourierTransform.spectral (1,1) logical {mustHaveSpectral(fourierTransform.spectral, MRSIStruct)}
    end

    % -- Match Brenden
    % 's defaulting logic exactly --
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




%% -------------------------------------------------------------------------
% applySlowFourierTranformMatrix: chunk over time blocks, same stuff
%% -------------------------------------------------------------------------
function image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, NPtemporal, kPtsPerCycle)
    yLength = length(getCoordinates(MRSIStruct, 'y'));
    xLength = length(getCoordinates(MRSIStruct, 'x'));

    %image dimensions after fourier tranform is time, y, x, and extras.
    imageDimensions = [NPtemporal, yLength, xLength, ...
                                        getSizeFromDimensions(MRSIStruct, {'extras'})];
    image = zeros(imageDimensions);
    
    for iPoint = 1:NPtemporal
        startingPoint = (iPoint - 1)*kPtsPerCycle + 1;
        endingPoint = iPoint * kPtsPerCycle;
        kSpaceSlice = data(startingPoint:endingPoint, :, :);
        vectorizedSlice = reshape(kSpaceSlice, [], size(kSpaceSlice,3));
        ftVectorizedSlice = sftOperator*vectorizedSlice;

        imageSlice = reshape(ftVectorizedSlice, [yLength, xLength, size(ftVectorizedSlice,2)]);
        image(iPoint, :, :, :) = imageSlice;
    end
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






