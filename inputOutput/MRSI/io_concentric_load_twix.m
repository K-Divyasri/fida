function concentric = io_concentric_load_twix(filename)
% io_concentric_load_twix_reordered.m
% Load Siemens TWIX data for concentric ring trajectory sequences
% MODIFIED: Reorders dimensions to [t, coils, averages, kshot, extras]
%
% USAGE:
% concentric = io_concentric_load_twix_reordered(filename)
%
% INPUTS:
% filename = Path to Siemens TWIX (.dat) file
%
% OUTPUTS:
% concentric = Structure with concentric ring data in FID-A format
%              Data dimension order: [t, coils, averages, kshot, extras]

fprintf('Loading concentric ring TWIX data (REORDERED)...\n');
fprintf('File: %s\n', filename);

try
    %% Step 1: Read TWIX file
    twix_obj = readTwixFile(filename);
    
    %% Step 2: Extract raw data
    dOut.data = twix_obj.image();
    data = squeeze(dOut.data);
    
    fprintf('Raw data size: [%s]\n', num2str(size(data)));
    
    %% Step 3: Get version and sequence info
    sqzDims = twix_obj.image.sqzDims;
    sequence = twix_obj.hdr.Config.SequenceFileName;
    
    fprintf('Squeezed dimensions: {%s}\n', strjoin(sqzDims, ', '));
    
    isRosette = false;
    isCartesian = false;
    isConcentric = contains(sequence, 'conc', 'IgnoreCase', true) || ...
                   contains(sequence, 'ring', 'IgnoreCase', true);
    
    fprintf('Sequence: %s\n', sequence);
    fprintf('Is Concentric: %d\n', isConcentric);
    
    %% Step 4: Fill and permute dimensions - MODIFIED FOR NEW ORDER
    spatialFT = false;
    spectralFT = false;
    
    dims = fillDimsField(sqzDims, spatialFT, spectralFT, isCartesian, isConcentric);
    [dims, data] = permuteDims_REORDERED(dims, data, spatialFT, spectralFT, isCartesian, isConcentric);
    
    sz = size(data);
    fprintf('\n=== AFTER REORDERED PERMUTE ===\n');
    fprintf('Data size: [%s]\n', num2str(sz));
    fprintf('Target order: [t, coils, averages, kshot, extras]\n');
    fprintf('Dimension mapping:\n');
    print_dims(dims, sz);
    fprintf('====================\n\n');
    
    %% Step 5: Calculate timing parameters
    [dwelltime, spectralwidth, leftshift] = calculateTimingParameters(twix_obj, dims, sz);
    
    n_time_points = sz(dims.t);
    adcTime = (0:(n_time_points-1)) * dwelltime;
    spectralTime = adcTime;
    spectralDwellTime = dwelltime;
    
    fprintf('Timing: dwell=%.3e s, SW=%.1f Hz, points=%d\n', ...
        dwelltime, spectralwidth, n_time_points);
    
    %% Step 6: Determine matrix dimensions
    [numX, numY, numZ] = determineMatrixDimensions(twix_obj, dims, sz, isCartesian, isConcentric);
    
    fprintf('Matrix dimensions: X=%d, Y=%d, Z=%d\n', numX, numY, numZ);
    
    %% Step 7: Determine subspecs and averages
    [subspecs, rawSubspecs, averages, rawAverages] = ...
        determineSubspecsAndAverages(dims, sz, twix_obj);
    
    fprintf('Averages: %d (raw: %d)\n', averages, rawAverages);
    fprintf('Subspecs: %d (raw: %d)\n', subspecs, rawSubspecs);
    
    %% Step 8: Get nucleus and gamma
    [nucleus, gamma] = getNucleusAndGamma(twix_obj);
    
    %% Step 9: Build the output structure
    concentric = struct();
    
    % Core data
    concentric.data = data;
    concentric.sz = sz;
    concentric.dims = dims;
    
    % Timing parameters
    concentric.adcDwellTime = dwelltime;
    concentric.adcTime = adcTime;
    concentric.spectralWidth = spectralwidth;
    concentric.spectralTime = spectralTime;
    concentric.spectralDwellTime = spectralDwellTime;
    
    % Scanner parameters
    concentric.txfrq = twix_obj.hdr.Meas.lFrequency;
    concentric.scanDate = findScanDate(twix_obj);
    concentric.Bo = twix_obj.hdr.Dicom.flMagneticFieldStrength;
    concentric.nucleus = nucleus;
    concentric.gamma = gamma;
    concentric.seq = sequence;
    concentric.te = twix_obj.hdr.MeasYaps.alTE{1}/1000;
    concentric.tr = twix_obj.hdr.MeasYaps.alTR{1}/1000;
    concentric.pointsToLeftshift = leftshift;
    
    % Spatial parameters
    concentric = findAndSetFov(concentric, twix_obj);
    concentric = calculateVoxelSize(concentric, numX, numY, numZ);
    
    % Acquisition parameters
    concentric.averages = averages;
    concentric.rawAverages = rawAverages;
    concentric.subspecs = subspecs;
    concentric.rawSubspecs = rawSubspecs;
    
    % Geometry
    concentric = findImageOrigin(concentric, twix_obj);
    concentric = calculateVoxelCoodinates(concentric);
    concentric = calculateAffineMatrix(concentric, twix_obj);
    
    % Flags
    concentric = setDefaultFlagValues(concentric, isCartesian);
    
    fprintf('\n=== COORDINATES CHECK ===\n');
    fprintf('Coordinates field exists: %d\n', isfield(concentric, 'coordinates'));
    if isfield(concentric, 'coordinates')
        fprintf('X coords: [%.2f to %.2f], length=%d\n', ...
            concentric.coordinates.x(1), concentric.coordinates.x(end), ...
            length(concentric.coordinates.x));
        fprintf('Y coords: [%.2f to %.2f], length=%d\n', ...
            concentric.coordinates.y(1), concentric.coordinates.y(end), ...
            length(concentric.coordinates.y));
    end
    fprintf('========================\n\n');
    
    fprintf('Data loading complete!\n');
    fprintf('✓ Dimension order: [t=%d, coils=%d, avg=%d, kshot=%d, extras=%d]\n', ...
            dims.t, dims.coils, dims.averages, dims.kshot, dims.extras);
    
catch ME
    fprintf('\n=== ERROR ===\n');
    fprintf('Message: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    fprintf('=============\n');
    rethrow(ME);
end

end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function twix_obj = readTwixFile(filename)
    if(~exist('filename', 'var'))
        twix_obj = mapVBVD;
    else
        twix_obj = mapVBVD(char(filename));
    end

    if isstruct(twix_obj)
        disp('single RAID file detected.');
    elseif iscell(twix_obj)
        disp('multi RAID file detected.');
        RaidLength = length(twix_obj);
        twix_obj = twix_obj{RaidLength};
    end
end

function dims = fillDimsField(sqzDims, isSpatialFT, isSpectralFT, isCartesian, isConcentric)
    % Initialize all possible fields to zero
    keys = {'t','f','coils','averages','timeinterleave','kx','kshot','kz','x','y','z','kpts','kshot','subspec','extras'};
    for k = 1:numel(keys)
        dims.(keys{k}) = 0; 
    end

    if ~isSpatialFT
        if isConcentric
            % CONCENTRIC RINGS: Map raw dimensions
            % Raw order typically: [Col, Cha, Lin, Par, Ave]
            map = containers.Map({'Col','Cha','Lin','Par','Rep','Ave','Seg','Eco'},...
                                 {'t','coils','kshot','extras','averages','averages','kshot','extras'});
        elseif isCartesian
            % Cartesian CSI
            map = containers.Map({'Col','Cha','Ave','Rep','Seg','Phs','Lin','Sli'},...
                                 {'t','coils','averages','timeinterleave','kx','kshot','kshot','kz'});
            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isKey(map,lbl)
                    dims.(map(lbl)) = i;
                elseif dims.extras == 0
                    dims.extras = i;
                else
                    dims.timeinterleave = i;
                end
            end
            return;
        else
            % Rosette/non-Cartesian
            fixed = struct('Col','t','Cha','coils','Ave','averages');
            nextIsKshot = true;
            nextIsExtras = false;

            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isfield(fixed,lbl)
                    dims.(fixed.(lbl)) = i;
                elseif strcmp(lbl,'Set')
                    dims.kshot = i;
                    nextIsKshot = false;
                    nextIsExtras = true;
                elseif nextIsKshot
                    dims.kshot = i;
                    nextIsKshot = false;
                    nextIsExtras = true;
                elseif nextIsExtras
                    dims.extras = i;
                    nextIsExtras = false;
                else
                    dims.timeinterleave = i;
                end
            end
            return;
        end
        
        % Map for concentric
        for i = 1:numel(sqzDims)
            lbl = sqzDims{i};
            if isKey(map, lbl)
                dims.(map(lbl)) = i;
            end
        end
    else
        % After spatial FT
        if isCartesian
            map = containers.Map({'Col','Cha','Ave','Rep','Phs','Lin','Sli'},...
                                 {'t','coils','averages','timeinterleave','x','y','z'});
        else
            map = containers.Map({'Col','Cha','Ave','Rep','X','Y'},...
                                 {'t','coils','averages','kshot','x','y'});
        end
        
        for i = 1:numel(sqzDims)
            lbl = sqzDims{i};
            if isKey(map, lbl)
                dims.(map(lbl)) = i;
            elseif dims.extras == 0
                dims.extras = i;
            else
                dims.timeinterleave = i;
            end
        end
    end
    
    if isSpectralFT && dims.t > 0
        dims.f = dims.t;
        dims.t = 0;
    end
end

function [dims, data] = permuteDims_REORDERED(dims, data, isSpatialFT, isSpectralFT, isCartesian, isConcentric)
    % MODIFIED: Reorder to [t, coils, averages, kshot, extras]
    
    if ~isSpatialFT
        if isConcentric
            % TARGET ORDER: [t, coils, averages, kshot, extras]
            target = {'t', 'coils', 'averages', 'kshot', 'extras'};
        elseif isCartesian
            target = {'t','coils','kx','kshot','kz','averages'};
        else
            target = {'t','coils','averages','kpts','kshot','extras'};
        end
    else
        if isCartesian
            target = {'t','coils','x','y','z','averages'};
        else
            target = {'t','coils','x','y','averages'};
        end
    end
    
    % Build current order and list of present dimensions
    currentOrder = [];
    presentDims = {};
    
    for i = 1:length(target)
        dimVal = dims.(target{i});
        if dimVal > 0
            currentOrder(end+1) = dimVal;
            presentDims{end+1} = target{i};
        end
    end
    
    fprintf('DEBUG: Current dims order in data: %s\n', mat2str(currentOrder));
    fprintf('DEBUG: Target dims order (1:%d): %s\n', length(currentOrder), mat2str(1:length(currentOrder)));
    
    % Permute if necessary
    if ~isequal(currentOrder, 1:length(currentOrder))
        fprintf('DEBUG: Permuting data from %s to %s\n', mat2str(currentOrder), mat2str(1:length(currentOrder)));
        data = permute(data, currentOrder);
    else
        fprintf('DEBUG: No permutation needed - already in correct order\n');
    end
    
    % Update dims indices to reflect new order
    for i = 1:length(presentDims)
        dims.(presentDims{i}) = i;
    end
    
    % Verify final order
    fprintf('DEBUG: Final dimension assignment:\n');
    for i = 1:length(presentDims)
        fprintf('  dims.%s = %d (size: %d)\n', presentDims{i}, i, size(data, i));
    end
end

function [dwelltime, spectralwidth, leftshift] = calculateTimingParameters(twix_obj, dims, sz)
    leftshift = 0;
    dwelltime = NaN;
    
    try
        if isfield(twix_obj.hdr.MeasYaps, 'sRXSPEC') && isfield(twix_obj.hdr.MeasYaps.sRXSPEC, 'alDwellTime')
            dwellArray = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime;
            if iscell(dwellArray) && ~isempty(dwellArray)
                dwelltime = double(dwellArray{1}) * 1e-9;
            elseif isnumeric(dwellArray) && ~isempty(dwellArray)
                dwelltime = double(dwellArray(1)) * 1e-9;
            end
        end
        
        if isnan(dwelltime) && isfield(twix_obj.hdr.MeasYaps, 'alDwellTime')
            dwellArray = twix_obj.hdr.MeasYaps.alDwellTime;
            if iscell(dwellArray) && ~isempty(dwellArray)
                dwelltime = double(dwellArray{1}) * 1e-9;
            elseif isnumeric(dwellArray) && ~isempty(dwellArray)
                dwelltime = double(dwellArray(1)) * 1e-9;
            end
        end
        
        if isnan(dwelltime) && isfield(twix_obj.hdr, 'Phoenix')
            dwelltime = getFieldSafe(twix_obj.hdr.Phoenix, 'alDwellTime', NaN);
            if ~isnan(dwelltime)
                dwelltime = dwelltime * 1e-9;
            end
        end
        
        if isnan(dwelltime) || dwelltime <= 0
            warning('Could not extract valid ADC dwell time from TWIX header. Using default.');
            dwelltime = 5e-6;
            fprintf('  Using dwell time: %.3e s\n', dwelltime);
        end
    catch ME
        warning('Failed to extract timing parameters: %s. Using defaults.', ME.message);
        dwelltime = 5e-6;
    end
    
    spectralwidth = 1 / dwelltime;
end

function [numX, numY, numZ] = determineMatrixDimensions(twix_obj, dims, sz, isCartesian, isConcentric)
    try
        numX = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
        numY = twix_obj.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
        numZ = twix_obj.hdr.MeasYaps.sKSpace.dSliceResolution;
        fprintf('  Matrix from header: %d x %d x %d\n', numX, numY, numZ);
    catch
        warning('Could not extract matrix dimensions from header. Using data dimensions.');
        if isCartesian
            numX = (dims.kx > 0) * sz(max(dims.kx, 1));
            numY = (dims.kshot > 0) * sz(max(dims.kshot, 1));
            numZ = (dims.kz > 0) * sz(max(dims.kz, 1));
        else
            numX = (dims.kshot > 0) * sz(max(dims.kshot, 1));
            numY = numX;
            numZ = 1;
        end
        fprintf('  Matrix from data: %d x %d x %d\n', numX, numY, numZ);
    end
    
    if numX == 0, numX = 1; end
    if numY == 0, numY = 1; end
    if numZ == 0, numZ = 1; end
end

function [subspecs, rawSubspecs, averages, rawAverages] = determineSubspecsAndAverages(dims, sz, twix_obj)
    if dims.subspec > 0
        subspecs = sz(dims.subspec);
        rawSubspecs = subspecs;
    else
        subspecs = 1;
        rawSubspecs = 1;
    end
    
    if dims.averages > 0
        averages = sz(dims.averages);
        rawAverages = averages;
    else
        try
            averages = twix_obj.hdr.MeasYaps.lAverages;
            if averages == 0, averages = 1; end
            rawAverages = averages;
        catch
            averages = 1;
            rawAverages = 1;
        end
    end
end

function [nucleus, gamma] = getNucleusAndGamma(twix_obj)
    nucleus = '';
    try
        if isfield(twix_obj.hdr.Dicom, 'tNucleus')
            nucleus = twix_obj.hdr.Dicom.tNucleus;
        end
        if isempty(nucleus) && isfield(twix_obj.hdr, 'Config')
            nucleus = getFieldSafe(twix_obj.hdr.Config, 'Nucleus', '');
        end
    catch
        nucleus = '';
    end
    
    if isempty(nucleus)
        warning('Could not determine nucleus from TWIX header. Using default 1H.');
        nucleus = '1H';
    end
    
    gammaTable = containers.Map({'1H','2H','3He','7Li','13C','15N','17O','19F','23Na','31P','129Xe'},...
                                {42.5760,6.5357,-32.4340,16.5470,10.7054,-4.3142,-5.7716,40.0541,11.2620,17.2350,-11.7770});
    
    if isKey(gammaTable, nucleus)
        gamma = gammaTable(nucleus);
    else
        warning('Unknown nucleus: %s. Using default gamma for 1H.', nucleus);
        gamma = 42.5760;
    end
end

function scanDate = findScanDate(twix_obj)
    try
        dateString = twix_obj.hdr.Dicom.tReferenceImage0;
        scanDate = datetime(dateString, 'InputFormat', 'yyyyMMdd');
    catch
        scanDate = datetime('today');
    end
end

function MRSIStruct = findAndSetFov(MRSIStruct, twix_obj)
    MRSIStruct.fov = struct();
    try
        if isfield(twix_obj.hdr.MeasYaps, 'sSliceArray') && isfield(twix_obj.hdr.MeasYaps.sSliceArray, 'asSlice')
            slice = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1};
            MRSIStruct.fov.x = slice.dReadoutFOV;
            MRSIStruct.fov.y = slice.dPhaseFOV;
            MRSIStruct.fov.z = slice.dThickness;
        elseif isfield(twix_obj.hdr, 'Config')
            MRSIStruct.fov.x = getFieldSafe(twix_obj.hdr.Config, 'ReadFoV', NaN);
            MRSIStruct.fov.y = getFieldSafe(twix_obj.hdr.Config, 'PhaseFoV', NaN);
            MRSIStruct.fov.z = getFieldSafe(twix_obj.hdr.Config, 'SliceThickness', NaN);
        else
            MRSIStruct.fov.x = NaN;
            MRSIStruct.fov.y = NaN;
            MRSIStruct.fov.z = NaN;
        end
        
        if isnan(MRSIStruct.fov.x) || isnan(MRSIStruct.fov.y) || isnan(MRSIStruct.fov.z) || ...
           MRSIStruct.fov.x <= 0 || MRSIStruct.fov.y <= 0 || MRSIStruct.fov.z <= 0
            warning('Could not extract FOV from TWIX header. Using default values.');
            if isnan(MRSIStruct.fov.x) || MRSIStruct.fov.x <= 0
                MRSIStruct.fov.x = 200;
            end
            if isnan(MRSIStruct.fov.y) || MRSIStruct.fov.y <= 0
                MRSIStruct.fov.y = 200;
            end
            if isnan(MRSIStruct.fov.z) || MRSIStruct.fov.z <= 0
                MRSIStruct.fov.z = 10;
            end
            fprintf('  Using FOV: %.1f x %.1f x %.1f mm\n', MRSIStruct.fov.x, MRSIStruct.fov.y, MRSIStruct.fov.z);
        end
    catch ME
        warning('Failed to extract FOV: %s. Using defaults.', ME.message);
        MRSIStruct.fov.x = 200;
        MRSIStruct.fov.y = 200;
        MRSIStruct.fov.z = 10;
    end
end

function MRSIStruct = calculateVoxelSize(MRSIStruct, numX, numY, numZ)
    MRSIStruct.voxelSize = struct();
    MRSIStruct.voxelSize.x = MRSIStruct.fov.x / numX;
    MRSIStruct.voxelSize.y = MRSIStruct.fov.y / numY;
    MRSIStruct.voxelSize.z = MRSIStruct.fov.z / numZ;
end

function MRSIStruct = findImageOrigin(MRSIStruct, twix_obj)
    MRSIStruct.imageOrigin = [0, 0, 0];
    try
        fields = {'VoI_Position_Sag', 'VoI_Position_Cor', 'VoI_Position_Tra'};
        for i = 1:3
            if isfield(twix_obj.hdr.Config, fields{i})
                MRSIStruct.imageOrigin(i) = twix_obj.hdr.Config.(fields{i});
            end
        end
    catch
    end
end

function MRSIStruct = calculateVoxelCoodinates(MRSIStruct)
    fovX = getFov(MRSIStruct, 'x');
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    xCoordinates = createCoordinates(fovX/2, voxSizeX);
    xCoordinates = xCoordinates - getImageOrigin(MRSIStruct, 'x');

    fovY = getFov(MRSIStruct, 'y');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    yCoordinates = createCoordinates(fovY/2, voxSizeY);
    yCoordinates = yCoordinates - getImageOrigin(MRSIStruct, 'y');

    MRSIStruct = setCoordinates(MRSIStruct, 'x', xCoordinates);
    MRSIStruct = setCoordinates(MRSIStruct, 'y', yCoordinates);
end

function MRSIStruct = calculateAffineMatrix(MRSIStruct, twix_obj)
    MRSIStruct.affineMatrix = eye(4);
    MRSIStruct.affineMatrix(1,1) = MRSIStruct.voxelSize.x;
    MRSIStruct.affineMatrix(2,2) = MRSIStruct.voxelSize.y;
    MRSIStruct.affineMatrix(3,3) = MRSIStruct.voxelSize.z;
end

function out = setDefaultFlagValues(out, isCartesian)
    out.flags.writtentostruct = 1;
    out.flags.gotparams = 1;
    out.flags.leftshifted = 0;
    out.flags.filtered = 0;
    out.flags.zeropadded = 0;
    out.flags.freqcorrected = 0;
    out.flags.phasecorrected = 0;
    out.flags.averaged = 0;
    out.flags.addedrcvrs = 0;
    out.flags.subtracted = 0;
    out.flags.writtentotext = 0;
    out.flags.downsampled = 0;
    out.flags.spatialFT = 0;
    out.flags.spectralFT = 0;
    out.flags.coilCombined = 0;
    out.flags.isFourSteps = 0;
    out.flags.isCartesian = isCartesian;
end

function value = getFov(MRSIStruct, dim)
    value = MRSIStruct.fov.(dim);
end

function value = getVoxSize(MRSIStruct, dim)
    value = MRSIStruct.voxelSize.(dim);
end

function value = getImageOrigin(MRSIStruct, dim)
    if strcmp(dim, 'x')
        value = MRSIStruct.imageOrigin(1);
    elseif strcmp(dim, 'y')
        value = MRSIStruct.imageOrigin(2);
    elseif strcmp(dim, 'z')
        value = MRSIStruct.imageOrigin(3);
    end
end

function MRSIStruct = setCoordinates(MRSIStruct, dim, coordinates)
    MRSIStruct.coordinates.(dim) = coordinates;
end

function coordinates = createCoordinates(halfFov, voxelSize)
    numVoxels = round(2 * halfFov / voxelSize);
    coordinates = linspace(-halfFov, halfFov, numVoxels);
end

function print_dims(dims, sz)
    fields = fieldnames(dims);
    for i = 1:length(fields)
        field = fields{i};
        idx = dims.(field);
        if idx > 0 && idx <= length(sz)
            fprintf('  dims.%s = %d (size: %d)\n', field, idx, sz(idx));
        elseif idx > 0
            fprintf('  dims.%s = %d (OUT OF BOUNDS)\n', field, idx);
        end
    end
end

function value = getFieldSafe(s, fieldName, defaultValue)
    if isfield(s, fieldName)
        value = s.(fieldName);
        if isempty(value)
            value = defaultValue;
        end
    else
        value = defaultValue;
    end
end

