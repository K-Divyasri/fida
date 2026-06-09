function [MRSIStruct, xShift_mm, yShift_mm] = io_MRSI_load_twix(filename, kFile)
% io_MRSI_load_twix - Universal Siemens TWIX loader for MRSI data
% Automatically detects sequence type and executes appropriate pipeline
%
% USAGE:
%   [MRSIStruct, xShift_mm, yShift_mm] = io_MRSI_load_twix(filename, kFile)
%
% INPUTS:
%   filename : Path to Siemens TWIX (.dat) file
%   kFile    : K-space trajectory file/data for ADC dwell time calculation
%              (optional for Cartesian, required for non-Cartesian)
%
% OUTPUTS:
%   MRSIStruct : Complete FID-A compatible structure with MRSI data
%   xShift_mm  : FOV shift in X direction (mm) - Sagittal
%   yShift_mm  : FOV shift in Y direction (mm) - Coronal
%
% SUPPORTED SEQUENCES:
%   - Concentric rings (detected by 'conc' or 'ring' in sequence name)
%   - Rosette (detected by 'ros' or 'selexc' in sequence name)
%   - Cartesian CSI (all other sequences)
%
% EXAMPLES:
%   % Concentric rings
%   [data, xShift, yShift] = io_MRSI_load_twix('scan_rings.dat', 'traj.mat');
%   
%   % Cartesian CSI
%   [data, xShift, yShift] = io_MRSI_load_twix('scan_csi.dat', '');
%
% Brenden Kadota, Jamie Near

%% Handle optional kFile argument
if nargin < 2
    kFile = '';
end

fprintf('\n========================================\n');
fprintf('  Universal MRSI TWIX Loader\n');
fprintf('========================================\n');
fprintf('File: %s\n', filename);
if ~isempty(kFile)
    fprintf('kFile: %s\n', mat2str(kFile));
end

try
    %% Step 1: Read TWIX file and detect sequence type
    twix_obj = readTwixFile(filename);
    sequence = twix_obj.hdr.Config.SequenceFileName;
    
    fprintf('\n--- Sequence Detection ---\n');
    fprintf('Sequence: %s\n', sequence);
    
    % Detect sequence type
    isConcentric = contains(sequence, 'conc', 'IgnoreCase', true) || ...
                   contains(sequence, 'ring', 'IgnoreCase', true);
    isRosette = contains(sequence, {'ros','selexc'}, 'IgnoreCase', true);
    isCartesian = ~isConcentric && ~isRosette;
    
    % Display detection results
    fprintf('Detected type: ');
    if isConcentric
        fprintf('CONCENTRIC RINGS\n');
        fprintf('Pipeline: io_concentric_load_twix\n');
    elseif isRosette
        fprintf('ROSETTE\n');
        fprintf('Pipeline: io_CSIload_twix3 (Rosette mode)\n');
    else
        fprintf('CARTESIAN CSI\n');
        fprintf('Pipeline: io_CSIload_twix3 (Cartesian mode)\n');
    end
    fprintf('-------------------------\n\n');
    
    %% Step 2: Execute appropriate pipeline
    if isConcentric
        % Use concentric rings pipeline
        fprintf('Executing CONCENTRIC RINGS pipeline...\n\n');
        MRSIStruct = loadConcentric(filename, kFile, twix_obj, sequence);
        
    else
        % Use CSI pipeline (handles both Cartesian and Rosette)
        fprintf('Executing CSI pipeline...\n\n');
        MRSIStruct = loadCSI(filename, kFile, twix_obj, sequence, isCartesian, isRosette);
    end
    
    %% Step 3: Compute FOV shifts (common for all pipelines)
    fprintf('\n--- Computing FOV Shifts ---\n');
    [xShift_mm, yShift_mm] = ComputeFOVShift(twix_obj);
    
    % Add shifts to structure
    MRSIStruct.xShift_mm = xShift_mm;
    MRSIStruct.yShift_mm = yShift_mm;
    
    fprintf('FOV Shifts: X=%.3f mm, Y=%.3f mm\n', xShift_mm, yShift_mm);
    fprintf('----------------------------\n');
    
    %% Step 4: Final summary
    fprintf('\n========================================\n');
    fprintf('  Loading Complete!\n');
    fprintf('========================================\n');
    fprintf('Data size: [%s]\n', num2str(size(MRSIStruct.data)));
    fprintf('Dwell time: %.3e s\n', MRSIStruct.adcDwellTime);
    fprintf('Spectral width: %.1f Hz\n', MRSIStruct.spectralWidth);
    fprintf('FOV: %.1f x %.1f x %.1f mm\n', ...
        MRSIStruct.fov.x, MRSIStruct.fov.y, MRSIStruct.fov.z);
    fprintf('FOV Shifts: X=%.3f mm, Y=%.3f mm\n', xShift_mm, yShift_mm);
    fprintf('========================================\n\n');
    
catch ME
    fprintf('\n========================================\n');
    fprintf('  ERROR DURING LOADING\n');
    fprintf('========================================\n');
    fprintf('Message: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    fprintf('========================================\n');
    rethrow(ME);
end

end

%% ========================================================================
%% PIPELINE FUNCTIONS
%% ========================================================================

function MRSIStruct = loadConcentric(filename, kFile, twix_obj, sequence)
% CONCENTRIC RINGS PIPELINE
% Data dimension order: [t, coils, averages, kshot, extras]

    fprintf('Loading concentric ring TWIX data (REORDERED)...\n');
    
    %% Extract raw data
    dOut.data = twix_obj.image();
    data = squeeze(dOut.data);
    
    fprintf('Raw data size: [%s]\n', num2str(size(data)));
    
    %% Get version and sequence info
    sqzDims = twix_obj.image.sqzDims;
    fprintf('Squeezed dimensions: {%s}\n', strjoin(sqzDims, ', '));
    
    isCartesian = false;
    isConcentric = true;
    spatialFT = false;
    spectralFT = false;
    
    %% Fill and permute dimensions - CONCENTRIC ORDER
    dims = fillDimsField(sqzDims, spatialFT, spectralFT, isCartesian, isConcentric);
    [dims, data] = permuteDims_REORDERED(dims, data, spatialFT, spectralFT, isCartesian, isConcentric);
    
    sz = size(data);
    fprintf('\n=== AFTER REORDERED PERMUTE ===\n');
    fprintf('Data size: [%s]\n', num2str(sz));
    fprintf('Target order: [t, coils, averages, kshot, extras]\n');
    fprintf('Dimension mapping:\n');
    print_dims(dims, sz);
    fprintf('====================\n\n');
    
    %% Calculate timing parameters using kFile
    [dwelltime, spectralwidth, leftshift] = calculateTimingParameters(twix_obj, dims, sz, kFile);
    
    n_time_points = sz(dims.t);
    adcTime = (0:(n_time_points-1)) * dwelltime;
    spectralTime = adcTime;
    spectralDwellTime = dwelltime;
    
    fprintf('Timing: dwell=%.3e s, SW=%.1f Hz, points=%d\n', ...
        dwelltime, spectralwidth, n_time_points);
    
    %% Determine matrix dimensions
    [numX, numY, numZ] = determineMatrixDimensions(twix_obj, dims, sz, isCartesian, isConcentric);
    
    fprintf('Matrix dimensions: X=%d, Y=%d, Z=%d\n', numX, numY, numZ);
    
    %% Determine subspecs and averages
    [subspecs, rawSubspecs, averages, rawAverages] = ...
        determineSubspecsAndAverages(dims, sz, twix_obj);
    
    fprintf('Averages: %d (raw: %d)\n', averages, rawAverages);
    fprintf('Subspecs: %d (raw: %d)\n', subspecs, rawSubspecs);
    
    %% Get nucleus and gamma
    [nucleus, gamma] = getNucleusAndGamma(twix_obj);
    
    %% Build the output structure
    MRSIStruct = struct();
    
    % Core data
    MRSIStruct.data = data;
    MRSIStruct.sz = sz;
    MRSIStruct.dims = dims;
    
    % Timing parameters
    MRSIStruct.adcDwellTime = dwelltime;
    MRSIStruct.adcTime = adcTime;
    MRSIStruct.spectralWidth = spectralwidth;
    MRSIStruct.spectralTime = spectralTime;
    MRSIStruct.spectralDwellTime = spectralDwellTime;
    
    % Scanner parameters
    MRSIStruct.txfrq = twix_obj.hdr.Meas.lFrequency;
    MRSIStruct.scanDate = findScanDate(twix_obj);
    MRSIStruct.Bo = twix_obj.hdr.Dicom.flMagneticFieldStrength;
    MRSIStruct.nucleus = nucleus;
    MRSIStruct.gamma = gamma;
    MRSIStruct.seq = sequence;
    MRSIStruct.te = twix_obj.hdr.MeasYaps.alTE{1}/1000;
    MRSIStruct.tr = twix_obj.hdr.MeasYaps.alTR{1}/1000;
    MRSIStruct.pointsToLeftshift = leftshift;
    
    % Spatial parameters
    MRSIStruct = findAndSetFov(MRSIStruct, twix_obj);
    MRSIStruct = calculateVoxelSize(MRSIStruct, numX, numY, numZ);
    
    % Acquisition parameters
    MRSIStruct.averages = averages;
    MRSIStruct.rawAverages = rawAverages;
    MRSIStruct.subspecs = subspecs;
    MRSIStruct.rawSubspecs = rawSubspecs;
    
    % Geometry
    MRSIStruct = findImageOrigin(MRSIStruct, twix_obj);
    MRSIStruct = calculateVoxelCoodinates(MRSIStruct);
    MRSIStruct = calculateAffineMatrix(MRSIStruct, twix_obj);
    
    % Flags
    MRSIStruct = setDefaultFlagValues(MRSIStruct, isCartesian);
    
    fprintf('Concentric data loading complete!\n');
end

function MRSIStruct = loadCSI(filename, kFile, twix_obj, sequence, isCartesian, isRosette)
% CSI PIPELINE (Cartesian and Rosette)
% Data dimension order: Standard CSI ordering

    fprintf('Loading CSI TWIX data...\n');
    
    %% Extract raw data
    rawData = squeeze(twix_obj.image());
    sqzDims = twix_obj.image.sqzDims;
    
    fprintf('Raw data size: [%s]\n', num2str(size(rawData)));
    fprintf('Squeezed dimensions: {%s}\n', strjoin(sqzDims, ', '));
    
    spatialFT = false;
    spectralFT = false;
    
    %% Build & permute DIMS
    dims = fillDimsField(sqzDims, spatialFT, spectralFT, isCartesian, false);
    [dims, rawData] = permuteDims(dims, rawData, spatialFT, spectralFT, isCartesian);
    
    sz = size(rawData);
    fprintf('\n=== AFTER PERMUTE ===\n');
    fprintf('Data size: [%s]\n', num2str(sz));
    print_dims(dims, sz);
    fprintf('====================\n\n');
    
    %% Determine dwell-time
    if (ischar(kFile) || isstring(kFile)) && strlength(kFile) > 0 && isfile(kFile)
        kdata = readmatrix(kFile);
        if size(kdata, 1) > 1 && size(kdata, 2) >= 4
            dwelltime = kdata(2,4) - kdata(1,4);
            dwellSource = 'k-file';
        else
            dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
            dwellSource = 'TWIX header';
        end
    else
        dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
        dwellSource = 'TWIX header';
    end
    
    assert(~isempty(dwelltime) && dwelltime > 0, ...
           'Unable to determine dwell-time – check kFile path or TWIX header.');
    
    try
        if dims.extras > 0
            adcTime = 0 : dwelltime : ((sz(dims.t)*sz(dims.extras))-1)*dwelltime;
        else
            adcTime = 0 : dwelltime : (sz(dims.t)-1)*dwelltime;
        end
    catch
        adcTime = 0 : dwelltime : (sz(dims.t)-1)*dwelltime;
    end
    spectralWidth = 1/dwelltime;
    
    fprintf('Timing: dwell=%.3e s (%s), SW=%.1f Hz\n', dwelltime, dwellSource, spectralWidth);
    
    %% Basic acquisition dimensions
    numX = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
    numY = twix_obj.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
    numZ = twix_obj.hdr.MeasYaps.sKSpace.dSliceResolution;
    
    fprintf('Matrix dimensions: X=%d, Y=%d, Z=%d\n', numX, numY, numZ);
    
    %% Averages / subspecs
    if dims.averages ~= 0
        averages = sz(dims.averages);
        rawAverages = averages;
    else
        averages = 1;
        rawAverages = 1;
    end
    if dims.subspec ~= 0
        subspecs = sz(dims.subspec);
        rawSubspecs = subspecs;
    else
        subspecs = 1;
        rawSubspecs = 1;
    end
    
    fprintf('Averages: %d (raw: %d)\n', averages, rawAverages);
    fprintf('Subspecs: %d (raw: %d)\n', subspecs, rawSubspecs);
    
    %% Populate FID-A struct
    MRSIStruct = struct();
    MRSIStruct.data = rawData;
    MRSIStruct.sz = size(rawData);
    MRSIStruct.dims = dims;
    
    % Timing
    MRSIStruct.adcDwellTime = dwelltime;
    MRSIStruct.adcTime = adcTime;
    MRSIStruct.spectralDwellTime = dwelltime;
    MRSIStruct.spectralWidth = spectralWidth;
    MRSIStruct.spectralTime = adcTime;
    
    % Scanner parameters
    MRSIStruct.txfrq = twix_obj.hdr.Meas.lFrequency;
    MRSIStruct.scanDate = findScanDate(twix_obj);
    MRSIStruct.Bo = twix_obj.hdr.Dicom.flMagneticFieldStrength;
    MRSIStruct.nucleus = '1H';
    MRSIStruct.gamma = 42.576;
    MRSIStruct.seq = sequence;
    MRSIStruct.te = twix_obj.hdr.MeasYaps.alTE{1}/1000;
    MRSIStruct.tr = twix_obj.hdr.MeasYaps.alTR{1}/1000;
    MRSIStruct.pointsToLeftshift = twix_obj.image.freeParam(1);
    
    % Spatial
    MRSIStruct = findAndSetFov(MRSIStruct, twix_obj);
    MRSIStruct = calculateVoxelSize(MRSIStruct, numX, numY, numZ);
    
    % Acquisition
    MRSIStruct.averages = averages;
    MRSIStruct.rawAverages = rawAverages;
    MRSIStruct.subspecs = subspecs;
    MRSIStruct.rawSubspecs = rawSubspecs;
    
    % Geometry
    MRSIStruct = findImageOrigin(MRSIStruct, twix_obj);
    MRSIStruct = calculateVoxelCoodinates(MRSIStruct);
    MRSIStruct = calculateAffineMatrix(MRSIStruct, twix_obj);
    
    % Flags
    MRSIStruct = setDefaultFlagValues(MRSIStruct, isCartesian);
    
    fprintf('CSI data loading complete!\n');
end

%% ========================================================================
%% HELPER FUNCTIONS - COMMON
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
    keys = {'t','f','coils','averages','timeinterleave','kx','ky','kz','x','y','z','kpts','kshot','subspec','extras'};
    for k = 1:numel(keys)
        dims.(keys{k}) = 0; 
    end

    if ~isSpatialFT
        if isConcentric
            % CONCENTRIC RINGS: Map raw dimensions
            map = containers.Map({'Col','Cha','Lin','Par','Rep','Ave','Seg','Eco'},...
                                 {'t','coils','kshot','extras','averages','averages','kshot','extras'});
            
            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isKey(map,lbl)
                    fld = map(lbl);
                    if dims.(fld) == 0
                        dims.(fld) = i;
                    else
                        if strcmp(fld,'kshot') || strcmp(fld,'extras') || strcmp(fld,'averages')
                            if dims.extras == 0
                                dims.extras = i;
                            elseif dims.timeinterleave == 0
                                dims.timeinterleave = i;
                            end
                        end
                    end
                else
                    if dims.extras == 0
                        dims.extras = i;
                    elseif dims.timeinterleave == 0
                        dims.timeinterleave = i;
                    end
                end
            end
            
        elseif isCartesian
            % Cartesian CSI
            map = containers.Map({'Col','Cha','Ave','Rep','Seg','Phs','Lin','Sli'},...
                                 {'t','coils','averages','timeinterleave','kx','ky','ky','kz'});
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
            
        else
            % Rosette/non-Cartesian
            fixed = struct('Col','t','Cha','coils','Ave','averages');
            nextIsKshot = true;

            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isfield(fixed,lbl)
                    dims.(fixed.(lbl)) = i;
                elseif strcmp(lbl,'Set')
                    dims.kshot = i;
                    nextIsKshot = false;
                elseif strcmp(lbl,'Eco')
                    dims.extras = i;
                else
                    if nextIsKshot
                        dims.kshot = i;
                        nextIsKshot = false;
                    else
                        if dims.extras == 0
                            dims.extras = i;
                        else
                            dims.timeinterleave = i;
                        end
                    end
                end
            end
        end
    else
        % Image space
        fixed = struct('Col','t','Cha','coils','Ave','averages','Lin','y','Par','z','Phs','x');
        for i = 1:numel(sqzDims)
            lbl = sqzDims{i};
            if isfield(fixed,lbl)
                dims.(fixed.(lbl)) = i;
            elseif dims.extras == 0
                dims.extras = i;
            else
                dims.timeinterleave = i;
            end
        end
        if isSpectralFT, dims.f = dims.t; dims.t = 0; end
    end
end

function [dims, data] = permuteDims_REORDERED(dims, data, isSpatialFT, isSpectralFT, isCartesian, isConcentric)
% PERMUTEDIMS_REORDERED - For concentric rings
% Re‑order to [t, coils, averages, kshot, extras]

    if isSpatialFT
        order = {'t','x','y','z','averages','coils','timeinterleave','extras'};
        if isSpectralFT, order{1} = 'f'; end
    elseif isCartesian
        order = {'t','ky','kx','kz','averages','coils','timeinterleave','extras'};
    else
        order = {'t','coils','averages','kshot','extras','timeinterleave'};
    end

    perm = [];
    for k = 1:numel(order)
        idx = dims.(order{k});
        if idx ~= 0, perm(end+1) = idx; end %#ok<AGROW>
    end
    perm = unique(perm,'stable');

    nd = ndims(data);
    if numel(perm) < nd
        perm = [perm, setdiff(1:nd, perm, 'stable')];
    end

    if ~isequal(perm, 1:numel(perm))
        data = permute(data, perm);
    end

    oldDims = dims;
    names = fieldnames(dims);
    for n = 1:numel(names), dims.(names{n}) = 0; end

    for newIdx = 1:numel(perm)
        oldIdx = perm(newIdx);
        for n = 1:numel(names)
            if oldDims.(names{n}) == oldIdx
                dims.(names{n}) = newIdx;
                break;
            end
        end
    end
end

function [dims, data] = permuteDims(dims, data, isSpatialFT, isSpectralFT, isCartesian)
% PERMUTEDIMS - Standard CSI permutation

    if isSpatialFT
        order = {'t','x','y','z','averages','coils','timeinterleave','extras'};
        if isSpectralFT, order{1} = 'f'; end
    elseif isCartesian
        order = {'t','ky','kx','kz','averages','coils','timeinterleave','extras'};
    else
        order = {'t','coils','averages','kpts','kshot','extras','timeinterleave'};
    end

    perm = [];
    for k = 1:numel(order)
        idx = dims.(order{k});
        if idx ~= 0, perm(end+1) = idx; end %#ok<AGROW>
    end
    perm = unique(perm,'stable');

    nd = ndims(data);
    if numel(perm) < nd
        perm = [perm, setdiff(1:nd, perm, 'stable')];
    end

    if ~isequal(perm, 1:numel(perm))
        data = permute(data, perm);
    end

    oldDims = dims;
    names = fieldnames(dims);
    for n = 1:numel(names), dims.(names{n}) = 0; end

    for newIdx = 1:numel(perm)
        oldIdx = perm(newIdx);
        for n = 1:numel(names)
            if oldDims.(names{n}) == oldIdx
                dims.(names{n}) = newIdx;
                break;
            end
        end
    end
end

function [dwelltime, spectralwidth, leftshift] = calculateTimingParameters(twix_obj, dims, sz, kFile)
    fprintf('\nCalculating timing parameters from kFile...\n');
    
    % Calculate ADC dwell time from kFile
    if isnumeric(kFile) && numel(kFile) > 1
        if size(kFile, 1) > 1
            if size(kFile, 2) >= 4
                dt = diff(kFile(:,4));
            else
                dt = diff(kFile(:,1));
            end
            dwelltime = mean(dt(dt > 0));
            fprintf('Dwell time calculated from kFile trajectory: %.3e s\n', dwelltime);
        else
            dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
            fprintf('Using dwell time from TWIX header: %.3e s\n', dwelltime);
        end
    elseif (ischar(kFile) || isstring(kFile)) && strlength(kFile) > 0
        try
            if isfile(kFile)
                try
                    kData = load(kFile);
                    if isstruct(kData)
                        fieldNames = fieldnames(kData);
                        if ismember('dwelltime', fieldNames)
                            dwelltime = kData.dwelltime;
                        elseif ismember('dt', fieldNames)
                            dwelltime = kData.dt;
                        elseif ismember('time', fieldNames)
                            dt = diff(kData.time);
                            dwelltime = mean(dt(dt > 0));
                        else
                            trajData = kData.(fieldNames{1});
                            if size(trajData, 2) >= 4
                                dt = diff(trajData(:,4));
                            else
                                dt = diff(trajData(:,1));
                            end
                            dwelltime = mean(dt(dt > 0));
                        end
                    else
                        if size(kData, 2) >= 4
                            dt = diff(kData(:,4));
                        else
                            dt = diff(kData(:,1));
                        end
                        dwelltime = mean(dt(dt > 0));
                    end
                catch
                    kData = readmatrix(kFile);
                    if size(kData, 2) >= 4
                        dt = diff(kData(:,4));
                    else
                        dt = diff(kData(:,1));
                    end
                    dwelltime = mean(dt(dt > 0));
                end
                fprintf('Dwell time loaded from kFile: %.3e s\n', dwelltime);
            else
                warning('kFile path does not exist, using TWIX header dwell time');
                dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
            end
        catch ME
            warning('Could not load kFile: %s. Using TWIX header dwell time', ME.message);
            dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
        end
    else
        dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;
        fprintf('Using dwell time from TWIX header: %.3e s\n', dwelltime);
    end
    
    spectralwidth = 1 / dwelltime;
    
    leftshift = 0;
    try
        if isfield(twix_obj.image, 'freeParam') && numel(twix_obj.image.freeParam) >= 1
            leftshift = twix_obj.image.freeParam(1);
        elseif isfield(twix_obj.hdr.MeasYaps, 'sSpecPara')
            if isfield(twix_obj.hdr.MeasYaps.sSpecPara, 'lAutoRefScanNo')
                leftshift = twix_obj.hdr.MeasYaps.sSpecPara.lAutoRefScanNo;
            end
        end
    catch
        leftshift = 0;
    end
    
    fprintf('Calculated timing: dwell=%.3e s, SW=%.1f Hz, leftshift=%d\n', ...
        dwelltime, spectralwidth, leftshift);
end

function [numX, numY, numZ] = determineMatrixDimensions(twix_obj, dims, sz, isCartesian, isConcentric)
    try
        numX = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution;
        numY = twix_obj.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
        numZ = twix_obj.hdr.MeasYaps.sKSpace.dSliceResolution;
    catch
        warning('Could not extract matrix dimensions from header. Using data dimensions.');
        if isCartesian
            numX = (dims.kx > 0) * sz(max(dims.kx, 1));
            numY = (dims.ky > 0) * sz(max(dims.ky, 1));
            numZ = (dims.kz > 0) * sz(max(dims.kz, 1));
        else
            numX = (dims.kshot > 0) * sz(max(dims.kshot, 1));
            numY = numX;
            numZ = 1;
        end
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

function [xShift_mm, yShift_mm] = ComputeFOVShift(twix_obj)
    fprintf('Computing FOV shifts...\n');
    
    twix = twix_obj;
    if iscell(twix)
        hdr = twix{2}.hdr;
    else
        hdr = twix.hdr;
    end

    normalStruct = hdr.MeasYaps.sSpecPara.sVoI.sNormal;
    nSag = 0; nCor = 0; nTra = 0;
    if isfield(normalStruct, 'dSag'); nSag = normalStruct.dSag; end
    if isfield(normalStruct, 'dCor'); nCor = normalStruct.dCor; end
    if isfield(normalStruct, 'dTra'); nTra = normalStruct.dTra; end
    normal = [nSag, nCor, nTra];

    [~, idx] = max(abs(normal));
    planes = {'Sagittal', 'Coronal', 'Axial (Transverse)'};
    closestPlane = planes{idx};

    switch closestPlane
        case 'Sagittal'
            rotation_angle_deg = acosd(abs(nSag));
            if abs(nCor) > abs(nTra)
                directionPlane = 'Sagittal to Coronal';
            else
                directionPlane = 'Sagittal to Axial';
            end
        case 'Coronal'
            rotation_angle_deg = acosd(abs(nCor));
            if abs(nSag) > abs(nTra)
                directionPlane = 'Coronal to Sagittal';
            else
                directionPlane = 'Coronal to Axial';
            end
        case 'Axial (Transverse)'
            rotation_angle_deg = acosd(abs(nTra));
            if abs(nCor) > abs(nSag)
                directionPlane = 'Axial to Coronal';
            else
                directionPlane = 'Axial to Sagittal';
            end
    end

    pSag = 0; pCor = 0; pTra = 0;
    if isfield(hdr.MeasYaps.sSpecPara.sVoI, 'sPosition')
        posStruct = hdr.MeasYaps.sSpecPara.sVoI.sPosition;
        if isfield(posStruct, 'dSag'); pSag = posStruct.dSag; end
        if isfield(posStruct, 'dCor'); pCor = posStruct.dCor; end
        if isfield(posStruct, 'dTra'); pTra = posStruct.dTra; end
        positionExists = true;
    else
        positionExists = false;
    end

    theta_rad = deg2rad(rotation_angle_deg);
    adjusted_Sagittal = 0;
    adjusted_Coronal = 0;

    if positionExists
        switch directionPlane
            case 'Axial to Coronal'
                adjusted_Coronal = pTra * sin(theta_rad);
                adjusted_Sagittal = pSag;
            case 'Axial to Sagittal'
                adjusted_Sagittal = pTra * sin(theta_rad);
                adjusted_Coronal = pCor;
            case 'Sagittal to Axial'
                adjusted_Sagittal = pSag * cos(theta_rad);
                adjusted_Coronal = pCor;
            case 'Sagittal to Coronal'
                adjusted_Sagittal = pSag * cos(theta_rad);
                adjusted_Coronal = pSag * sin(theta_rad);
            case 'Coronal to Axial'
                adjusted_Coronal = pCor * cos(theta_rad);
                adjusted_Sagittal = pSag;
            case 'Coronal to Sagittal'
                adjusted_Coronal = pCor * cos(theta_rad);
                adjusted_Sagittal = pCor * sin(theta_rad);
            otherwise
                adjusted_Sagittal = pSag;
                adjusted_Coronal = pCor;
        end
    end

    xShift_mm = adjusted_Sagittal;
    yShift_mm = adjusted_Coronal;
end

function [MRSIStruct] = findAndSetFov(MRSIStruct, twix_obj)
    MRSIStruct.fov = struct();
    try
        if isfield(twix_obj.hdr.MeasYaps, 'sSliceArray') && isfield(twix_obj.hdr.MeasYaps.sSliceArray, 'asSlice')
            slice = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1};
            MRSIStruct.fov.x = slice.dReadoutFOV;
            MRSIStruct.fov.y = slice.dPhaseFOV;
            MRSIStruct.fov.z = slice.dThickness;
        else
            MRSIStruct.fov.x = 200;
            MRSIStruct.fov.y = 200;
            MRSIStruct.fov.z = 10;
        end
    catch
        MRSIStruct.fov.x = 200;
        MRSIStruct.fov.y = 200;
        MRSIStruct.fov.z = 10;
    end
end

function scanDate = findScanDate(twix_obj)
    try
        dateString = twix_obj.hdr.MeasYaps.tReferenceImage0;
        scanDate = regexp(dateString, '\.(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})', 'names');
        scanDate = datetime(str2double(scanDate.year), str2double(scanDate.month), str2double(scanDate.day));
    catch
        try
            dateString = twix_obj.hdr.Dicom.tReferenceImage0;
            scanDate = datetime(dateString, 'InputFormat', 'yyyyMMdd');
        catch
            scanDate = datetime('today');
        end
    end
end

function MRSIStruct = calculateVoxelSize(MRSIStruct, numX, numY, numZ)
    MRSIStruct.voxelSize = struct();
    MRSIStruct.voxelSize.x = MRSIStruct.fov.x / numX;
    MRSIStruct.voxelSize.y = MRSIStruct.fov.y / numY;
    MRSIStruct.voxelSize.z = MRSIStruct.fov.z / numZ;
end

function MRSIStruct = findImageOrigin(MRSIStruct, twix_obj)
    MRSIStruct.imageOrigin = zeros(1,3);
    try
        if isfield(twix_obj.hdr.Config, 'VoI_Position_Sag')
            fields = {'VoI_Position_Sag', 'VoI_Position_Cor', 'VoI_Position_Tra'};
        elseif isfield(twix_obj.hdr.Config, 'Voi_Position_Sag')
            fields = {'Voi_Position_Sag', 'Voi_Position_Cor', 'Voi_Position_Tra'};
        else
            fields = {};
        end
        
        for i = 1:length(fields)
            if isfield(twix_obj.hdr.Config, fields{i}) && ~isempty(twix_obj.hdr.Config.(fields{i}))
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
    try
        [zVector, theta] = getZVectorAndTheta(twix_obj);
        rotationMatrix = getRotationMatrixFromVector(zVector(1), zVector(2), zVector(3), theta);
        rotationMatrix(4,4) = 1;

        yVector = [0, 1, -zVector(2)/zVector(3)];
        yVector = yVector/norm(yVector);
        yVector = rotationMatrix*[yVector'; 1];
        yVector = yVector(1:3);
        xVector = cross(yVector, zVector);

        affineRotationMatrix = [xVector', yVector, zVector', [0,0,0]'; 0 0 0 1];

        affineScaleMatrix = eye(4);
        affineScaleMatrix(1,1) = MRSIStruct.voxelSize.x;
        affineScaleMatrix(2,2) = MRSIStruct.voxelSize.y;
        affineScaleMatrix(3,3) = MRSIStruct.voxelSize.z;

        affineTranslateMatrix = eye(4);
        affineTranslateMatrix(1,4) = -MRSIStruct.imageOrigin(1) - MRSIStruct.fov.x/2;
        affineTranslateMatrix(2,4) = -MRSIStruct.imageOrigin(2) - MRSIStruct.fov.y/2;
        affineTranslateMatrix(3,4) = MRSIStruct.imageOrigin(3) - MRSIStruct.fov.z/2;

        affineMatrix = affineRotationMatrix*affineTranslateMatrix*affineScaleMatrix;
        MRSIStruct.affineMatrix = affineMatrix;
    catch
        MRSIStruct.affineMatrix = eye(4);
        MRSIStruct.affineMatrix(1,1) = MRSIStruct.voxelSize.x;
        MRSIStruct.affineMatrix(2,2) = MRSIStruct.voxelSize.y;
        MRSIStruct.affineMatrix(3,3) = MRSIStruct.voxelSize.z;
    end
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

function [z_vect, theta] = getZVectorAndTheta(twix_obj)
    z_vect = [0, 0, 1];
    theta = 0;
    try
        if isfield(twix_obj.hdr.Meas, 'VoI_Normal_Sag')
            z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Sag);
            z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Cor);
            z_vect(3) = initalizeOneIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Tra);
            theta = initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_InPlaneRotAngle);
        elseif isfield(twix_obj.hdr.Meas, 'VoiNormalSag')
            z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalSag);
            z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalCor);
            z_vect(3) = initalizeOneIfEmpty(twix_obj.hdr.Meas.VoiNormalTra);
            theta = initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiInPlaneRot);
        end
    catch
    end
end

function rotation_matrix = getRotationMatrixFromVector(x, y, z, theta)
    vect = [x,y,z];
    rotation_matrix = zeros(3,3);
    rotation_matrix(1,1) = cos(theta)+vect(1)^2*(1-cos(theta));
    rotation_matrix(1,2) = vect(1)*vect(2)*(1-cos(theta))-vect(3)*sin(theta);
    rotation_matrix(1,3) = vect(1)*vect(3)*(1-cos(theta))+vect(2)*sin(theta);
    rotation_matrix(2,1) = vect(2)*vect(1)*(1-cos(theta))+vect(3)*sin(theta);
    rotation_matrix(2,2) = cos(theta)+vect(2)^2*(1-cos(theta));
    rotation_matrix(2,3) = vect(2)*vect(3)*(1-cos(theta))-vect(1)*sin(theta);
    rotation_matrix(3,1) = vect(3)*vect(1)*(1-cos(theta))-vect(2)*sin(theta);
    rotation_matrix(3,2) = vect(3)*vect(2)*(1-cos(theta))+vect(2)*sin(theta);
    rotation_matrix(3,3) = cos(theta) + vect(3)^2*(1-cos(theta));
end

function out = initalizeZeroIfEmpty(value)
    if(isempty(value))
        out = 0;
    else
        out = value;
    end
end

function out = initalizeOneIfEmpty(value)
    if(isempty(value))
        out = 1;
    else
        out = value;
    end
end