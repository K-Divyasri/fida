%======================================================================
% op_CSICombineCoils1 – combine receive-coils for CSI data
% Enhanced version that works with both Cartesian and non-Cartesian data
%
% [out, phaseMap, weightMap] = op_CSICombineCoils1(in, ...
% samplePoint, phaseMap, weightMap, DEBUG)
%
% DEBUG == true prints intermediate shapes.
%======================================================================
function [MRSIStruct, phaseMap, weightMap] = op_CSICombineCoils1( ...
 MRSIStruct, samplePoint, phaseMap, weightMap, DEBUG)

% ---------------- defaults -------------------------------------------
if nargin < 2 || isempty(samplePoint), samplePoint = 1; end
if nargin < 3, phaseMap = []; end
if nargin < 4, weightMap = []; end
if nargin < 5, DEBUG = true; end % Enable debug by default for troubleshooting

assert(~isfield(MRSIStruct.flags,'addedrcvrs') || ...
 ~MRSIStruct.flags.addedrcvrs, ...
'op_CSICombineCoils1: data already coil-combined');

checkSpatialFT(MRSIStruct);

% ---------- Detect data type based on actual dimensions ---------------
dims = MRSIStruct.dims;
originalSize = MRSIStruct.sz;

% Check if averages dimension exists AND has size > 1
hasAverages = (dims.averages > 0) && (dims.averages <= length(originalSize)) && (originalSize(dims.averages) > 1);

if DEBUG
    fprintf('\n[CoilCombine1 DEBUG] Initial Analysis:\n');
    fprintf(' Input data size: %s\n', mat2str(originalSize));
    fprintf(' Averages dim index: %d\n', dims.averages);
    if hasAverages
        fprintf(' Averages size: %d (HAS AVERAGES)\n', originalSize(dims.averages));
    else
        fprintf(' No valid averages dimension (NO AVERAGES)\n');
    end
end

% Choose appropriate dimension ordering
if hasAverages
    % For non-Cartesian data with averages: include averages in reshape
    targetDims = {'t','coils','averages','y','x'};
    fprintf('Mode: Non-Cartesian with averages\n');
else
    % For standard Cartesian data or non-Cartesian without averages
    targetDims = {'t','coils','y','x'};
    fprintf('Mode: Cartesian or non-Cartesian without averages\n');
end

% ---------- reshape to canonical order --------------------------------
[MRSIStruct, permInfo, origSize] = reshapeDimensions(MRSIStruct, targetDims);
data = getData(MRSIStruct);

% Get dimension indices after reshaping
coilDim = getDimension(MRSIStruct,'coils');
yDim = getDimension(MRSIStruct,'y');
xDim = getDimension(MRSIStruct,'x');

% For non-Cartesian with averages, also get averages dimension
if hasAverages
    avgDim = getDimension(MRSIStruct,'averages');
else
    avgDim = 0;
end

if DEBUG
 fprintf(' Size after reshapeDimensions: %s\n', mat2str(size(data)));
 fprintf(' coilDim=%d yDim=%d xDim=%d', coilDim, yDim, xDim);
 if hasAverages
     fprintf(' avgDim=%d', avgDim);
 end
 fprintf(' samplePoint=%d\n', samplePoint);
end

% ==================== (1) phase map ================================
if isempty(phaseMap)
 % Extract reference signal for phase correction
 if hasAverages
     % For non-Cartesian with averages: average across averages dimension first
     refSignal = squeeze(data(samplePoint,:,:,:,:)); % (Coils,Avg,Ny,Nx)
     refSignal = squeeze(mean(refSignal, avgDim-1)); % Average over averages, (Coils,Ny,Nx)
     phaseMap = angle(refSignal);
 else
     % For Cartesian or non-Cartesian without averages: direct extraction
     refData = data(samplePoint,:,:,:);
     if ndims(refData) == 4
         phaseMap = squeeze(angle(refData)); % Remove singleton time dimension
     else
         phaseMap = angle(refData);
     end
     phaseMap = reshape(phaseMap,[size(data,coilDim) size(data,yDim) size(data,xDim)]);
 end
end

if DEBUG, fprintf(' phaseMap size: %s\n', mat2str(size(phaseMap))); end

% Apply phase correction
if hasAverages
    data = applySpatialMapWithAverages(data, exp(-1i*phaseMap), coilDim, avgDim, yDim, xDim);
else
    data = applySpatialMap(data, exp(-1i*phaseMap), coilDim, yDim, xDim);
end

% ==================== (2) weight map ================================
if isempty(weightMap)
 if hasAverages
     % For non-Cartesian with averages: get magnitude and average across averages
     mag = squeeze(abs(data(samplePoint,:,:,:,:))); % (Coils,Avg,Ny,Nx)
     mag = squeeze(mean(mag, avgDim-1)); % Average over averages, (Coils,Ny,Nx)
     weightMap = mag ./ sqrt(sum(mag.^2,1,'omitnan'));
 else
     % For Cartesian or non-Cartesian without averages: direct calculation
     magData = abs(data(samplePoint,:,:,:));
     if ndims(magData) == 4
         mag = squeeze(magData); % Remove singleton time dimension
     else
         mag = magData;
     end
     mag = reshape(mag, size(phaseMap));
     weightMap = mag ./ sqrt(sum(mag.^2,1,'omitnan'));
 end
end

if DEBUG, fprintf(' weightMap size: %s\n', mat2str(size(weightMap))); end

% Apply weight correction
if hasAverages
    data = applySpatialMapWithAverages(data, weightMap, coilDim, avgDim, yDim, xDim);
else
    data = applySpatialMap(data, weightMap, coilDim, yDim, xDim);
end

% put weighted, phase-corrected data back into struct
MRSIStruct = setData(MRSIStruct,data);

% ---------- reshape **back** BEFORE summing over coils ---------------
MRSIStruct = reshapeBack(MRSIStruct, permInfo, origSize);

% ---------- now safely collapse the coil dimension -------------------
coilDimOrig = getDimension(MRSIStruct,'coils');
data = sum(getData(MRSIStruct), coilDimOrig); % size==1 along coil axis

% **NEW**: Squeeze out the singular coil dimension
data = squeeze(data);
if DEBUG
    fprintf(' data size after sum: %s\n', mat2str(size(getData(MRSIStruct))));
    fprintf(' data size after squeeze: %s\n', mat2str(size(data)));
end

MRSIStruct = setData(MRSIStruct,data);
MRSIStruct = removeDimension(MRSIStruct,'coils');

% **NEW**: Update the size array to reflect the squeezed dimensions
MRSIStruct.sz = size(data);

MRSIStruct = setFlags(MRSIStruct,'addedrcvrs',true);

if DEBUG
 fprintf(' final data size: %s\n', mat2str(size(getData(MRSIStruct))));
 fprintf(' final MRSIStruct.sz: %s\n', mat2str(MRSIStruct.sz));
 fprintf(' final dims after coil removal:\n');
 disp(MRSIStruct.dims);
 fprintf('[CoilCombine1 DEBUG] done.\n\n');
end

end

%======================================================================
%----------------------------------------------------------------------
% applySpatialMap – broadcast (Coils×Ny×Nx) map onto N-D data
% (Original function for Cartesian data)
%----------------------------------------------------------------------
function out = applySpatialMap(data, map3D, coilDim, yDim, xDim)
nd = ndims(data);
shape = ones(1,nd);
shape(coilDim) = size(map3D,1);
shape(yDim) = size(map3D,2);
shape(xDim) = size(map3D,3);
mapND = reshape(map3D, shape); % singleton-expansion ready
out = data .* mapND; % implicit broadcasting
end

%----------------------------------------------------------------------
% applySpatialMapWithAverages – enhanced version for non-Cartesian data
% Handles averages dimension by broadcasting the map across it
%----------------------------------------------------------------------
function out = applySpatialMapWithAverages(data, map3D, coilDim, avgDim, yDim, xDim)
nd = ndims(data);
shape = ones(1,nd);
shape(coilDim) = size(map3D,1);
shape(yDim) = size(map3D,2);
shape(xDim) = size(map3D,3);
% The averages dimension should be broadcast (size 1 in the map)
mapND = reshape(map3D, shape); % singleton-expansion ready
out = data .* mapND; % implicit broadcasting across all dimensions including averages
end

%======================================================================
% op_CSIAverage1.m - Enhanced version that works with dimension tracking
%
% USAGE:
% [out]=op_CSIAverage1(in);
% 
% DESCRIPTION:
% Averages CSI FID-A structure along average dimension with better error handling.
% Works correctly after coil combination for non-Cartesian data.
% 
% INPUTS:
% in   = Input CSI FID-A structure
%
% OUTPUTS:
% out = Averaged CSI FID-A structure
%======================================================================

function MRSIStruct = op_CSIAverage1(MRSIStruct)
    flag = checkArgumentsEnhanced(MRSIStruct);
    %if check arguments returns true, abort
    if(flag); return; end

    fprintf('Starting averaging operation...\n');
    fprintf('Input data size: %s\n', mat2str(MRSIStruct.sz));
    fprintf('Input dims:\n');
    disp(MRSIStruct.dims);

    data = getData(MRSIStruct);

    averageDimension = getDimension(MRSIStruct, 'averages');
    fprintf('Averaging along dimension %d (size=%d)\n', averageDimension, MRSIStruct.sz(averageDimension));
    
    data = squeeze(mean(data, averageDimension));
    fprintf('Data size after averaging: %s\n', mat2str(size(data)));

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = setFlags(MRSIStruct, 'averaged',  1);
    MRSIStruct = removeDimension(MRSIStruct, 'averages');
    
    fprintf('Final data size: %s\n', mat2str(MRSIStruct.sz));
    fprintf('Final dims:\n');
    disp(MRSIStruct.dims);
    fprintf('Averaging completed successfully.\n\n');
end

% Enhanced check function with better error reporting
function flag = checkArgumentsEnhanced(in)
    flag = false;
    
    % Check if averages dimension exists and is valid
    if ~isfield(in.dims, 'averages') || in.dims.averages == 0
        fprintf('\n=== AVERAGING SKIPPED ===\n');
        fprintf('No averages dimension found in data structure.\n');
        fprintf('This is normal for water-suppressed data or already-averaged data.\n');
        fprintf('Current dims structure:\n');
        disp(in.dims);
        fprintf('Data size: %s\n', mat2str(in.sz));
        fprintf('========================\n\n');
        flag = true;
        return;
    end
    
    % Check if already averaged
    if isfield(in.flags, 'averaged') && in.flags.averaged
        fprintf('Data already averaged! Skipping averaging operation.\n');
        flag = true;
        return;
    end
    
    % Check if averages dimension has valid size
    avgDim = in.dims.averages;
    if avgDim > length(in.sz)
        fprintf('ERROR: Averages dimension index (%d) exceeds data dimensions (%d)\n', ...
            avgDim, length(in.sz));
        flag = true;
        return;
    end
    
    if in.sz(avgDim) <= 1
        fprintf('WARNING: Averages dimension has size %d (≤1), cannot average.\n', in.sz(avgDim));
        fprintf('This may indicate the data was already averaged or has no repetitions.\n');
        flag = true;
        return;
    end
    
    fprintf('Averages dimension validation passed: dim=%d, size=%d\n', avgDim, in.sz(avgDim));
end