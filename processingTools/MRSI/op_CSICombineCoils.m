%======================================================================
% op_CSICombineCoils – combine receive-coils for CSI data
%
% [out, phaseMap, weightMap] = op_CSICombineCoils(in, ...
%                       samplePoint, phaseMap, weightMap, DEBUG)
%
% DEBUG == true prints intermediate shapes.
%======================================================================

function [MRSIStruct, phaseMap, weightMap] = op_CSICombineCoils( ...
            MRSIStruct, samplePoint, phaseMap, weightMap, DEBUG)

% ---------------- defaults -------------------------------------------
if nargin < 2 || isempty(samplePoint), samplePoint = 1;    end
if nargin < 3,                         phaseMap    = [];   end
if nargin < 4,                         weightMap   = [];   end
if nargin < 5,                         DEBUG       = false;end

assert(~isfield(MRSIStruct.flags,'addedrcvrs') || ...
       ~MRSIStruct.flags.addedrcvrs, ...
       'op_CSICombineCoils: data already coil-combined');
checkSpatialFT(MRSIStruct);

% ---------- reshape to canonical order  t,coils,y,x,...  --------------
[MRSIStruct, permInfo, origSize] = reshapeDimensions( ...
                               MRSIStruct, {'t','coils','y','x'});
data    = getData(MRSIStruct);
coilDim = getDimension(MRSIStruct,'coils');
yDim    = getDimension(MRSIStruct,'y');
xDim    = getDimension(MRSIStruct,'x');

if DEBUG
    fprintf('\n[CoilCombine DEBUG]\n');
    fprintf('  size after reshapeDimensions : %s\n', mat2str(size(data)));
    fprintf('  coilDim=%d  yDim=%d  xDim=%d  samplePoint=%d\n', ...
            coilDim,yDim,xDim,samplePoint);
end

% ====================  (1) phase map  ================================
if isempty(phaseMap)
    phaseMap = squeeze( angle( data(samplePoint,:,:,:,:) ) );  % (Coils,Ny,Nx)
    phaseMap = reshape(phaseMap,[size(data,coilDim) size(data,yDim) size(data,xDim)]);
end
if DEBUG, fprintf('  phaseMap size               : %s\n', mat2str(size(phaseMap))); end
data = applySpatialMap(data, exp(-1i*phaseMap), coilDim, yDim, xDim);

% ====================  (2) weight map ================================
if isempty(weightMap)
    mag = squeeze( abs( data(samplePoint,:,:,:,:) ) );
    mag = reshape(mag, size(phaseMap));
    weightMap = mag ./ sqrt(sum(mag.^2,1,'omitnan'));
end
if DEBUG, fprintf('  weightMap size              : %s\n', mat2str(size(weightMap))); end
data = applySpatialMap(data, weightMap, coilDim, yDim, xDim);

% put weighted, phase-corrected data back into struct
MRSIStruct = setData(MRSIStruct,data);

% ---------- reshape **back** BEFORE summing over coils ---------------
MRSIStruct = reshapeBack(MRSIStruct, permInfo, origSize);

% ---------- now safely collapse the coil dimension -------------------
coilDimOrig = getDimension(MRSIStruct,'coils');
data        = sum(getData(MRSIStruct), coilDimOrig);   % size==1 along coil axis
MRSIStruct  = setData(MRSIStruct,data);
MRSIStruct  = removeDimension(MRSIStruct,'coils');
MRSIStruct  = setFlags(MRSIStruct,'addedrcvrs',true);

if DEBUG
    fprintf('  final data size             : %s\n', mat2str(size(MRSIStruct.data)));
    fprintf('[CoilCombine DEBUG] done.\n\n');
end
end
%======================================================================


%----------------------------------------------------------------------
% applySpatialMap  – broadcast (Coils×Ny×Nx) map onto N-D data
%----------------------------------------------------------------------
function out = applySpatialMap(data, map3D, coilDim, yDim, xDim)
nd   = ndims(data);

shape           = ones(1,nd);
shape(coilDim)  = size(map3D,1);
shape(yDim)     = size(map3D,2);
shape(xDim)     = size(map3D,3);

mapND = reshape(map3D, shape);   % singleton-expansion ready
out   = data .* mapND;           % implicit broadcasting
end
