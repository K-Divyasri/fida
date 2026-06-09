function [MRSIStruct, xShift_mm, yShift_mm] = io_CSIload_twix3(filename, kFile)
% io_CSIload_twix3  – Siemens-TWIX → FID-A loader (XA-safe)
%                      Cartesian + Non-Cartesian (Rosette / Concentric)
%
%   [out, xShift_mm, yShift_mm] = io_CSIload_twix3(datFile, kFile)
%
% Non-Cartesian layout (like your previous reader):
%   sz   = [Nt  Ncoils  Navg  Nkshot  Nextras]
%   dims = t:1, coils:2, averages:3, kshot:4, extras:5, kpts:0, ky:0, kx:0
%
% Cartesian layout:
%   dims = t, coils, averages, ky, kx, (kz), extras/timeinterleave as present
%
% Brenden Kadota & Jamie Near (Sunnybrook) — original
% Revised — XA-safe, rosette layout parity 

% -------------------------------------------------------------------------
arguments
    filename (1,:) {mustBeFile}
    kFile string   = ""
end
% -------------------------------------------------------------------------

%% 0) Load raw data / classify sequence -----------------------------------
twix_obj        = readTwixFile(filename);
rawData         = squeeze(twix_obj.image());      % remove singleton dims
sqzDims         = twix_obj.image.sqzDims;         % Siemens labels (cellstr)
sequence        = string(twix_obj.hdr.Config.SequenceFileName);

isRosette       = contains(sequence, {'ros','rosette'}, 'IgnoreCase', true);
isRings         = contains(sequence, {'conc','ring'},   'IgnoreCase', true);
isNonCart       = isRosette || isRings;
isCartesian     = ~isNonCart;

spatialFT       = false;  % RAW k-space
spectralFT      = false;  % RAW time-domain

%% 1) Build & permute DIMS (non-cart gets kshot/extras as in your rose) ---
dims                = buildDims(sqzDims, isCartesian);
[dims, rawData]     = applyCanonicalPermutation(dims, rawData, isCartesian);

%% 2) Dwell-time: kFile overrides if present ------------------------------
[dwelltime, dwellSource] = chooseDwellTime(twix_obj, kFile);

% adcTime length = Nt * Nextras (so op_CSICombineTime can fold it later)
Nt = size(rawData, dims.t);
Ne = 1; if dims.extras ~= 0, Ne = size(rawData, dims.extras); end
adcTime         = 0 : dwelltime : (Nt*Ne - 1)*dwelltime;
spectralWidth   = 1/dwelltime;

%% 3) Spatial matrix sizes (safe getters for XA) --------------------------
numX = getFieldOrDefault(twix_obj.hdr.MeasYaps.sKSpace, 'lBaseResolution',     1);
numY = getFieldOrDefault(twix_obj.hdr.MeasYaps.sKSpace, 'lPhaseEncodingLines', 1);
numZ = getFieldOrDefault(twix_obj.hdr.MeasYaps.sKSpace, 'dSliceResolution',    1);

%% 4) Averages / subspecs -------------------------------------------------
sz              = size(rawData);
if dims.averages ~= 0
    averages    = sz(dims.averages); rawAverages = averages;
else
    averages    = 1;                  rawAverages = 1;
end
if isfield(dims,'subspec') && dims.subspec ~= 0
    subspecs    = sz(dims.subspec);   rawSubspecs = subspecs;
else
    subspecs    = 1;                  rawSubspecs = 1;
end

%% 5) Populate FID-A struct ----------------------------------------------
MRSIStruct                      = struct();
MRSIStruct.data                 = rawData;
MRSIStruct.sz                   = size(rawData);
MRSIStruct.dims                 = dims;

% spectral/ADC timing (always present)
MRSIStruct.adcDwellTime         = dwelltime;
MRSIStruct.adcTime              = adcTime;
MRSIStruct.spectralDwellTime    = dwelltime;
MRSIStruct.spectralWidth        = spectralWidth;

% frequency / nucleus / meta
MRSIStruct.txfrq                = getTxFrequency(twix_obj);
MRSIStruct.scanDate             = findScanDate(twix_obj);
MRSIStruct.Bo                   = twix_obj.hdr.Dicom.flMagneticFieldStrength;
MRSIStruct.nucleus              = '1H';
MRSIStruct.gamma                = 42.576;                          % MHz/T
MRSIStruct.seq                  = char(sequence);
MRSIStruct.te                   = twix_obj.hdr.MeasYaps.alTE{1}/1000;
MRSIStruct.tr                   = twix_obj.hdr.MeasYaps.alTR{1}/1000;

% left-shift (TWIX freeParam is robust across VD/VE/XA)
if isprop(twix_obj.image, 'freeParam') && ~isempty(twix_obj.image.freeParam)
    MRSIStruct.pointsToLeftshift = twix_obj.image.freeParam(1);
else
    MRSIStruct.pointsToLeftshift = 0;
end

% FoV / geometry
MRSIStruct                      = findAndSetFov(MRSIStruct, twix_obj);
MRSIStruct                      = calcualteVoxelSize(MRSIStruct, numX, numY, numZ);
MRSIStruct                      = findImageOrigin(MRSIStruct, twix_obj);
MRSIStruct                      = calculateVoxelCoodinates(MRSIStruct);
MRSIStruct                      = calculateAffineMatrix(MRSIStruct, twix_obj);

% counters / flags
MRSIStruct.averages             = averages;
MRSIStruct.rawAverages          = rawAverages;
MRSIStruct.subspecs             = subspecs;
MRSIStruct.rawSubspecs          = rawSubspecs;
MRSIStruct                      = setDefaultFlagValues(MRSIStruct, isCartesian);

%% 6) VOI translation due to rotation (unchanged) -------------------------
[xShift_mm, yShift_mm]          = ComputeFOVShift(twix_obj);

fprintf('  dwell-time : %.3g µs   (%s)\n', dwelltime*1e6, dwellSource);
fprintf('  spectral-width : %.0f Hz\n',   spectralWidth);
end


% ========================================================================
%                            HELPERS (XA-SAFE)
% ========================================================================

function twix_obj = readTwixFile(filename)
    twix_obj = mapVBVD(char(filename));
    if iscell(twix_obj), twix_obj = twix_obj{end}; end
end

function [dwelltime, source] = chooseDwellTime(twix_obj, kFile)
    if strlength(kFile)>0 && isfile(kFile)
        kdata = readmatrix(kFile);
        % time column autodetect: prefer 4th, else 3rd, else error
        if size(kdata,2) >= 4
            dwelltime = kdata(3,4) - kdata(2,4);
        elseif size(kdata,2) >= 3
            dwelltime = kdata(3,3) - kdata(2,3);
        else
            error('kFile does not contain a recognizable time column.');
        end
        source = 'k-file';
    else
        dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9; % ns → s
        source = 'TWIX header';
    end
    assert(~isempty(dwelltime) && dwelltime>0, 'Unable to determine dwell-time.');
end

function val = getFieldOrDefault(S, fieldname, defaultVal)
    try
        if ~isempty(S) && isstruct(S) && isfield(S, fieldname)
            v = S.(fieldname);
            if ~isempty(v), val = v; return; end
        end
    catch
    end
    val = defaultVal;
end

function f = getTxFrequency(twix_obj)
    try
        f = twix_obj.hdr.Meas.lFrequency;
        if isempty(f), error('empty'); end
    catch
        Bo = twix_obj.hdr.Dicom.flMagneticFieldStrength;
        f  = Bo * 42.576 * 1e6; % Hz
    end
end

% ------------------------------------------------------------------------
% DIMS: make non-cart look like your old ROSETTE layout
% ------------------------------------------------------------------------
function dims = buildDims(sqzDims, isCartesian)
    keys = {'t','f','coils','averages','timeinterleave', ...
            'kx','ky','kz','x','y','z','kpts','kshot','subspec','extras'};
    for k = 1:numel(keys), dims.(keys{k}) = 0; end

    if isCartesian
        % Typical CSI Cartesian mapping (VD/VE/XA stable)
        mp = containers.Map( ...
            {'Col','Cha','Ave','Rep','Seg','Phs','Lin','Sli'}, ...
            {'t'  ,'coils','averages','timeinterleave', ...
             'kx' ,'ky'   ,'ky'      ,'kz'});
        for i = 1:numel(sqzDims)
            lbl = sqzDims{i};
            if isKey(mp,lbl)
                dims.(mp(lbl)) = i;
            elseif dims.extras == 0
                dims.extras = i;
            else
                dims.timeinterleave = i;
            end
        end
    else
        % === Non-Cartesian target: [t,coils,averages,kshot,extras] ===
        % Fixed parts common to all
        fixed = struct('Col','t','Cha','coils','Ave','averages');
        unknownCount = 0;

        for i = 1:numel(sqzDims)
            lbl = sqzDims{i};
            if isfield(fixed, lbl)
                dims.(fixed.(lbl)) = i;

            elseif strcmpi(lbl,'Set')
                % Siemens often uses Set for “shot”/“interleave”
                dims.kshot = i;

            else
                % First unknown → kshot (if not already assigned)
                % Second unknown → extras
                % Any others → timeinterleave (fallback)
                if dims.kshot == 0
                    dims.kshot = i;
                elseif dims.extras == 0
                    dims.extras = i;
                else
                    dims.timeinterleave = i;
                end
                unknownCount = unknownCount + 1;
            end
        end

        % Keep kpts unused (0) to mirror your previous data layout
        dims.kpts = 0;
        % No ky/kx in non-cart:
        dims.kx = 0; dims.ky = 0; dims.kz = 0;
    end
end

function [dims, data] = applyCanonicalPermutation(dims, data, isCartesian)
    if isCartesian
        order = {'t','coils','averages','ky','kx','kz','timeinterleave','extras'};
    else
        % Match your old rosette layout exactly:
        order = {'t','coils','averages','kshot','extras','timeinterleave'};
    end

    perm = [];
    for k = 1:numel(order)
        idx = dims.(order{k}); if idx ~= 0, perm(end+1) = idx; end %#ok<AGROW>
    end
    perm = unique(perm, 'stable');
    nd   = ndims(data);
    if numel(perm) < nd
        perm = [perm, setdiff(1:nd, perm, 'stable')];
    end

    if ~isequal(perm, 1:numel(perm))
        data = permute(data, perm);
    end

    % Rebuild dims to new indices
    old = dims; names = fieldnames(dims);
    for n = 1:numel(names), dims.(names{n}) = 0; end
    for newIdx = 1:numel(perm)
        oldIdx = perm(newIdx);
        for n = 1:numel(names)
            if old.(names{n}) == oldIdx
                dims.(names{n}) = newIdx; break;
            end
        end
    end
end

% ------------------------------------------------------------------------
% Geometry / meta (unchanged from your earlier versions)
% ------------------------------------------------------------------------
function [MRSIStruct] = findAndSetFov(MRSIStruct, twix_obj)
    fovX = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
    fovY = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
    fovZ = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    MRSIStruct.fov.x = fovX;
    MRSIStruct.fov.y = fovY;
    MRSIStruct.fov.z = fovZ;
end

function scanDate = findScanDate(twix_obj)
    scanDate = regexp(twix_obj.hdr.MeasYaps.tReferenceImage0,...
        '\.(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})', 'names');
    scanDate = datetime(str2double(scanDate.year), str2double(scanDate.month), str2double(scanDate.day));
end

function MRSIStruct = calcualteVoxelSize(MRSIStruct, numX, numY, numZ)
    MRSIStruct = setVoxelSize(MRSIStruct, 'x', getFov(MRSIStruct, 'x')/max(numX,1));
    MRSIStruct = setVoxelSize(MRSIStruct, 'y', getFov(MRSIStruct, 'y')/max(numY,1));
    MRSIStruct = setVoxelSize(MRSIStruct, 'z', getFov(MRSIStruct, 'z')/max(numZ,1));
end

function MRSIStruct = calculateVoxelCoodinates(MRSIStruct)
    fovX = getFov(MRSIStruct, 'x');
    fovY = getFov(MRSIStruct, 'y');
    voxX = getVoxSize(MRSIStruct, 'x');
    voxY = getVoxSize(MRSIStruct, 'y');

    xCoordinates = createCoordinates(fovX/2, voxX);
    yCoordinates = createCoordinates(fovY/2, voxY);

    xCoordinates = xCoordinates - getImageOrigin(MRSIStruct, 'x');
    yCoordinates = yCoordinates - getImageOrigin(MRSIStruct, 'y');

    MRSIStruct = setCoordinates(MRSIStruct, 'x', xCoordinates);
    MRSIStruct = setCoordinates(MRSIStruct, 'y', yCoordinates);
end

function MRSIStruct = calculateAffineMatrix(MRSIStruct, twixObj)
    [zVector, theta] = getZVectorAndTheta(twixObj);
    R = getRotationMatrixFromVector(zVector(1), zVector(2), zVector(3), theta);
    R(4,4) = 1;

    yVector = [0, 1, -zVector(2)/max(zVector(3),eps)];
    yVector = yVector / norm(yVector);
    yVector = (R * [yVector'; 1]); yVector = yVector(1:3);
    xVector = cross(yVector, zVector);

    Arot = [xVector', yVector, zVector', [0,0,0]'; 0 0 0 1];
    Ascl = eye(4);
    Ascl(1,1) = MRSIStruct.voxelSize.x;
    Ascl(2,2) = MRSIStruct.voxelSize.y;
    Ascl(3,3) = MRSIStruct.voxelSize.z;

    Atrs = eye(4);
    Atrs(1,4) = -MRSIStruct.imageOrigin(1) - MRSIStruct.fov.x/2;
    Atrs(2,4) = -MRSIStruct.imageOrigin(2) - MRSIStruct.fov.y/2;
    Atrs(3,4) =  MRSIStruct.imageOrigin(3) - MRSIStruct.fov.z/2;

    MRSIStruct.affineMatrix = Arot * Atrs * Ascl;
end

function [z_vect, theta] = getZVectorAndTheta(twix_obj)
    try
        z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalSag);
        z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalCor);
        z_vect(3) =  initalizeOneIfEmpty(twix_obj.hdr.Meas.VoiNormalTra);
        theta     =  initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiInPlaneRot);
    catch
        z_vect = [0 0 1]; theta = 0; % fallback
    end
end

function rotation_matrix = getRotationMatrixFromVector(x, y, z, theta)
    v = [x,y,z]; v = v / max(norm(v), eps);
    c = cos(theta); s = sin(theta); C = 1 - c;
    rotation_matrix = [ ...
        c+v(1)^2*C,       v(1)*v(2)*C - v(3)*s, v(1)*v(3)*C + v(2)*s; ...
        v(2)*v(1)*C + v(3)*s, c+v(2)^2*C,       v(2)*v(3)*C - v(1)*s; ...
        v(3)*v(1)*C - v(2)*s, v(3)*v(2)*C + v(1)*s, c+v(3)^2*C ];
end

function out = initalizeZeroIfEmpty(value)
    if isempty(value), out = 0; else, out = value; end
end

function out = initalizeOneIfEmpty(value)
    if isempty(value), out = 1; else, out = value; end
end

function MRSIStruct = findImageOrigin(MRSIStruct, twix_obj)
    MRSIStruct.imageOrigin = zeros(1,3);
    if isfield(twix_obj.hdr.Config, 'VoI_Position_Sag')
        fields = ["VoI_Position_Sag", "VoI_Position_Cor", "VoI_Position_Tra"];
    elseif isfield(twix_obj.hdr.Config, 'Voi_Position_Sag')
        fields = ["Voi_Position_Sag", "Voi_Position_Cor", "Voi_Position_Tra"];
    else
        fields = [];
    end
    for i = 1:length(fields)
        v = twix_obj.hdr.Config.(fields(i));
        if ~isempty(v), MRSIStruct.imageOrigin(i) = v; end
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

function val = getFov(S, axisName)
    val = S.fov.(axisName);
end

function S = setVoxelSize(S, axisName, val)
    if ~isfield(S,'voxelSize'), S.voxelSize = struct(); end
    S.voxelSize.(axisName) = val;
end

function val = getVoxSize(S, axisName)
    val = S.voxelSize.(axisName);
end

function coords = createCoordinates(halfFov, voxSize)
    if voxSize <= 0, coords = 0; return; end
    N = round((2*halfFov)/voxSize);
    start = -halfFov + voxSize/2;
    coords = start + (0:N-1)*voxSize;
end

function S = setCoordinates(S, axisName, vec)
    if ~isfield(S,'coordinates'), S.coordinates = struct(); end
    S.coordinates.(axisName) = vec(:).';
end

function val = getImageOrigin(S, axisName)
    switch axisName
        case 'x', val = S.imageOrigin(1);
        case 'y', val = S.imageOrigin(2);
        case 'z', val = S.imageOrigin(3);
        otherwise, val = 0;
    end
end

function [xShift_mm, yShift_mm] = ComputeFOVShift(twix_obj)
    if iscell(twix_obj)
        hdr = twix_obj{2}.hdr;
    else
        hdr = twix_obj.hdr;
    end

    nSag = 0; nCor = 0; nTra = 0;
    if isfield(hdr.MeasYaps,'sSpecPara') && isfield(hdr.MeasYaps.sSpecPara,'sVoI') ...
       && isfield(hdr.MeasYaps.sSpecPara.sVoI,'sNormal')
        ns = hdr.MeasYaps.sSpecPara.sVoI.sNormal;
        if isfield(ns,'dSag'), nSag = ns.dSag; end
        if isfield(ns,'dCor'), nCor = ns.dCor; end
        if isfield(ns,'dTra'), nTra = ns.dTra; end
    end
    normal = [nSag, nCor, nTra];
    [~, idx] = max(abs(normal));
    planes = {'Sagittal', 'Coronal', 'Axial'};
    closestPlane = planes{idx};
    switch closestPlane
        case 'Sagittal', rot_deg = acosd(abs(nSag));
        case 'Coronal',  rot_deg = acosd(abs(nCor));
        case 'Axial',    rot_deg = acosd(abs(nTra));
    end
    theta = deg2rad(rot_deg);

    pSag = 0; pCor = 0; pTra = 0;
    if isfield(hdr.MeasYaps,'sSpecPara') && isfield(hdr.MeasYaps.sSpecPara,'sVoI') ...
       && isfield(hdr.MeasYaps.sSpecPara.sVoI,'sPosition')
        ps = hdr.MeasYaps.sSpecPara.sVoI.sPosition;
        if isfield(ps,'dSag'), pSag = ps.dSag; end
        if isfield(ps,'dCor'), pCor = ps.dCor; end
        if isfield(ps,'dTra'), pTra = ps.dTra; end
    end

    switch closestPlane
        case 'Axial'
            xShift_mm = pTra * sin(theta);
            yShift_mm = pTra * sin(theta);
        case 'Sagittal'
            xShift_mm = pSag * cos(theta);
            yShift_mm = pSag * sin(theta);
        case 'Coronal'
            xShift_mm = pCor * sin(theta);
            yShift_mm = pCor * cos(theta);
    end
end