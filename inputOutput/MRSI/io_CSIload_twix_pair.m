function [met, ref] = io_CSIload_twix_pair(metFile, refFile, kFile)
% IO_CSILOAD_TWIX_PAIR  Load a metabolite + water-reference TWIX pair for
%                       MRSI (rosette / concentric-rings / Cartesian).
%
% Replaces these three-stage chains:
%   --- ROSETTE ---
%       [a, b] = load_twix2(metFile, refFile, kFile);    % io_CSIload_twix3 + combine + shift
%       met    = reshape_twix_data(a, kFile);
%       ref    = reshape_twix_data(b, kFile);
%
%   --- CONCENTRIC RINGS ---
%       cc       = io_MRSI_load_twix(metFile, kFile);
%       cc_w     = io_MRSI_load_twix(refFile, kFile);
%       a        = op_CSICombineTime(cc,   'extras');
%       b        = op_CSICombineTime(cc_w, 'extras');
%       met      = reshape_twix_data(a, kFile);
%       ref      = reshape_twix_data(b, kFile);
%
% USAGE:
%   [met, ref] = io_CSIload_twix_pair(metFile, refFile, kFile)
%
% INPUTS:
%   metFile  – Siemens TWIX (.dat), water-suppressed acquisition
%   refFile  – Siemens TWIX (.dat), water-reference acquisition
%   kFile    – CSV trajectory file ("TR,Kx,Ky,time"); pass '' for Cartesian
%
% OUTPUTS:
%   met, ref – FID-A MRSI structs.  For non-Cartesian:
%                  data = [t, coils, (avg), kpts, kshot]
%              VOI offset is applied for ROSETTE only (matches the previous
%              behaviour of load_twix2; concentric pipeline never applied it).
%
% Sequence detection is automatic, based on
% twix_obj.hdr.Config.SequenceFileName:
%   'ros' / 'rosette' / 'selexc'  → rosette
%   'conc' / 'ring'               → concentric rings
%   anything else                 → Cartesian CSI

arguments
    metFile (1,:) char {mustBeFile}
    refFile (1,:) char {mustBeFile}
    kFile   (1,:) char = ''
end

% Trajectory points-per-cycle (only needed for non-Cartesian).
kPtsPerCycle = [];
if ~isempty(kFile) && isfile(kFile)
    [kTable, ~]  = readKFile(kFile);
    kPtsPerCycle = getKPtsPerCycle(kTable);
end

fprintf('=== Loading metabolite: %s ===\n', metFile);
[met, x_mm, y_mm, seqType] = loadOneTwix(metFile, kFile);

fprintf('=== Loading water-ref:  %s ===\n', refFile);
[ref, ~, ~, ~]             = loadOneTwix(refFile, kFile);   % VOI offset reused from met

met = postProcess(met, seqType, x_mm, y_mm, kFile, kPtsPerCycle);
ref = postProcess(ref, seqType, x_mm, y_mm, kFile, kPtsPerCycle);
end


% =========================================================================
%   POST-PROCESS  (per-sequence)
% =========================================================================
function s = postProcess(s, seqType, x_mm, y_mm, kFile, kPtsPerCycle)
    switch seqType
        case 'cartesian'
            % Loader output is the final result.
            return

        case 'rosette'
            s = op_CSICombineTime1(s, 'extras');
            s = op_CSIShift(s, x_mm, y_mm, kFile);    % VOI offset (rosette only)
            s = splitReadoutKpts(s, kPtsPerCycle);

        case 'concentric'
            s = op_CSICombineTime1(s, 'extras');
            s = splitReadoutKpts(s, kPtsPerCycle);
    end
end


% =========================================================================
%   SPLIT READOUT  →  [t, kpts]            (was: reshape_twix_data)
%
%   In:   4D [t*kpts, coils, avg,  kshot]   →   5D [t, coils, avg, kpts, kshot]
%         3D [t*kpts, coils,       kshot]   →   4D [t, coils,      kpts, kshot]
% =========================================================================
function s = splitReadoutKpts(s, kPtsPerCycle)
    if isempty(kPtsPerCycle), return; end

    sz   = s.sz;
    data = s.data;
    kp   = kPtsPerCycle;

    if numel(sz) == 4
        nT   = sz(1) / kp;
        data = reshape(data, [kp, nT, sz(2), sz(3), sz(4)]);
        data = permute(data, [2, 3, 4, 1, 5]);
        s.dims.kpts  = 4;
        s.dims.kshot = 5;

    elseif numel(sz) == 3
        nT   = sz(1) / kp;
        data = reshape(data, [kp, nT, sz(2), sz(3)]);
        data = permute(data, [2, 3, 1, 4]);
        s.dims.kpts  = 3;
        s.dims.kshot = 4;
    end

    s.data = data;
    s.sz   = size(data);
end


% =========================================================================
%   LOAD ONE TWIX FILE
%   (was: io_CSIload_twix3 / loadCSI / loadConcentric in io_MRSI_load_twix)
% =========================================================================
function [s, xShift_mm, yShift_mm, seqType] = loadOneTwix(filename, kFile)
    twix_obj = readTwixFile(filename);
    rawData  = squeeze(twix_obj.image());
    sqzDims  = twix_obj.image.sqzDims;
    sequence = string(twix_obj.hdr.Config.SequenceFileName);

    seqType  = detectSequence(sequence);
    fprintf('  sequence: %s  →  type: %s\n', sequence, seqType);

    % Build & permute dims into canonical layout.
    dims              = buildDims(sqzDims, seqType);
    [dims, rawData]   = applyCanonicalPermutation(dims, rawData, seqType);

    % Dwell time — k-file (column 4 = time) overrides the TWIX header.
    [dwelltime, src]  = chooseDwellTime(twix_obj, kFile);
    Nt = size(rawData, dims.t);
    Ne = 1;  if dims.extras ~= 0, Ne = size(rawData, dims.extras); end
    adcTime         = 0 : dwelltime : (Nt*Ne - 1)*dwelltime;
    spectralWidth   = 1/dwelltime;

    % Spatial matrix sizes.
    numX = getFieldOrDefault(twix_obj.hdr.MeasYaps.sKSpace, 'lBaseResolution',     1);
    numY = getFieldOrDefault(twix_obj.hdr.MeasYaps.sKSpace, 'lPhaseEncodingLines', 1);
    numZ = getFieldOrDefault(twix_obj.hdr.MeasYaps.sKSpace, 'dSliceResolution',    1);

    % Averages / subspecs.
    sz = size(rawData);
    if dims.averages ~= 0,                  averages = sz(dims.averages); else, averages = 1; end
    if isfield(dims,'subSpecs') && dims.subSpecs ~= 0
        subspecs = sz(dims.subSpecs);
    else
        subspecs = 1;
    end

    % Nucleus / gamma.
    [nucleus, gamma] = getNucleusAndGamma(twix_obj);

    % Populate FID-A struct.
    s = struct();
    s.data              = rawData;
    s.sz                = size(rawData);
    s.dims              = dims;

    s.adcDwellTime      = dwelltime;
    s.adcTime           = adcTime;
    s.spectralDwellTime = dwelltime;
    s.spectralWidth     = spectralWidth;
    s.spectralTime      = adcTime;

    s.txfrq             = getTxFrequency(twix_obj);
    s.scanDate          = findScanDate(twix_obj);
    s.Bo                = twix_obj.hdr.Dicom.flMagneticFieldStrength;
    s.nucleus           = nucleus;
    s.gamma             = gamma;
    s.seq               = char(sequence);
    s.te                = twix_obj.hdr.MeasYaps.alTE{1}/1000;
    s.tr                = twix_obj.hdr.MeasYaps.alTR{1}/1000;

    if isprop(twix_obj.image, 'freeParam') && ~isempty(twix_obj.image.freeParam)
        s.pointsToLeftshift = twix_obj.image.freeParam(1);
    else
        s.pointsToLeftshift = 0;
    end

    s = findAndSetFov(s, twix_obj);
    s = calculateVoxelSize(s, numX, numY, numZ);
    s = findImageOrigin(s, twix_obj);
    s = calculateVoxelCoordinates(s);
    s = calculateAffineMatrix(s, twix_obj);

    s.averages    = averages;
    s.rawAverages = averages;
    s.subspecs    = subspecs;
    s.rawSubspecs = subspecs;
    s             = setDefaultFlagValues(s, strcmp(seqType,'cartesian'));

    [xShift_mm, yShift_mm] = computeFOVShift(twix_obj);
    s.xShift_mm = xShift_mm;
    s.yShift_mm = yShift_mm;

    fprintf('  dwell-time : %.3g µs   (%s)\n', dwelltime*1e6, src);
    fprintf('  spectral-width : %.0f Hz\n',   spectralWidth);
end


% =========================================================================
%   SEQUENCE DETECTION
% =========================================================================
function seqType = detectSequence(sequence)
    if contains(sequence, {'conc','ring'},          'IgnoreCase', true)
        seqType = 'concentric';
    elseif contains(sequence, {'ros','rosette','selexc'}, 'IgnoreCase', true)
        seqType = 'rosette';
    else
        seqType = 'cartesian';
    end
end


% =========================================================================
%   TWIX HELPERS
% =========================================================================
function twix_obj = readTwixFile(filename)
    twix_obj = mapVBVD(char(filename));
    if iscell(twix_obj), twix_obj = twix_obj{end}; end
end

function [dwelltime, source] = chooseDwellTime(twix_obj, kFile)
    if ~isempty(kFile) && isfile(kFile)
        kdata = readmatrix(kFile);
        if size(kdata,2) >= 4
            dwelltime = kdata(3,4) - kdata(2,4);
        elseif size(kdata,2) >= 3
            dwelltime = kdata(3,3) - kdata(2,3);
        else
            error('kFile does not contain a recognisable time column.');
        end
        source = 'k-file';
    else
        dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9;   % ns → s
        source = 'TWIX header';
    end
    assert(~isempty(dwelltime) && dwelltime>0, 'Unable to determine dwell-time.');
end

function val = getFieldOrDefault(S, fieldname, defaultVal)
    val = defaultVal;
    try
        if ~isempty(S) && isstruct(S) && isfield(S, fieldname)
            v = S.(fieldname);
            if ~isempty(v), val = v; end
        end
    catch
    end
end

function f = getTxFrequency(twix_obj)
    try
        f = twix_obj.hdr.Meas.lFrequency;
        if isempty(f), error('empty'); end
    catch
        Bo = twix_obj.hdr.Dicom.flMagneticFieldStrength;
        f  = Bo * 42.576 * 1e6;   % Hz
    end
end

function [nucleus, gamma] = getNucleusAndGamma(twix_obj)
    nucleus = '';
    try
        if isfield(twix_obj.hdr.Dicom, 'tNucleus')
            nucleus = twix_obj.hdr.Dicom.tNucleus;
        end
    catch
    end
    if isempty(nucleus), nucleus = '1H'; end

    gammaTable = containers.Map( ...
        {'1H','2H','3He','7Li','13C','15N','17O','19F','23Na','31P','129Xe'}, ...
        {42.5760,6.5357,-32.4340,16.5470,10.7054,-4.3142,-5.7716,40.0541,11.2620,17.2350,-11.7770});

    if isKey(gammaTable, nucleus)
        gamma = gammaTable(nucleus);
    else
        gamma = 42.5760;
    end
end


% -------------------------------------------------------------------------
%   DIMENSION MAPPING
% -------------------------------------------------------------------------
function dims = buildDims(sqzDims, seqType)
    keys = {'t','f','coils','averages','timeinterleave', ...
            'kx','ky','kz','x','y','z','kpts','kshot','subSpecs','extras'};
    for k = 1:numel(keys), dims.(keys{k}) = 0; end

    switch seqType
        case 'concentric'
            % CONCENTRIC RINGS — inherits the io_MRSI_load_twix mapping.
            map = containers.Map( ...
                {'Col','Cha','Lin','Par','Rep','Ave','Seg','Eco'}, ...
                {'t', 'coils','kshot','extras','averages','averages','kshot','extras'});
            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isKey(map, lbl)
                    fld = map(lbl);
                    if dims.(fld) == 0
                        dims.(fld) = i;
                    elseif dims.extras == 0
                        dims.extras = i;
                    elseif dims.timeinterleave == 0
                        dims.timeinterleave = i;
                    end
                else
                    if dims.extras == 0
                        dims.extras = i;
                    elseif dims.timeinterleave == 0
                        dims.timeinterleave = i;
                    end
                end
            end

        case 'rosette'
            % ROSETTE — inherits the io_CSIload_twix3 mapping.
            % First unknown label becomes kshot; second becomes extras.
            fixed = struct('Col','t','Cha','coils','Ave','averages');
            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isfield(fixed, lbl)
                    dims.(fixed.(lbl)) = i;
                elseif strcmpi(lbl,'Set')
                    dims.kshot = i;            % Siemens often uses Set for shot/interleave
                elseif strcmpi(lbl,'Eco')
                    dims.extras = i;
                else
                    if dims.kshot == 0
                        dims.kshot = i;
                    elseif dims.extras == 0
                        dims.extras = i;
                    else
                        dims.timeinterleave = i;
                    end
                end
            end

        case 'cartesian'
            mp = containers.Map( ...
                {'Col','Cha','Ave','Rep','Seg','Phs','Lin','Sli'}, ...
                {'t'  ,'coils','averages','timeinterleave', ...
                 'kx' ,'ky'   ,'ky'      ,'kz'});
            for i = 1:numel(sqzDims)
                lbl = sqzDims{i};
                if isKey(mp, lbl)
                    dims.(mp(lbl)) = i;
                elseif dims.extras == 0
                    dims.extras = i;
                else
                    dims.timeinterleave = i;
                end
            end
    end
end

function [dims, data] = applyCanonicalPermutation(dims, data, seqType)
    if strcmp(seqType, 'cartesian')
        order = {'t','coils','averages','ky','kx','kz','timeinterleave','extras'};
    else
        % Both rosette and concentric land in the same canonical layout:
        order = {'t','coils','averages','kshot','extras','timeinterleave'};
    end

    perm = [];
    for k = 1:numel(order)
        idx = dims.(order{k});
        if idx ~= 0, perm(end+1) = idx; end %#ok<AGROW>
    end
    perm = unique(perm, 'stable');
    nd   = ndims(data);
    if numel(perm) < nd
        perm = [perm, setdiff(1:nd, perm, 'stable')];
    end

    if ~isequal(perm, 1:numel(perm))
        data = permute(data, perm);
    end

    % Re-index dims to their new positions.
    old   = dims;
    names = fieldnames(dims);
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


% -------------------------------------------------------------------------
%   GEOMETRY / METADATA
% -------------------------------------------------------------------------
function S = findAndSetFov(S, twix_obj)
    S.fov = struct();
    try
        asSlice = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1};
        S.fov.x = asSlice.dReadoutFOV;
        S.fov.y = asSlice.dPhaseFOV;
        S.fov.z = asSlice.dThickness;
    catch
        S.fov.x = 200;  S.fov.y = 200;  S.fov.z = 10;   % defensive defaults
    end
end

function scanDate = findScanDate(twix_obj)
    try
        m = regexp(twix_obj.hdr.MeasYaps.tReferenceImage0, ...
            '\.(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})', 'names');
        scanDate = datetime(str2double(m.year), str2double(m.month), str2double(m.day));
    catch
        scanDate = datetime('today');
    end
end

function S = calculateVoxelSize(S, numX, numY, numZ)
    S = setVoxelSize(S, 'x', S.fov.x / max(numX,1));
    S = setVoxelSize(S, 'y', S.fov.y / max(numY,1));
    S = setVoxelSize(S, 'z', S.fov.z / max(numZ,1));
end

function S = calculateVoxelCoordinates(S)
    xCoord = createCoords(S.fov.x/2, S.voxelSize.x) - S.imageOrigin(1);
    yCoord = createCoords(S.fov.y/2, S.voxelSize.y) - S.imageOrigin(2);
    S = setCoordinates(S, 'x', xCoord);
    S = setCoordinates(S, 'y', yCoord);
end

function coords = createCoords(halfFov, voxSize)
    if voxSize <= 0, coords = 0; return; end
    N      = round((2*halfFov)/voxSize);
    start  = -halfFov + voxSize/2;
    coords = start + (0:N-1)*voxSize;
end

function S = setVoxelSize(S, axisName, val)
    if ~isfield(S,'voxelSize'), S.voxelSize = struct(); end
    S.voxelSize.(axisName) = val;
end

function S = setCoordinates(S, axisName, vec)
    if ~isfield(S,'coordinates'), S.coordinates = struct(); end
    S.coordinates.(axisName) = vec(:).';
end

function S = findImageOrigin(S, twix_obj)
    % Populate imageOrigin from the prescribed slice CENTRE in MeasYaps.
    %
    % Why this matters: op_NUFFTSpatial1 computes its n_shift from
    % xCoords/yCoords, and xCoords = createCoords(...) - imageOrigin(1:2).
    % If imageOrigin stays at (0,0,0) (as it did with the original code on
    % XA60/VE+ where hdr.Config.VoI_Position_* is empty), the NUFFT never
    % compensates the off-centre slab prescription and the recon sits at
    % scanner isocentre instead of the slab centre.
    %
    % Sign convention: we negate the Cor component.  Empirically, the
    % xCoords/yCoords pipeline feeding into n_shift treats the dSag value
    % directly (positive pSag -> correct X shift) but expects dCor with the
    % opposite sign (negative pCor still wants a positive shift in FID-A
    % yC).  This matches the DICOM -> RAS flip that the affineMatrix in
    % calculateAffineMatrix already applies on the second row.
    %
    % The result has been verified against the spec2nii sform: with the
    % flip the recon centroid lands within ~0.6 voxels of the spec2nii
    % slab centre.
    S.imageOrigin = zeros(1,3);
    yaps = twix_obj.hdr.MeasYaps;
    paths_primary = {{'sSliceArray','asSlice',1,'sPosition','dSag'}, ...
                     {'sSliceArray','asSlice',1,'sPosition','dCor'}, ...
                     {'sSliceArray','asSlice',1,'sPosition','dTra'}};
    paths_fall    = {{'sSpecPara','sVoI','sPosition','dSag'}, ...
                     {'sSpecPara','sVoI','sPosition','dCor'}, ...
                     {'sSpecPara','sVoI','sPosition','dTra'}};
    for i = 1:3
        v = walk_path(yaps, paths_primary{i});
        if isempty(v), v = walk_path(yaps, paths_fall{i}); end
        if ~isempty(v), S.imageOrigin(i) = double(v); end
    end
    % large-bore table offsets (rare, present on whole-body scanners)
    tbl = {'lScanRegionPosSag','lScanRegionPosCor','lScanRegionPosTra'};
    for i = 1:3
        v = walk_path(yaps, {tbl{i}});
        if ~isempty(v), S.imageOrigin(i) = S.imageOrigin(i) + double(v); end
    end
    % DICOM Cor -> FID-A yC sign flip (see header note)
    S.imageOrigin(2) = -S.imageOrigin(2);
end

function S = calculateAffineMatrix(S, twix_obj)
    try
        [zVec, theta] = getZVectorAndTheta(twix_obj);
        R = getRotationMatrixFromVector(zVec(1), zVec(2), zVec(3), theta);
        R(4,4) = 1;

        yVec = [0, 1, -zVec(2)/max(zVec(3), eps)];
        yVec = yVec / norm(yVec);
        yVec = R * [yVec'; 1];  yVec = yVec(1:3);
        xVec = cross(yVec, zVec);

        Arot = [xVec', yVec, zVec', [0;0;0]; 0 0 0 1];
        Ascl = diag([S.voxelSize.x, S.voxelSize.y, S.voxelSize.z, 1]);
        Atrs = eye(4);
        Atrs(1,4) = -S.imageOrigin(1) - S.fov.x/2;
        Atrs(2,4) = -S.imageOrigin(2) - S.fov.y/2;
        Atrs(3,4) =  S.imageOrigin(3) - S.fov.z/2;

        S.affineMatrix = Arot * Atrs * Ascl;
    catch
        S.affineMatrix = diag([S.voxelSize.x, S.voxelSize.y, S.voxelSize.z, 1]);
    end
end

function [z_vect, theta] = getZVectorAndTheta(twix_obj)
    z_vect = [0 0 1]; theta = 0;
    try
        if isfield(twix_obj.hdr.Meas, 'VoI_Normal_Sag')
            z_vect(1) = -nz(twix_obj.hdr.Meas.VoI_Normal_Sag);
            z_vect(2) = -nz(twix_obj.hdr.Meas.VoI_Normal_Cor);
            z_vect(3) =  no(twix_obj.hdr.Meas.VoI_Normal_Tra);
            theta     =  nz(twix_obj.hdr.Meas.VoI_InPlaneRotAngle);
        elseif isfield(twix_obj.hdr.Meas, 'VoiNormalSag')
            z_vect(1) = -nz(twix_obj.hdr.Meas.VoiNormalSag);
            z_vect(2) = -nz(twix_obj.hdr.Meas.VoiNormalCor);
            z_vect(3) =  no(twix_obj.hdr.Meas.VoiNormalTra);
            theta     =  nz(twix_obj.hdr.Meas.VoiInPlaneRot);
        end
    catch
    end
end

function R = getRotationMatrixFromVector(x, y, z, theta)
    v = [x,y,z]; v = v / max(norm(v), eps);
    c = cos(theta); s = sin(theta); C = 1 - c;
    R = [ ...
        c+v(1)^2*C,           v(1)*v(2)*C - v(3)*s, v(1)*v(3)*C + v(2)*s; ...
        v(2)*v(1)*C + v(3)*s, c+v(2)^2*C,           v(2)*v(3)*C - v(1)*s; ...
        v(3)*v(1)*C - v(2)*s, v(3)*v(2)*C + v(1)*s, c+v(3)^2*C ];
end

function v = nz(x), if isempty(x), v = 0; else, v = x; end, end
function v = no(x), if isempty(x), v = 1; else, v = x; end, end


function S = setDefaultFlagValues(S, isCartesian)
    S.flags.writtentostruct = 1;
    S.flags.gotparams       = 1;
    S.flags.leftshifted     = 0;
    S.flags.filtered        = 0;
    S.flags.zeropadded      = 0;
    S.flags.freqcorrected   = 0;
    S.flags.phasecorrected  = 0;
    S.flags.averaged        = 0;
    S.flags.addedrcvrs      = 0;
    S.flags.subtracted      = 0;
    S.flags.writtentotext   = 0;
    S.flags.downsampled     = 0;
    S.flags.spatialFT       = 0;
    S.flags.spectralFT      = 0;
    S.flags.coilCombined    = 0;
    S.flags.isFourSteps     = 0;
    S.flags.isCartesian     = isCartesian;
end


% -------------------------------------------------------------------------
%   FOV (VOI) OFFSET — used by op_CSIShift (rosette path only)
% -------------------------------------------------------------------------
function [xShift_mm, yShift_mm] = computeFOVShift(twix_obj)
    if iscell(twix_obj), hdr = twix_obj{2}.hdr; else, hdr = twix_obj.hdr; end

    nSag = 0; nCor = 0; nTra = 0;
    if isfield(hdr.MeasYaps,'sSpecPara') && isfield(hdr.MeasYaps.sSpecPara,'sVoI') ...
       && isfield(hdr.MeasYaps.sSpecPara.sVoI,'sNormal')
        ns = hdr.MeasYaps.sSpecPara.sVoI.sNormal;
        if isfield(ns,'dSag'), nSag = ns.dSag; end
        if isfield(ns,'dCor'), nCor = ns.dCor; end
        if isfield(ns,'dTra'), nTra = ns.dTra; end
    end
    normal     = [nSag, nCor, nTra];
    [~, idx]   = max(abs(normal));
    planes     = {'Sagittal', 'Coronal', 'Axial'};
    plane      = planes{idx};
    switch plane
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

    switch plane
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


% Walk a MeasYaps path like {'sSliceArray','asSlice',1,'sPosition','dSag'}.
% Returns the value (numeric scalar) or [] if any step is missing.
function v = walk_path(s, path)
    v = []; cur = s;
    for k = 1:numel(path)
        p = path{k};
        if isnumeric(p)
            if iscell(cur) && numel(cur) >= p,    cur = cur{p};
            elseif isstruct(cur) && numel(cur) >= p, cur = cur(p);
            else, return; end
        else
            if isstruct(cur) && isfield(cur, p), cur = cur.(p);
            else, return; end
        end
        if isempty(cur), return; end
    end
    if isnumeric(cur) && isscalar(cur), v = double(cur); end
end
