function MRS = op_CSItoMRS(MRSIStruct, xCoordinate, yCoordinate, index)
    % Accept a positional struct for 'index' (or nothing)
    if nargin < 4 || isempty(index), index = struct; end
    index = normalizeIndexStruct(index);

    checkArguments(MRSIStruct, index);
    MRS = MRSIStruct;

    % Build voxel index (x,y plus optional average/coil/subspec/extra)
    MRSIndex = buildIndex(MRSIStruct, xCoordinate, yCoordinate, index);

    % --- Pick spectral axis: prefer 't', otherwise fall back to 'f'
    dim_t = getDimension(MRSIStruct, 't');
    dim_f = getDimension(MRSIStruct, 'f');
    tfDim = pickTFDim(dim_t, dim_f);

    % Slice the voxel
    data = getData(MRSIStruct);
    if getFlags(MRSIStruct, 'spectralft') == true
        MRS.specs = data(MRSIndex{:});
        % Keep your original direction; swap to ifft/ifftshift if you prefer.
        MRS.fids  = fft(fftshift(MRS.specs, tfDim), [], tfDim);
    else
        MRS.fids  = data(MRSIndex{:});
        MRS.specs = fftshift(ifft(MRS.fids, [], tfDim), tfDim);
    end

    % Trim and squeeze
    if isfield(MRS,'data'), MRS = rmfield(MRS,'data'); end
    MRS.specs = squeeze(MRS.specs);

    % --- Axes/metadata
    Ntf = size(MRS.specs, tfDim);

    % Time axis
    if isfield(MRSIStruct,'spectralTime') && ~isempty(MRSIStruct.spectralTime)
        MRS.t = MRSIStruct.spectralTime;
    else
        if isfield(MRSIStruct,'spectralDwellTime') && ~isempty(MRSIStruct.spectralDwellTime)
            dt = MRSIStruct.spectralDwellTime;
        elseif isfield(MRSIStruct,'spectralWidth') && ~isempty(MRSIStruct.spectralWidth) && MRSIStruct.spectralWidth>0
            dt = 1/MRSIStruct.spectralWidth;
        else
            error('Cannot infer spectral time axis: missing spectralDwellTime and spectralWidth.');
        end
        MRS.t = (0:Ntf-1)*dt;
    end

    % Dwell time / spectral width
    if isfield(MRSIStruct,'spectralDwellTime') && ~isempty(MRSIStruct.spectralDwellTime)
        MRS.dwelltime = MRSIStruct.spectralDwellTime;
    else
        if isfield(MRSIStruct,'spectralWidth') && ~isempty(MRSIStruct.spectralWidth) && MRSIStruct.spectralWidth>0
            MRS.dwelltime = 1/MRSIStruct.spectralWidth;
        else
            error('Missing spectralDwellTime and spectralWidth: cannot set dwelltime.');
        end
    end

    if isfield(MRSIStruct,'spectralWidth') && ~isempty(MRSIStruct.spectralWidth)
        MRS.spectralwidth = MRSIStruct.spectralWidth;
    else
        MRS.spectralwidth = 1/MRS.dwelltime;
    end

    % ppm axis
    if isfield(MRSIStruct,'ppm')
        MRS.ppm = MRSIStruct.ppm;
    else
        stepHz    = MRS.spectralwidth / Ntf;
        freqBound = MRS.spectralwidth/2 - stepHz/2;
        freqHz    = -freqBound:stepHz:freqBound;
        if ~(isfield(MRSIStruct,'Bo') && isfield(MRSIStruct,'gamma') ...
             && ~isempty(MRSIStruct.Bo) && ~isempty(MRSIStruct.gamma))
            error('Missing Bo or gamma to construct ppm axis.');
        end
        MRS.ppm = freqHz / (MRSIStruct.Bo * MRSIStruct.gamma);
    end

    MRS.sz = size(MRS.specs);

    % Remove spatial dims and any explicitly indexed dims
    MRS = removeDimension(MRS, 'x');
    MRS = removeDimension(MRS, 'y');
    if index.averageIndex >= 1, MRS = removeDimension(MRS, 'averages'); end
    if index.coilIndex    >= 1, MRS = removeDimension(MRS, 'coils');    end
    if index.subSpecIndex >= 1, MRS = removeDimension(MRS, 'subspec');  end
    if index.extraIndex   >= 1, MRS = removeDimension(MRS, 'extras');   end

    % Relabel the spectral axis of the single-voxel output as 't' and drop the
    % MRSI 'f' label.  After squeezing one voxel, the spectral data lies along
    % dimension 1, so this matches the standard FID-A SVS convention (dims.t)
    % and lets SVS tools (op_plotspec, io_writelcm, ...) consume the voxel
    % directly.  The data itself is unchanged; only the dimension label moves
    % from 'f' to 't'.
    MRS = setDimension(MRS, 't', 1);
    MRS = setDimension(MRS, 'f', 0);
end

function s = normalizeIndexStruct(s)
    def = struct('averageIndex',0,'coilIndex',0,'subSpecIndex',0,'extraIndex',0,'linearIndex',0);
    fn  = fieldnames(def);
    for i=1:numel(fn)
        f = fn{i};
        if ~isfield(s,f) || isempty(s.(f)), s.(f) = def.(f); end
        s.(f) = double(s.(f));
    end
end

function tfDim = pickTFDim(dim_t, dim_f)
    if dim_t > 0
        tfDim = dim_t;
    elseif dim_f > 0
        tfDim = dim_f;
    else
        error('Neither ''t'' nor ''f'' dimension found in the input structure.');
    end
end

function checkArguments(in, index)
    if (getFlags(in, 'spatialft') == false)
        error('Please Fourier transform along the spatial and spectral dimensions before using this function.');
    end
    if (index.linearIndex > 0 && (index.averageIndex || index.coilIndex || index.subSpecIndex || index.extraIndex))
        % no-op retained for compatibility
    end
end

function MRSIndex = buildIndex(MRSIStruct, xCoordinate, yCoordinate, index)
    MRSIndex = repmat({':'}, 1, numel(MRSIStruct.sz));
    MRSIndex{MRSIStruct.dims.x} = xCoordinate;
    MRSIndex{MRSIStruct.dims.y} = yCoordinate;

    dA = getDimension(MRSIStruct, 'averages');
    if dA>0 && index.averageIndex>0, MRSIndex{dA} = index.averageIndex; end

    dC = getDimension(MRSIStruct, 'coils');
    if dC>0 && index.coilIndex>0,    MRSIndex{dC} = index.coilIndex;    end

    dS = getDimension(MRSIStruct, 'subspec');     % singular first
    if dS==0, dS = getDimension(MRSIStruct, 'subspecs'); end
    if dS>0 && index.subSpecIndex>0, MRSIndex{dS} = index.subSpecIndex; end

    dE = getDimension(MRSIStruct, 'extras');
    if dE>0 && index.extraIndex>0,   MRSIndex{dE} = index.extraIndex;   end
end
