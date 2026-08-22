function MRSIStruct = op_CSItimeZeroFill(MRSIStruct, Ntarget)
%OP_CSITIMEZEROFILL  Zero-fill an MRSI FID along the TIME dimension to Ntarget pts.
%
%   MRSIStruct = op_CSItimeZeroFill(MRSIStruct, 4096)
%
%   Pads the FID with zeros (from its current length up to Ntarget) along dims.t
%   for every voxel. Interpolates the spectrum (finer ppm sampling, smoother
%   lineshapes) without adding information. Call on the TIME-domain struct
%   (ccav / ccav_w) BEFORE op_CSIFourierTransform. Dwell time is unchanged, so
%   spectralWidth is unchanged; only the number of points grows.
%
%   (Distinct from op_CSIZeroFill, which zero-fills the SPATIAL kx/ky dims.)

    td  = MRSIStruct.dims.t;
    assert(td > 0, 'op_CSItimeZeroFill needs a time dimension (dims.t).');
    Nt0 = size(MRSIStruct.data, td);
    if Ntarget <= Nt0
        fprintf('op_CSItimeZeroFill: target %d <= current %d, nothing to do.\n', Ntarget, Nt0);
        return
    end

    sz      = size(MRSIStruct.data);
    sz(td)  = Ntarget;
    newdata = zeros(sz, 'like', MRSIStruct.data);
    idx     = repmat({':'}, 1, ndims(MRSIStruct.data));
    idx{td} = 1:Nt0;
    newdata(idx{:}) = MRSIStruct.data;

    MRSIStruct.data = newdata;
    MRSIStruct.sz   = size(newdata);

    Nt = Ntarget;
    if isfield(MRSIStruct,'spectralDwellTime') && ~isempty(MRSIStruct.spectralDwellTime)
        MRSIStruct.spectralTime = (0:Nt-1) * MRSIStruct.spectralDwellTime;
    end
    if isfield(MRSIStruct,'adcDwellTime') && ~isempty(MRSIStruct.adcDwellTime)
        MRSIStruct.adcTime = (0:Nt-1) * MRSIStruct.adcDwellTime;
    end
    if isfield(MRSIStruct,'flags'), MRSIStruct.flags.zerofilled = 1; end

    fprintf('op_CSItimeZeroFill: %d -> %d time points (dwell unchanged).\n', Nt0, Nt);
end
