function MRSIStruct = op_CSIleftshift(MRSIStruct, ls)
%OP_CSILEFTSHIFT  Drop leading FID points from an MRSI struct (1st-order phase).
%
%   MRSIStruct = op_CSIleftshift(MRSIStruct)       % ls = MRSIStruct.pointsToLeftshift
%   MRSIStruct = op_CSIleftshift(MRSIStruct, ls)   % drop ls points explicitly
%
%   The FID's first point isn't at t=0 (ADC/echo delay), which puts a large
%   first-order phase across the spectrum and distorts lineshapes. Removing the
%   leading ls points before the spectral FT flattens that phase. MRSI version
%   of op_leftshift (which is single-voxel and uses .fids); this works on
%   .data along dims.t for all voxels at once.
%
%   Call it on the TIME-domain struct (ccav / ccav_w), BEFORE op_CSIFourierTransform.

    if nargin < 2 || isempty(ls)
        ls = 0;
        if isfield(MRSIStruct,'pointsToLeftshift') && ~isempty(MRSIStruct.pointsToLeftshift)
            ls = double(MRSIStruct.pointsToLeftshift);
        end
    end
    ls = round(ls);
    if ls <= 0
        fprintf('op_CSIleftshift: ls=%d, nothing to do.\n', ls); return
    end
    if isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'leftshifted') ...
            && MRSIStruct.flags.leftshifted
        fprintf('op_CSIleftshift: already left-shifted, skipping.\n'); return
    end

    td = MRSIStruct.dims.t;
    assert(td > 0, 'op_CSIleftshift needs a time dimension (dims.t).');
    Nt0 = size(MRSIStruct.data, td);
    assert(ls < Nt0, 'ls (%d) must be < number of time points (%d).', ls, Nt0);

    % drop the first ls samples along the time dimension, keep all voxels
    idx      = repmat({':'}, 1, ndims(MRSIStruct.data));
    idx{td}  = (ls+1):Nt0;
    MRSIStruct.data = MRSIStruct.data(idx{:});
    MRSIStruct.sz   = size(MRSIStruct.data);

    Nt = MRSIStruct.sz(td);
    if isfield(MRSIStruct,'spectralDwellTime') && ~isempty(MRSIStruct.spectralDwellTime)
        MRSIStruct.spectralTime = (0:Nt-1) * MRSIStruct.spectralDwellTime;
    end
    if isfield(MRSIStruct,'adcDwellTime') && ~isempty(MRSIStruct.adcDwellTime)
        MRSIStruct.adcTime = (0:Nt-1) * MRSIStruct.adcDwellTime;
    end
    MRSIStruct.flags.leftshifted = 1;

    fprintf('op_CSIleftshift: dropped %d leading FID point(s); t dim %d -> %d.\n', ...
            ls, Nt0, Nt);
end
