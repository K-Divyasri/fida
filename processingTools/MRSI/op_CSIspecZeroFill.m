function MRSIStruct = op_CSIspecZeroFill(MRSIStruct, Ntarget)
%OP_CSISPECZEROFILL  Zero-fill a spectral (frequency-domain) MRSI struct to
%   Ntarget points by padding the underlying FID, then re-transforming.
%
%   MRSIStruct = op_CSIspecZeroFill(ftSpec, 4096)
%
%   Use this AFTER the spectral FT (e.g. right after B0 correction) instead of
%   zero-filling the raw FID before the FT. That way the expensive per-voxel
%   steps (lipid removal, B0 correction) run on the short NATIVE FID (fast,
%   and no zero-padded region to distort), and only the final spectrum is
%   interpolated to the finer grid for LCModel / the viewer.
%
%   It inverts op_CSIFourierTransform's spectral transform
%   (spec = fftshift(ifft(fid))), pads the FID with zeros along the spectral
%   dimension, then re-applies the forward transform. This is a pure
%   interpolation of the spectrum: it adds no information, matching a
%   time-domain zero-fill of the native data.
%
%   (Distinct from op_CSItimeZeroFill, which pads a TIME-domain struct, and
%    from op_CSIZeroFill, which pads the SPATIAL kx/ky dims.)

    fdim = MRSIStruct.dims.f;
    if fdim == 0, fdim = MRSIStruct.dims.t; end
    assert(fdim > 0, 'op_CSIspecZeroFill needs a spectral (f) or time (t) dimension.');

    N0 = size(MRSIStruct.data, fdim);
    if Ntarget <= N0
        fprintf('op_CSIspecZeroFill: target %d <= current %d, nothing to do.\n', Ntarget, N0);
        return
    end

    % --- frequency -> time (inverse of fftshift(ifft(fid)) used by the FT) ---
    fid = fft(ifftshift(MRSIStruct.data, fdim), [], fdim);

    % --- pad the FID with zeros to Ntarget along the spectral dimension ------
    sz        = size(fid);
    sz(fdim)  = Ntarget;
    fidPad    = zeros(sz, 'like', fid);
    idx       = repmat({':'}, 1, ndims(fid));
    idx{fdim} = 1:N0;
    fidPad(idx{:}) = fid;

    % --- time -> frequency (forward transform, same convention as FT) -------
    MRSIStruct.data = fftshift(ifft(fidPad, [], fdim), fdim);
    MRSIStruct.sz   = size(MRSIStruct.data);

    % --- update spectral / time axes for the new point count ----------------
    Nt = Ntarget;
    if isfield(MRSIStruct,'spectralDwellTime') && ~isempty(MRSIStruct.spectralDwellTime)
        MRSIStruct.spectralTime = (0:Nt-1) * MRSIStruct.spectralDwellTime;
    end
    if isfield(MRSIStruct,'adcDwellTime') && ~isempty(MRSIStruct.adcDwellTime)
        MRSIStruct.adcTime = (0:Nt-1) * MRSIStruct.adcDwellTime;
    end

    % ppm axis: same formula as op_CSIFourierTransform/calculatePPM
    if isfield(MRSIStruct,'spectralWidth') && ~isempty(MRSIStruct.spectralWidth) ...
            && isfield(MRSIStruct,'Bo') && isfield(MRSIStruct,'gamma') ...
            && ~isempty(MRSIStruct.Bo) && ~isempty(MRSIStruct.gamma)
        sw   = MRSIStruct.spectralWidth;
        step = sw / Nt;
        f    = (-sw/2 + step/2) : step : (sw/2 - step/2);
        ppm  = -f / (MRSIStruct.Bo * MRSIStruct.gamma);
        if isfield(MRSIStruct,'nucleus') && strcmp(MRSIStruct.nucleus,'1H')
            ppm = ppm + 4.65;
        end
        MRSIStruct.ppm = ppm;
    end

    if isfield(MRSIStruct,'flags'), MRSIStruct.flags.zerofilled = 1; end

    fprintf('op_CSIspecZeroFill: %d -> %d spectral points (interpolated).\n', N0, Nt);
end
