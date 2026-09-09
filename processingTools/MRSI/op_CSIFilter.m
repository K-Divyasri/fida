% op_CSIFilter.m
%
% USAGE:
% [out,lor] = op_CSIFilter(in,lb);
%
% DESCRIPTION:
% Time-domain exponential (Lorentzian) line broadening for MRSI data — the CSI
% analogue of op_filter. Multiplies the FID along dims.t by exp(-t*pi*lb), i.e.
% a Lorentzian with T2 = 1/(pi*lb). Broadens every voxel's lineshape by lb Hz.
% Call on the TIME-domain struct (ccav / ccav_w) BEFORE op_CSIFourierTransform.
%
% INPUTS:
% in     = input MRSI FID-A structure (time domain).
% lb     = line broadening factor in Hz.
%
% OUTPUTS:
% out    = MRSI structure with the filter applied to in.data.
% lor    = the exponential (time-domain) filter envelope [Nt x 1] that was used.

function [out, lor] = op_CSIFilter(in, lb)

    if nargin < 2 || isempty(lb), lb = 0; end
    out = in;  lor = [];
    if lb == 0, return; end

    % --- must be time domain (the filter multiplies the FID) ---
    if isfield(in,'flags') && isfield(in.flags,'spectralFT') && in.flags.spectralFT
        error(['op_CSIFilter: data is in the spectral domain. Apply the filter ' ...
               'to the FID BEFORE op_CSIFourierTransform.']);
    end
    if isfield(in,'flags') && isfield(in.flags,'filtered') && in.flags.filtered
        warning('op_CSIFilter: data already filtered; applying line broadening again.');
    end

    td = in.dims.t;
    assert(td > 0, 'op_CSIFilter needs a time dimension (dims.t).');
    Nt = size(in.data, td);

    % --- spectral dwell time (s) ---
    if isfield(in,'spectralWidth') && ~isempty(in.spectralWidth) && in.spectralWidth > 0
        dt = 1 / in.spectralWidth;
    elseif isfield(in,'spectralDwellTime') && ~isempty(in.spectralDwellTime) && in.spectralDwellTime > 0
        dt = in.spectralDwellTime;
    else
        error('op_CSIFilter: need spectralWidth or spectralDwellTime to build the time axis.');
    end

    % --- Lorentzian (exponential) envelope, T2 = 1/(pi*lb) ---
    t   = (0:Nt-1).' * dt;              % [Nt x 1] s
    lor = exp(-pi * lb * t);           % same convention as op_filter

    % --- broadcast along the time dimension and apply ---
    shape       = ones(1, ndims(in.data));
    shape(td)   = Nt;
    out.data    = in.data .* reshape(lor, shape);

    if isfield(out,'flags'), out.flags.filtered = 1; end
    fprintf('op_CSIFilter: applied %g Hz exponential line broadening (%d time pts).\n', lb, Nt);
end
