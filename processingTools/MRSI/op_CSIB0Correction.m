function [Sout, phaseMap, freqMap] = op_CSIB0Correction(Sin, varargin)
% op_CSIB0Correction
% Compute and/or apply B0 (off-resonance) and static phase corrections for CSI.
%
% USAGE
%   [Sout, phaseMap, freqMap] = op_CSIB0Correction(Swater)
%     -> builds maps from Swater (preferably water-unsuppressed) and applies them to Swater
%
%   Sout = op_CSIB0Correction(Sin, phaseMap, freqMap)
%     -> applies precomputed maps to Sin
%
% INPUT
%   Sin:  FID-A style CSI struct with fields:
%         - data (complex) in either time or frequency domain
%         - sz, dims (with dims.x and dims.y nonzero)
%         - time info in adcTime/spectralTime or dwell times / spectralWidth
%
% OPTIONAL NAME-VALUE (when computing maps)
%   'Mask'          : logical [Ny,Nx], restrict fit/application to mask (default: all true)
%   'MaxTimePoints' : integer, number of earliest time points for the linear fit (default: auto)
%   'RefChoice'     : 'center' or 'maxmag' (only used if phase map is also desired by alignment)
%                     NOTE: For this implementation, phase map comes from intercept term, so
%                     we do not need a separate reference voxel. 'RefChoice' kept for parity.
%
% OUTPUT
%   Sout    : B0-corrected struct (same domain as input; flags updated)
%   phaseMap: [Ny,Nx] static phase (radians), estimated at t=0
%   freqMap : [Ny,Nx] off-resonance (Hz), slope/2π from unwrapped FID phase vs time
%
% NOTES
% - If Sin is in frequency domain, the routine internally IFFT->time, corrects, then FFT back.
% - Correction applied is: fid_corr(t) = fid(t) * exp(-1j*(phase + 2*pi*freq*t)).
% - If maps are provided, their shape must be [Ny,Nx] matching dims.y/x of Sin.

% -------- Parse inputs
p = inputParser;
addRequired(p,'Sin',@(s)isstruct(s) && isfield(s,'data') && isfield(s,'sz') && isfield(s,'dims'));
addOptional(p,'phaseMap',[],@(a)isnumeric(a)||islogical(a));
addOptional(p,'freqMap',[],@(a)isnumeric(a));
addParameter(p,'Mask',[],@(m)islogical(m)||isempty(m));
addParameter(p,'MaxTimePoints',[],@(n)isempty(n)||(isscalar(n)&&n>10));
addParameter(p,'RefChoice','center',@(s)ischar(s)||isstring(s));
parse(p, Sin, varargin{:});
Sin       = p.Results.Sin;
phaseMapI = p.Results.phaseMap;
freqMapI  = p.Results.freqMap;
mask      = p.Results.Mask;
maxTPts   = p.Results.MaxTimePoints; %#ok<NASGU> (kept for API clarity)
% RefChoice kept for compatibility (not used here explicitly)

% -------- Validate dims
dims = Sin.dims;
if ~isfield(dims,'x') || ~isfield(dims,'y') || dims.x==0 || dims.y==0
    error('op_CSIB0Correction: dims.x and dims.y must be defined and non-zero.');
end
sz  = Sin.sz;
Nx  = sz(dims.x);
Ny  = sz(dims.y);

% -------- Build/verify mask
if isempty(mask)
    mask = true(Ny,Nx);
else
    if ~isequal(size(mask),[Ny,Nx])
        error('Mask must be size [Ny,Nx] = [%d,%d].',Ny,Nx);
    end
end

% -------- Determine spectral/time dimension and domain
tDim = getdim(dims,'t');
fDim = getdim(dims,'f');
hasT = tDim>0;
hasF = fDim>0;

if ~(hasT || hasF)
    % assume first dimension is spectral if nothing given
    fDim = 1; hasF = true;
end

% -------- Construct time vector (seconds)
[tvec, dt] = get_time_vector(Sin, tDim, fDim, sz); %#ok<ASGLU>
Nt = numel(tvec);

% -------- Prepare permutation so time is dimension 1 and [y,x] are 2,3
permToTyx  = make_perm_order(numel(sz), tDim, dims.y, dims.x);
invPerm    = invert_perm(permToTyx);

% -------- Bring data to time domain, shape [Nt, Ny, Nx, ...rest...]
data = Sin.data;
% Move dims -> [t, y, x, rest...]
if hasT
    data_tyxl = permute(data, permToTyx);
    data_t = data_tyxl;
else
    % Frequency -> time along fDim after permute
    data_fyxl = permute(data, make_perm_order(numel(sz), fDim, dims.y, dims.x));
    data_t = ifftshift(data_fyxl, 1);
    data_t = ifft(data_t, [], 1); % spectrum -> fid
end

% Ensure complex
if ~isfloat(data_t)
    data_t = double(data_t);
end
if ~isreal(data_t)
    % ok
else
    % some pipelines store real spectrum; we still allow but warn
    warning('Input data appears real. B0 estimation expects complex data.');
    data_t = complex(data_t, 0);
end

% Record trailing shape (...rest...)
sz_t = size(data_t);
if numel(sz_t)<3, sz_t(3)=1; end
NyN  = sz_t(2);
NxN  = sz_t(3);
if NyN~=Ny || NxN~=Nx
    error('Unexpected spatial shape after permutation. Got [%d,%d], expected [%d,%d].',NyN,NxN,Ny,Nx);
end
trailingSz = sz_t(4:end);
nVox  = Ny*Nx;
nTail = prod(trailingSz);
data_t = reshape(data_t, [Nt, Ny, Nx, nTail]);

% -------- If maps provided: apply and return
if ~isempty(phaseMapI) && ~isempty(freqMapI)
    [phaseMap, freqMap] = validate_maps(phaseMapI, freqMapI, Ny, Nx);
    Sout = Sin;
    % broadcast maps across trailing dims
    phaseMapB = repmat(reshape(phaseMap,[1,Ny,Nx,1]), [Nt,1,1,nTail]);
    freqMapB  = repmat(reshape(freqMap, [1,Ny,Nx,1]), [Nt,1,1,nTail]);
    corr = exp(-1j*(phaseMapB + 2*pi*freqMapB.*tvec(:)));
    data_t_corr = data_t .* corr;

    % reshape back
    data_t_corr = reshape(data_t_corr, [Nt, Ny, Nx, trailingSz]);

    % Return to input domain & original dim order
    if hasT
        data_out = ipermute(data_t_corr, permToTyx);
    else
        spec = fft(data_t_corr, [], 1);
        spec = fftshift(spec, 1);
        data_out = ipermute(spec, make_perm_order(numel(sz), fDim, dims.y, dims.x));
    end
    Sout.data = data_out;
    Sout = set_flags_corrected(Sout);
    return;
end

% -------- Otherwise: compute maps from current Sin (assumed water ref)
% We estimate per-voxel linear phase: unwrap(angle(fid(t))) ~ phi0 + 2π*freq*t
phaseMap = zeros(Ny,Nx);
freqMap  = zeros(Ny,Nx);

% Choose a reasonable number of early points for fitting if not provided
% Heuristic: use up to min( max(64, ceil(Nt*0.25)), Nt )
max_tpt = min(max(64, ceil(0.25*Nt)), Nt);
tfit = tvec(1:max_tpt);

% Flatten trailing dims while computing maps: average across coils/subspec/etc if present
if nTail>1
    % coil/subspec averaging on the FID for more robust fit
    data_t_mean = mean(data_t, 4, 'omitnan');
else
    data_t_mean = data_t;
end

% Compute maps only inside mask
for y = 1:Ny
    for x = 1:Nx
        if ~mask(y,x), continue; end
        fid = data_t_mean(:,y,x);             % Nt x 1
        if all(abs(fid)<eps)
            % empty voxel
            phaseMap(y,x) = 0;
            freqMap(y,x)  = 0;
            continue;
        end

        ph = unwrap(angle(fid(1:max_tpt)));
        % Robust linear fit: ph ~ a*t + b
        % Use simple least squares; could switch to robustfit if available
        tt = tfit(:);
        A  = [tt, ones(numel(tt),1)];
        % Solve in double
        ab = A \ double(ph(:));
        a  = ab(1);  % slope (rad/s)
        b  = ab(2);  % intercept at t=0 (rad)

        freqMap(y,x)  = a/(2*pi);  % Hz
        phaseMap(y,x) = b;         % rad
    end
end

% -------- Apply maps to Sin and return
Sout = Sin;
phaseMapB = repmat(reshape(phaseMap,[1,Ny,Nx,1]), [Nt,1,1,nTail]);
freqMapB  = repmat(reshape(freqMap, [1,Ny,Nx,1]), [Nt,1,1,nTail]);
corr = exp(-1j*(phaseMapB + 2*pi*freqMapB.*tvec(:)));

data_t_corr = data_t .* corr;
data_t_corr = reshape(data_t_corr, [Nt, Ny, Nx, trailingSz]);

if hasT
    data_out = ipermute(data_t_corr, permToTyx);
else
    spec = fft(data_t_corr, [], 1);
    spec = fftshift(spec, 1);
    data_out = ipermute(spec, make_perm_order(numel(sz), fDim, dims.y, dims.x));
end

Sout.data = data_out;
Sout = set_flags_corrected(Sout);

end % ===== main =====


% ---------- helpers ----------
function d = getdim(dims, name)
if isfield(dims,name)
    d = dims.(name);
else
    d = 0;
end
end

function [tvec, dt] = get_time_vector(S, tDim, fDim, sz)
% Returns tvec in seconds, length = Nt (the spectral length)
% Prefers explicit time vectors; falls back to dwell times or spectralWidth.
% Works whether current domain is time or frequency.
Nt = sz(max([tDim,fDim,1]));
% 1) explicit time vectors
if isfield(S,'spectralTime') && ~isempty(S.spectralTime) && numel(S.spectralTime)==Nt
    tvec = S.spectralTime(:); dt = tvec(2)-tvec(1); return;
end
if isfield(S,'adcTime') && ~isempty(S.adcTime)
    % adcTime might be length of raw acquisition; use first Nt points
    if numel(S.adcTime) >= Nt
        tvec = S.adcTime(1:Nt).'; dt = tvec(2)-tvec(1); return;
    end
end
% 2) dwell times
if isfield(S,'spectralDwellTime') && ~isempty(S.spectralDwellTime) && S.spectralDwellTime>0
    dt = S.spectralDwellTime;
    tvec = (0:Nt-1)' * dt; return;
end
if isfield(S,'adcDwellTime') && ~isempty(S.adcDwellTime) && S.adcDwellTime>0
    dt = S.adcDwellTime;
    tvec = (0:Nt-1)' * dt; return;
end
% 3) spectralWidth
if isfield(S,'spectralWidth') && ~isempty(S.spectralWidth) && S.spectralWidth>0
    dt = 1/S.spectralWidth;
    tvec = (0:Nt-1)' * dt; return;
end
% 4) fallback
dt = 1e-4;
tvec = (0:Nt-1)' * dt;
warning('No time info found; using default dwell time of 0.1 ms.');
end

function p = make_perm_order(nd, firstDim, yDim, xDim)
% Create permutation so dims become [firstDim, yDim, xDim, rest...]
all = 1:nd;
use = [firstDim, yDim, xDim];
use = use(use>0);
% Remove duplicates (if firstDim==yDim etc.)
use = unique(use,'stable');
rest = setdiff(all, use, 'stable');
p = [use, rest];
% Ensure firstDim is first; if it was 0, force 1 as first (already handled upstream)
if isempty(firstDim) || firstDim==0
    if p(1)~=1
        % ensure something occupies first spot (won't happen with upstream guards)
    end
end
end

function invp = invert_perm(p)
invp = zeros(size(p));
invp(p) = 1:numel(p);
end

function [phaseMap, freqMap] = validate_maps(P, F, Ny, Nx)
if ~isequal(size(P),[Ny,Nx]) || ~isequal(size(F),[Ny,Nx])
    error('Provided phaseMap/freqMap must both be size [Ny,Nx]=[%d,%d].',Ny,Nx);
end
phaseMap = double(P);
freqMap  = double(F);
end

function S = set_flags_corrected(S)
if ~isfield(S,'flags') || ~isstruct(S.flags)
    S.flags = struct();
end
S.flags.phasecorrected = 1;
S.flags.freqcorrected  = 1;
end
