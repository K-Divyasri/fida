function MRSIStruct = op_CSIaddPhase(MRSIStruct, ph0, ph1, ppm0, suppressPlot)
% op_CSIaddPhase  Zero- and first-order phase correction for CSI/MRSI data.
%
% USAGE:
%   out = op_CSIaddPhase(in, ph0, ph1, ppm0, suppressPlot)
%
% INPUTS:
%   ph0          = zero-order phase (degrees)
%   ph1          = first-order phase (SECONDS) about ppm0
%   ppm0         = pivot ppm for the first-order phase (default 4.65 = water)
%   suppressPlot = suppress the before/after plot (default true)
%
% Robust: applies the phase along the ACTUAL spectral dimension (dims.f, else
% dims.t) — not a hard-coded dim — and converts ppm->Hz with txfrq (falls back
% to Bo*gamma). Works whether the struct is in the time or spectral domain.

    arguments
        MRSIStruct   (1,1) struct
        ph0          (1,1) double
        ph1          (1,1) double = 0
        ppm0         (1,1) double = 4.65
        suppressPlot (1,1) logical = true
    end

    dims = MRSIStruct.dims;
    wasFreqDomain = isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'spectralFT') ...
                    && MRSIStruct.flags.spectralFT;

    % --- spectral dimension ---
    if wasFreqDomain
        sdim = dims.f;  if sdim == 0, sdim = dims.t; end
    else
        sdim = dims.t;                 % freq lands on the time dim after the FFT
    end
    assert(sdim > 0, 'op_CSIaddPhase: no spectral/time dimension found.');

    data = MRSIStruct.data;
    if ~wasFreqDomain
        data = fftshift(fft(data, [], sdim), sdim);
    end

    % --- ppm axis ---
    if isfield(MRSIStruct,'ppm') && ~isempty(MRSIStruct.ppm)
        ppm = MRSIStruct.ppm(:);
    else
        N   = size(data, sdim);  sw = MRSIStruct.spectralWidth;  df = sw/N;
        ppm = ((-sw/2 : df : sw/2 - df).') / (MRSIStruct.Bo * MRSIStruct.gamma);
    end
    assert(numel(ppm) == size(data, sdim), ...
        'op_CSIaddPhase: ppm length %d ~= spectral dim %d', numel(ppm), size(data, sdim));

    % --- Hz per ppm: prefer txfrq (reliable), else Bo*gamma ---
    if isfield(MRSIStruct,'txfrq') && ~isempty(MRSIStruct.txfrq) && MRSIStruct.txfrq > 0
        hzppm = MRSIStruct.txfrq / 1e6;
    else
        hzppm = MRSIStruct.Bo * MRSIStruct.gamma;
    end

    % --- phase ramp: ph0 (deg) + first-order 2*pi*ph1*(f - fpivot) ---
    phase = ph0 * pi/180 + 2*pi*ph1 * (ppm - ppm0) * hzppm;   % [N x 1]
    sh        = ones(1, ndims(data));
    sh(sdim)  = numel(phase);
    dataPhased = data .* reshape(exp(1i*phase), sh);

    % --- optional before/after plot at a central column ---
    if ~suppressPlot
        toCol = @(x) reshape(permute(x, [sdim, setdiff(1:ndims(x), sdim)]), size(x, sdim), []);
        B = toCol(data);  A = toCol(dataPhased);  c = max(round(size(B,2)/2), 1);
        figure; plot(ppm, real(B(:,c)), 'b', ppm, real(A(:,c)), 'r--');
        set(gca,'XDir','reverse'); grid on;
        legend('before','after'); xlabel('ppm'); ylabel('real');
        title(sprintf('op_CSIaddPhase: ph0=%g deg, ph1=%g s (pivot %.2f ppm)', ph0, ph1, ppm0));
    end

    % --- back to time domain if the input was time domain ---
    if ~wasFreqDomain
        dataPhased = ifft(ifftshift(dataPhased, sdim), [], sdim);
    end
    MRSIStruct.data = dataPhased;
end
