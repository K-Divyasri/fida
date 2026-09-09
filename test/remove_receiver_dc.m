function [obj, dc, ratio] = remove_receiver_dc(obj, kFile, varargin)
%REMOVE_RECEIVER_DC  Estimate + remove the receiver DC offset from raw rosette k-space.
%   A constant complex bias per coil (ADC/receiver DC) is uniform across all k, so
%   the rosette k=0 revisit concentrates it into the FOV-CENTER voxel -> center
%   dropout / inverted water / LCModel CRLB 999. This estimates it from the
%   high-|k| samples (the object's k-space energy is near k=0; out at high |k| it
%   averages to ~0 over thousands of samples, so the surviving mean is the DC).
%
%   [obj,dc,ratio] = remove_receiver_dc(obj, kFile, 'frac',0.8, 'apply',true)
%     frac  : use |k| > frac*max(|k|) as the DC-estimation region (default 0.8)
%     apply : subtract it (true) or just report (false)
%
%   Prints |DC|/|k0 signal| per coil. If that is tiny (<~2%), DC is NOT your
%   center problem -> it is coil sensitivity / B1 at the center of a large phantom
%   (not a recon fix). If it is several %+, DC was real and is now removed.

    ip = inputParser;
    ip.addParameter('frac', 0.8);
    ip.addParameter('apply', true);
    ip.parse(varargin{:});
    frac = ip.Results.frac;  doApply = ip.Results.apply;

    % --- |k| per sample from the k-file (kpt-fastest ordering) ---
    T = readtable(kFile); vn = lower(string(T.Properties.VariableNames));
    ix = find(vn=="kx",1); iy = find(vn=="ky",1);
    if isempty(ix)||isempty(iy), ix=2; iy=3; end
    kr = hypot(double(T{:,ix}), double(T{:,iy}));   % [Nk x 1]
    Nk = numel(kr);

    % --- data + coils dim ---
    d  = obj.data;
    s  = size(d);
    try, cDim = getDimension(obj,'coils'); catch, cDim = 2; end
    assert(cDim==2, 'expected coils in dim 2, got %d - adjust reshape.', cDim);

    % merge trailing k-dims into one Nk axis: [t coils rest Nk]
    assert(mod(prod(s(3:end)), Nk)==0 && prod(s(end-1:end))==Nk || prod(s(3:end))/Nk==round(prod(s(3:end))/Nk), ...
           'cannot factor k dims to Nk=%d from size %s', Nk, mat2str(s));
    rest = prod(s(3:end)) / Nk;
    d2   = reshape(d, [s(1) s(2) rest Nk]);          % [t coils rest Nk]

    hiK = kr > frac*max(kr);
    loK = kr < 0.1*max(kr);
    dc  = squeeze(mean(d2(:,:,:,hiK), [1 3 4]));      % [coils x 1] receiver DC
    sig = squeeze(mean(abs(d2(:,:,:,loK)),[1 3 4]));  % [coils x 1] near-k0 magnitude
    ratio = abs(dc(:)) ./ max(sig(:), eps);

    fprintf('\n[remove_receiver_dc] frac=%.2f  Nk=%d  hiK=%d loK=%d samples\n', ...
            frac, Nk, nnz(hiK), nnz(loK));
    fprintf('  |DC|/|k0 signal| per coil (median=%.3f max=%.3f):\n    ', ...
            median(ratio), max(ratio));
    fprintf('%.3f ', ratio); fprintf('\n');
    if median(ratio) < 0.02
        fprintf('  -> DC is small. NOT the center problem (look at coil sensitivity / B1).\n');
    else
        fprintf('  -> DC is significant. Removing it should fill the FOV-center dropout.\n');
    end

    if doApply
        dcShape = ones(1, numel(s)); dcShape(2) = s(2);
        obj.data = d - reshape(dc, dcShape);          % broadcast over t/avg/k
        fprintf('  DC subtracted from all samples.\n');
    end
end
