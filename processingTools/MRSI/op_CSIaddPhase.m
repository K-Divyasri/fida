function MRSIStruct = op_CSIaddPhase(MRSIStruct, ph0, ph1, ppm0, suppressPlot)
% Correct and verified phase correction (0th and 1st order) for CSI data.
arguments
    MRSIStruct (1,1) struct
    ph0 (1,1) double
    ph1 (1,1) double = 0
    ppm0 (1,1) double = 4.65
    suppressPlot (1,1) logical = true
end

% Extract dimensions and data
dims = MRSIStruct.dims;
tDim = dims.t;
data = MRSIStruct.data;

% Check if spectralFT is already done
wasFreqDomain = MRSIStruct.flags.spectralFT;

% FFT to spectral domain if needed (correct handling)
if ~wasFreqDomain
    data = fftshift(fft(data, [], tDim), tDim);
end

% Generate ppm axis if not present
if isfield(MRSIStruct, 'ppm')
    ppm = MRSIStruct.ppm;
else
    Nt = size(data, tDim);
    sw = MRSIStruct.spectralWidth;
    df = sw/Nt;
    freq = linspace(-sw/2, sw/2 - df, Nt);
    ppm = freq / (MRSIStruct.Bo * MRSIStruct.gamma);
end

% Apply zeroth and first-order phase correctly in spectral domain
ppm_shifted = ppm - ppm0;
phase0_rad = ph0 * pi/180; % degrees to radians
phase1_rad = 2 * pi * ph1 * ppm_shifted * MRSIStruct.Bo * MRSIStruct.gamma;

% Reshape for broadcasting
phase_corr = reshape(exp(1i * (phase0_rad + phase1_rad)), [], 1, 1);

% Apply the phase correction
data = bsxfun(@times, data, phase_corr);

% Inverse FFT to time domain if original data was in time domain
if ~wasFreqDomain
    data = ifft(ifftshift(data, tDim), [], tDim);
end

% Save corrected data back
MRSIStruct.data = data;

% Optional plotting for verification
if ~suppressPlot
    ix = round(size(data, dims.x)/2);
    iy = round(size(data, dims.y)/2);
    voxel_before = squeeze(MRSIStruct.data(:, iy, ix));
    voxel_after = squeeze(data(:, iy, ix));

    % Ensure frequency-domain for plotting
    if wasFreqDomain
        spec_before = voxel_before;
        spec_after = voxel_after;
    else
        spec_before = fftshift(fft(voxel_before));
        spec_after = fftshift(fft(voxel_after));
    end

    figure;
    plot(ppm, real(spec_before), 'b', 'DisplayName', 'Before Phase');
    hold on;
    plot(ppm, real(spec_after), 'r--', 'DisplayName', 'After Phase');
    set(gca, 'XDir', 'reverse');
    xlabel('ppm'); ylabel('Magnitude');
    legend;
    title(sprintf('Spectrum at voxel (%d,%d)', ix, iy));
    grid on;
end

end
