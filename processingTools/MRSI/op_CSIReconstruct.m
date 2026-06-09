function ftSpatial = op_CSIReconstruct(dComp, kFile, method)
%OP_CSIRECONSTRUCT  Unified spatial-reconstruction dispatcher for CSI.
%
%   ftSpatial = op_CSIReconstruct(dComp, kFile, method)
%
%   method:
%     'nufft'    -> op_NUFFTSpatial1            (NUFFT, gridded)
%     'dft'      -> op_CSIFourierTransform_dft  (slow DFT, sft2 operator)
%     'tikhonov' -> op_CSIFourierTransform      (Tikhonov inverse)
%
%   Only the SPATIAL transform is performed here — spectral FT is left to
%   downstream code so each backend behaves identically at this stage.

    if nargin < 3 || isempty(method)
        method = 'nufft';
    end

    methodLower = lower(string(method));
    switch methodLower
        case "nufft"
            ftSpatial = op_NUFFTSpatial1(dComp, kFile);

        case "dft"
            ftSpatial = op_CSIFourierTransform_dft(dComp, kFile, ...
                'spatial', true, 'spectral', false);

        case {"tikhonov", "tikh"}
            ftSpatial = op_CSIFourierTransform(dComp, kFile, ...
                'spatial', true, 'spectral', false);

        otherwise
            error('op_CSIReconstruct:badMethod', ...
                'Unknown FT method "%s". Use ''nufft'', ''dft'', or ''tikhonov''.', method);
    end
end
