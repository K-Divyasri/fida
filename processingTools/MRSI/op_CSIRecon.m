function ftSpatial = op_CSIRecon(MRSIStruct, kFile, dcfMethod, ftMethod, varargin)
%OP_CSIRECON  CSI reconstruction = density compensation + spatial FT.
%
%   ftSpatial = op_CSIRecon(MRSIStruct, kFile, dcfMethod, ftMethod)
%   ftSpatial = op_CSIRecon(MRSIStruct, kFile, dcfMethod, ftMethod, Name, Value, ...)
%
%   dcfMethod : 'nn' | 'voronoi' | 'pipe_menon' | 'none'
%   ftMethod  : 'nufft' | 'dft' | 'tikhonov'
%
%   When ftMethod is 'tikhonov', density compensation is skipped — the
%   Tikhonov-regularised inverse already absorbs the role of the DCF.  Pass
%   dcfMethod='none' to skip explicitly with any FT method.
%
%   CARTESIAN DATA: if the trajectory file is empty OR the struct carries
%   flags.isCartesian == true, the data is reconstructed with a plain
%   FFT (op_CSIFourierTransform).  There is no non-uniform sampling, so
%   density compensation and NUFFT/DFT/Tikhonov are all skipped — matching
%   the manual Cartesian path used in cart.m.
%
%   Extra Name/Value pairs are forwarded to the DCF backend (e.g.
%   'numIterations' for pipe_menon, 'modelType' for voronoi, etc.).
%
%   Returns the spatially-reconstructed MRSI struct.

    if nargin < 2, kFile = ''; end
    if nargin < 3 || isempty(dcfMethod), dcfMethod = 'nn';    end
    if nargin < 4 || isempty(ftMethod),  ftMethod  = 'nufft'; end

    % --- Cartesian detection: empty k-file, or loader flagged it Cartesian ---
    kfChar = '';
    if ~isempty(kFile), kfChar = char(kFile); end
    isCartesian = isempty(kfChar) || ...
        (isfield(MRSIStruct,'flags') && isfield(MRSIStruct.flags,'isCartesian') ...
         && MRSIStruct.flags.isCartesian);

    if isCartesian
        fprintf('=== op_CSIRecon: Cartesian -> FFT (no DCF, no NUFFT) ===\n');
        ftSpatial = op_CSIFourierTransform(MRSIStruct, "", ...
                        'spatial', true, 'spectral', false);
        return;
    end

    skipDcf = strcmpi(ftMethod, 'tikhonov') || strcmpi(ftMethod, 'tikh') ...
              || strcmpi(dcfMethod, 'none');

    if skipDcf
        fprintf('=== op_CSIRecon: ft=%s  (DCF skipped) ===\n', ftMethod);
        ftSpatial = op_CSIReconstruct(MRSIStruct, kFile, ftMethod);
    else
        fprintf('=== op_CSIRecon: dcf=%s  ft=%s ===\n', dcfMethod, ftMethod);
        dComp     = op_CSIPSFCorrection(MRSIStruct, kFile, dcfMethod, varargin{:});
        ftSpatial = op_CSIReconstruct(dComp, kFile, ftMethod);
    end
end
