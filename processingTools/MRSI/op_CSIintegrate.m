%op_CSIIntegrate.m
% Brenden Kadota, Sunnybrook 2021.
%
% USAGE:
% in=op_CSIintegrate(in);
%
% DESCRIPTION:
% Integrates the spectrum from ppmmin to ppmmax in all csi voxels.
%
% INPUTS:
% MRSIStruct        = CSI FID-A data structure
% pmmmin    = lower integration bounds
% ppmmax    = upper integration bounds
% mode      = mode (optional):
%                        -'re' (integral performed on real part (default)).
%                        -'im' (integral performed on imaginary part).
%                        -'mag' (integral performed on magnitude part).
%
% OUTPUTS:
% map       = map of integrated area
function map = op_CSIintegrate(MRSIStruct, ppmmin, ppmmax, mode)
    arguments
        MRSIStruct (1, 1) struct
        ppmmin (1, 1) double
        ppmmax (1, 1) double
        mode (1, :) char {mustBeMember(mode, {'re', 'im', 'mag'})} = 're'
    end
    % check MRSI Struct
    checkArguments(MRSIStruct);
    
    % resahpe to time, y and x dimensions
    MRSIStruct = reshapeDimensions(MRSIStruct, {'f', 'y', 'x'});
    % intalize map size
    %map = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));

    nY      = getSizeFromDimensions(MRSIStruct, {'y'});   % 1-based Y
    nX      = getSizeFromDimensions(MRSIStruct, {'x'});   % 1-based X
    nExtras = getSizeFromDimensions(MRSIStruct, {'extras'});
    if nExtras == 0, nExtras = 1; end                     % treat “no extras” as 1
    
    map = zeros(nY, nX, nExtras);
    
    hasExtras = (getDimension(MRSIStruct,'extras') > 0);
    
    for e = 1:nExtras
        for y = 1:nY          % FIRST loop over Y  (dims.y)
            for x = 1:nX      % then   loop over X  (dims.x)
                if hasExtras
                    voxel = op_CSItoMRS(MRSIStruct, x, y, struct('extraIndex',e));
                else
                    voxel = op_CSItoMRS(MRSIStruct, x, y);
                end
                map(y, x, e) = op_integrate(voxel, ppmmin, ppmmax, mode);
            end
        end
    end

end

% argument checks. Check if spectral and spatial fourier transform have been done.
function checkArguments(in)
    if(in.flags.spectralFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spectral dimension');
    end
    if(in.flags.spatialFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spacial dimension');
    end
end
