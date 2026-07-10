%op_CSIRemoveWater.m
%Brenden Kadota, Sunnybrook 2021.
%
% USAGE:
% in=op_CSIRemoveWater(in);
%
% DESCRIPTION:
% This function uses op_removeWater to remove water signal from spectrum.
% op_removeWater ruses the HSVD method to remove the water which is
% described by H. BARKHUIJSEN et al. 1987.
%
% INPUTS:
% in        = CSI FID-A data structure
%
% OUTPUTS:
% out       = output
function MRSIStruct = op_CSIRemoveWater(MRSIStruct, wlim, kInitial, M)
    arguments
        MRSIStruct (1, 1) struct
        wlim (1, 2) double = [4.4 5]
        kInitial (1, 1) double = 30
        M (1, 1) double = floor(MRSIStruct.sz(MRSIStruct.dims.f)*.75);
    end
    residual_error = 0;
    [MRSIStruct, prevPermute, prevShape]= reshapeDimensions(MRSIStruct, {'f', 'y', 'x'});
    data = zeros(getSizeFromDimensions(MRSIStruct, {'f', 'y', 'x', 'extras'}));
    % threshold to skip empty/noise voxels (HSVD blows up on near-zero fids)
    allData = getData(MRSIStruct);
    voxThr = 0.02 * max(abs(allData(:)));
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                voxel = op_CSItoMRS(MRSIStruct, x, y, struct('extraIndex', e));
                if max(abs(voxel.specs(:))) < voxThr    % no real signal -> pass through
                    data(:, y, x, e) = voxel.specs;
                    continue;
                end
                voxel_supressed = op_removeWater(voxel, wlim, kInitial, M);
                data(:, y, x, e) = voxel_supressed.specs;
                residual_error = residual_error + voxel_supressed.watersupp.residual_error;
            end
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);


    MRSIStruct.watersupp.residual_error = residual_error;
end
