%op_CSILipidPerformanceMap.m
% Lubna Burki, Sunnybrook 2024.
%
% USAGE:
% in=op_CSILipidPerformanceMap(in);
%
% DESCRIPTION:
% Show ratio between metabolite of interest and lipid region
%
% INPUTS:
% MRSIStruct = CSI FID-A data structure
% ppm = approximate ppm value of metabolite of interest
% 
% OUTPUTS:
% map = map of ratios
function [map_cr_lip,map_snr] = op_CSILipidPerformanceMap(MRSIStruct, ppm)
    arguments
        MRSIStruct (1, 1) struct
        ppm (1, 1) double
    end
    % check MRSI Struct
    format long
    checkArguments(MRSIStruct);
    
    % resahpe to time, y and x dimensions
    MRSIStruct = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
    % intalize map size
    map_cr_lip = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x'}));
    map_snr = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x'}));


    
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'y'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'x'})
                warning('off')
                voxel = op_CSItoMRS(MRSIStruct, x, y, "Extra", e);
                %vox_cc=op_complexConj(voxel);
                %[out, K, wppm, amp, alpha, ph, model]=op_removeWater(voxel,[0,5]);
                [~, ~, wppm, amp, ~, ~, ~]=op_removeWater(voxel,[0,5]);
                if size(wppm)==1
                    metabolite_amp =0;
                else
                    nearest_ppm = interp1(wppm, wppm, ppm, 'nearest');
                    metabolite_index = find(wppm == nearest_ppm);
                    metabolite_amp = amp(metabolite_index);
                end
                if isnan(nearest_ppm)
                    metabolite_amp = 0;
                end
                lipid_amp = op_integrate(voxel, 0.89, 1.41, 'mag');
                map_cr_lip(y, x, e) = metabolite_amp / lipid_amp;
                map_snr(y,x,e) = op_getSNR(voxel,2.4,2.8,10,11);
            end
        end
    end

end
function checkArguments(in)
    if(in.flags.spectralFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spectral dimension');
    end
    if(in.flags.spatialFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spacial dimension');
    end
end
%end
