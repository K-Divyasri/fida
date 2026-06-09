%op_CSIShift.m
%Andres Melani, Jamie Near, Sunnybrook Research Institute, 2023
%
%USAGE:
%[MRSIStruct, weightMatrix] = op_CSIShift(obj_in, delta_x, delta_y, k_file)
%
%DESCRIPTION:
%Takes a MRSI struct on k-space and a k-space file and applies a spatial
%shift 
%function [MRSIStruct, weightMatrix] = op_CSIShift(obj_in, delta_x, delta_y, k_file)

    fileName='.\dataForAndres\Human\spinecho_rosette\meas_MID01202_FID45793_fireRosette_lowRes_TE30.dat';
    fileName_w='.\dataForAndres\Human\spinecho_rosette\meas_MID01203_FID45794_fireRosette_lowRes_TE30_w.dat';
    kFile='.\dataForAndres\RosetteKspace_48x48_fov240_sw1587.txt';
    %[kSpaceX, kSpaceY] = getKSpaceReshaped(kFile);
    %[kSpX, kSpY] = getKSpaceReshaped_v2(kFile);
    
    %Load data
    rose = io_CSIload_twix(fileName);
    rose_w = io_CSIload_twix(fileName_w);
   
    %Input correct ADC dwell time (**This needs to be automated):
    rose.adcDwellTime = 5e-6; 
    rose_w.adcDwellTime = 5e-6;
    
    %Combine the ADC blocks
    timeCombined = op_CSICombineTime(rose, 'extras');
    timeCombined_w = op_CSICombineTime(rose_w, 'extras');
    clear rose
    clear rose_w

    %Re-calculate the adcTime vector, using the correct dwell time (**This
    %needs to be automated). 
    timeCombined.adcTime = timeCombined.adcDwellTime*[0:(timeCombined.sz(1)-1)];
    timeCombined_w.adcTime = timeCombined_w.adcDwellTime*[0:(timeCombined_w.sz(1)-1)];

    %Do density compensation
    dComp = op_CSIDensityCompensation(timeCombined,kFile);
    dComp_w = op_CSIDensityCompensation(timeCombined_w,kFile);
    clear timeCombined
    clear timeCombined_w

    %Spatial shift on k-space

    delta_x=0;
    delta_y=-80;

    dComp = op_CSIShift(dComp,delta_x,delta_y,kFile);
    dComp_w = op_CSIShift(dComp_w,delta_x,delta_y,kFile);

%     shifted = dComp.data .* exp(-1i*2*pi*(kSpaceX * delta_x + kSpaceY * delta_y));
%     shifted_w = dComp_w.data .* exp(-1i*2*pi*(kSpX * delta_x + kSpY * delta_y));
%     dComp.data = shifted;
%     dComp_w.data = shifted_w;

    ftSpatial = op_CSIFourierTransform(dComp, kFile, 'spatial', true, 'spectral', false);
    ftSpatial_w = op_CSIFourierTransform(dComp_w, kFile, 'spatial', true, 'spectral', false);

    [coilCombined_w, phase, weights] = op_CSICombineCoils(ftSpatial_w);
    [coilCombined] = op_CSICombineCoils(ftSpatial, 1, phase, weights);
    clear ftSpatial
    clear ftSpatial_w

    ccav = op_CSIAverage(coilCombined);
    ccav_w = op_CSIAverage(coilCombined_w);

    ftSpec = op_CSIFourierTransform(ccav);
    ftSpec_w = op_CSIFourierTransform(ccav_w);

% function [kSpX, kSpY] = getKSpaceReshaped(filename)
%     % Read the data from the CSV file into a matrix
%     kSpaceT = readtable(filename);
%     kx = kSpaceT.Kx;
%     ky = kSpaceT.Ky;
%     %Organize in 126 points per 76 TR
%     kx=reshape(kx,126,1,1,76);
%     ky=reshape(ky,126,1,1,76);
%     %Replicate the points to match dComp data size
%     kSpX=repmat(kx,576,16,6);
%     kSpY=repmat(ky,576,16,6);
%     %Add next dimension (16 coil channels)
% %     kx=repmat(kx,1,16);
% %     ky=repmat(ky,1,16);
% %     kx=reshape(kx,72576,16,76);
% %     ky=reshape(ky,72576,16,76);
% %     %Add next dimension (6 avgs)
% %     kx=repmat(kx,1,6);
% %     ky=repmat(ky,1,6);
% %     %Reshape to match k-space signal data
% %     kSpX=reshape(kx,72576,16,6,76);
% %     kSpY=reshape(ky,72576,16,6,76);
% end
% 
% function [kSpX, kSpY] = getKSpaceReshaped_v2(filename)
%     % Read the data from the CSV file into a matrix
%     kSpaceT = readtable(filename);
%     kx = kSpaceT.Kx;
%     ky = kSpaceT.Ky;
%     %Organize in 126 points per 76 TR
%     kx=reshape(kx,126,1,76);
%     ky=reshape(ky,126,1,76);
%     %Replicate the points to match dComp data size
%     kSpX=repmat(kx,576,16);
%     kSpY=repmat(ky,576,16);
%     
% end