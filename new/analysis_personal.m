% %%Create spectra for central voxel
% 
% vox = op_CSItoMRS(ftSpec,24,24);
% op_plotspec(vox,0,4);
% 
% 
% format long
% lipidmap=op_CSIintegrate(ftSpec,0.89,1.41,'mag'); %-2 tp 1.5
% 
% cr_lip = op_CSILipidPerformanceMap(ftSpec,2.6);
% brain_area = createMask(ftSpec,24,24);
% cr_lip=cr_lip.*brain_area;
% crlip_mean = mean2(cr_lip(cr_lip~=0));
% disp(crlip_mean)
% 
% [snr, signal, noised] = op_CSIgetSNR(ftSpec,2.4,2.8,10,12);
% snr = snr.*brain_area;
% snr_mean = mean2(snr(snr~=0));
% disp(snr_mean)
% 
% io_CSIwritelcm(ftSpec_smooth,'ftSpec_smooth_lcm');
% io_CSIwritelcm(ftSpec_smooth_w,'ftSpec_smooth_lcm_w');

%Enter the .mat file starting in ftSpec_raw... should load into workspace
%as ftSpec
[met,crlb] = op_CSI_ManyMetMap("CrPCr GPCPCh NAANAAG GluGln Ins", ...%look at csv files and look at header naming
    ftSpec ... %ftSpec 
    ,lcm_file_name, ... %location of these coord table files - 
    figure_folder_name=figurefoldername); %can get rid of this
%lcm_file_name = [lcm_file_name '_out'];
LW = op_CSIlw_map_lcm(ftSpec,lcm_file_name,figure_folder_name=figurefoldername);
SNR = op_CSIsnr_map_lcm(ftSpec,lcm_file_name,figure_folder_name=figurefoldername);