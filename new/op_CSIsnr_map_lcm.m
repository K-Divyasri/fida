%%Creating SNR map
%Change these parameters

% op_CSIlcm_snrmap.m
% Lubna Burki, SunnyBrook Hosptial 2025.
%
% USAGE:
% [SNR]=op_CSIsnr_map_lcm(ftSpec,lcmodel_file_name);
% 
% DESCRIPTION:
% Creation of SNR map based on LCM output
% 
% INPUTS:
% ftSpec = FID-A CSI structure
% lcmodel_file_name = Name of LCModel file string
%
% OUTPUTS:
% SNR = map of SNR in relevant area

function [SNR] = op_CSIsnr_map_lcm(MRSIStructspectral,lcmodel_name,options)
    arguments
        MRSIStructspectral (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        lcmodel_name (1,:) char {mustBeFile}
        options.isPlotMask (1,1) logical = false;
        options.MRSIStructspatial (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        options.figure_folder_name (1,:) char = "false"
        % phaseMap double = []
        % weightMap double = []
        
    end
% size_map = 24;
% cd ftSpec_smooth_24x24_8p5_lcm_out/
% lcmodel_name = 'ftSpec_smooth_24x24_8p5_lcm';
size_map= MRSIStructspectral.sz(MRSIStructspectral.dims.x);
lcmodel_directory = [lcmodel_name,'_out'];
cd(lcmodel_directory);

%Rest will run
a = dir('*.coord');
size_folder = numel(a);

SNR=zeros(size_map);

for n=1:size_folder
    text = fileread(a(n).name);
    SNR_integer = strfind(text,'S/N');
    
    coordinate_oi = strcat(lcmodel_name,'_sl1_');
    redoa = erase(a(n).name,coordinate_oi);
    redoa = erase(redoa,'.coord');
    redoa = split(redoa,'-');
    x=str2num(redoa{1});
    y=str2num(redoa{2});
    if size(str2num(text(SNR_integer+5:SNR_integer+9))) ==0
        SNR(y,x) = 0;
    else
        SNR(y,x) = str2num(text(SNR_integer+5:SNR_integer+9));
    end
end

if options.isPlotMask == true
    brain_area_raw = createBrainArea(MRSIStructspectral,size_map/2,size_map/2,isPlotMask=options.isPlotMask,MRSIStructspatial=options.MRSIStructspatial);
    brain_area = imfill(MRSIStructspectral.mask.brainmasks,'holes');
    brain_area(brain_area == 0) = brain_area_raw(brain_area==0);
else
    brain_area_raw = createBrainArea(MRSIStructspectral,size_map/2,size_map/2);
    brain_area = imfill(MRSIStructspectral.mask.brainmasks,'holes');
    brain_area(brain_area == 0) = brain_area_raw(brain_area==0);
end

SNR = SNR .* brain_area;

SNR_rs = reshape(SNR,[1,size(SNR,1)*size(SNR,2)]);
SNR_mean=((size(SNR_rs(SNR_rs > 3))) / size(SNR_rs(SNR_rs>0)) )*100;
SNR_mean = ['Percent SNR that passes is ' num2str(SNR_mean)];

%Final SNR map generation
f=figure;
imagesc(SNR)
colorbar('Color','w','FontSize',20)
colormap hot
clim([0 10]);
title("SNR Map, from LCModel",'Color','w')
subtitle(SNR_mean,'Color','w')
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[0 0 0]);
axis = gca;
axis.YColor = 'w';
axis.XColor = 'w';

cd ..

if(options.figure_folder_name ~= "false")
    saveas(f,fullfile(pwd,options.figure_folder_name,'SNR_map_lcm'),'jpg');
end


end