%%Creates LW Maps

% op_CSIlcm_lwmap.m
% Lubna Burki, SunnyBrook Hosptial 2025.
%
% USAGE:
% [LW]=op_CSIlw_map_lcm(ftSpec,lcmodel_file_name);
% 
% DESCRIPTION:
% Creation of LW map based on LCM output
% 
% INPUTS:
% ftSpec = FID-A CSI structure
% lcmodel_file_name = Name of LCModel file string
%
% OUTPUTS:
% LW = map of LW in relevant area


function [LW] = op_CSIlw_map_lcm(MRSIStructspectral,lcmodel_name,options)
    arguments
        MRSIStructspectral (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        lcmodel_name (1,:) char {mustBeFile}
        options.isPlotMask (1,1) logical = false %Plot brain mask
        options.MRSIStructspatial (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        options.figure_folder_name (1,:) char = "false"
        % phaseMap double = []
        % weightMap double = []
        
    end

%size_map = 24;
size_map= MRSIStructspectral.sz(MRSIStructspectral.dims.x);
%cd ftSpec_smooth_24x24_8p5_lcm_out/
%lcmodel_name = 'ftSpec_smooth_24x24_8p5_lcm';
lcmodel_directory = [lcmodel_name,'_out'];
cd(lcmodel_directory);

%Rest will run
a = dir('*.coord');
size_folder = numel(a);

LW=zeros(size_map);

for n=1:size_folder
    text = fileread(a(n).name);
    LW_integer = strfind(text,'FWHM');
    
    coordinate_oi = strcat(lcmodel_name,'_sl1_');
    redoa = erase(a(n).name,coordinate_oi);
    redoa = erase(redoa,'.coord');
    redoa = split(redoa,'-');
    x=str2num(redoa{1});
    y=str2num(redoa{2});
    if size(str2num(text(LW_integer+6:LW_integer+12))) ==0
        LW(y,x) = 0;
    else
        LW(y,x) = str2num(text(LW_integer+6:LW_integer+12));
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
    imagesc(brain_area)
end
LW = LW .* brain_area;
% 
% disp("Percent LW map that passes is")
LW_rs = reshape(LW,[1,size(LW,1)*size(LW,2)]);
% disp( ((size(LW_rs(LW_rs < 0.1 & LW_rs>0))) / size(LW_rs(LW_rs>0)) )*100)
LW_mean=((size(LW_rs(LW_rs < 0.1 & LW_rs>0))) / size(LW_rs(LW_rs>0)) )*100;
LW_mean = ['Percent LW that passes is ' num2str(LW_mean)];

%Final LW Map Generation
f=figure;
imagesc(LW)
colorbar('Color','w','FontSize',20)
colormap hot
clim([0 0.5]);
title("LW Map in ppm, from LCModel",'Color','w')
subtitle(LW_mean,'Color','w')
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[0 0 0]);
axis = gca;
axis.YColor = 'w';
axis.XColor = 'w';


cd ..

if(options.figure_folder_name ~= "false")
    saveas(f,fullfile(pwd,options.figure_folder_name,'LW_map_lcm'),'jpg');
end


end
