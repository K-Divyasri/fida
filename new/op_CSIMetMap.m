%%Creating metabolite maps
%Change these parameters

% op_CSIMetMap.m
% Lubna Burki, SunnyBrook Hosptial 2025.
%
% USAGE:
% [met_map,crlb_map]=op_CSIMetMap(ftSpec,lcmodel_file_name);
% 
% DESCRIPTION:
% Creation of metabolite maps based on LCM fitting
% 
% INPUTS:
% ftSpec = FID-A CSI structure
% lcmodel_file_name = Name of LCModel file string
%
% OUTPUTS:
% met_map = metabolite map of given metabolite
% crlb_map = CRLB map of a given metabolite


function [met_map,crlb_map] = op_CSIMetMap(metabolite,MRSIStructspectral,lcmodel_name,options)
    arguments
        metabolite (1,:) char;
        MRSIStructspectral (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        lcmodel_name (1,:) char
        options.isPlotMask (1,1) logical = false;
        options.MRSIStructspatial (1,1) struct %For mask generation, PRE-LIPID REMOVAL!
        % phaseMap double = []
        % weightMap double = []
        
    end
%%Add csv


file_name = [lcmodel_name, '.csv'];
file = readtable(file_name);
%file = readtable('lcm_24x24_te8p5_uniform_ssp.csv');
%%Initialize metabolite and CRLB maps
x_dim = MRSIStructspectral.sz(MRSIStructspectral.dims.x);
y_dim = MRSIStructspectral.sz(MRSIStructspectral.dims.y);

met_map = zeros(x_dim,y_dim);
crlb_map = zeros(x_dim,y_dim);

%Identify metabolites of interest
%metabolites = "Ala Asp PCh Cr PCr GABA Gln Glu Gly Ins Lac NAA Tau Glc NAAG GPC PE Ser Asc";
%metabolites = "NAA";




metabolite_sd = [metabolite, '_SD'];


if options.isPlotMask == true
    brain_area_raw = createBrainArea(MRSIStructspectral,x_dim/2,y_dim/2,isPlotMask=options.isPlotMask,MRSIStructspatial=options.MRSIStructspatial);
    brain_area = imfill(MRSIStructspectral.mask.brainmasks,'holes');
    brain_area(brain_area == 0) = brain_area_raw(brain_area==0);
else
    [brain_area_raw,~] = createBrainArea(MRSIStructspectral,x_dim/2,y_dim/2);
    brain_area = imfill(MRSIStructspectral.mask.brainmasks,'holes');
    brain_area(brain_area == 0) = brain_area_raw(brain_area==0);
end


for i = 1:size(file.(metabolite),1)
    met_map(file.Col(i),file.Row(i)) = file.(metabolite)(i);
    met_map = met_map.*brain_area;
    crlb_map(file.Col(i),file.Row(i)) = file.(metabolite_sd)(i);
    crlb_map = crlb_map.*brain_area;
end

%crlb_rs = reshape(crlb_map,[1,size(crlb_map,1)*size(crlb_map,2)]);
%disp(["Percent CRLB that passes for", string(metabolite), 'is', num2str((size(crlb_rs(crlb_rs > 0 & crlb_rs < 20)) / size(crlb_rs(crlb_rs > 0)))*100)]) 

end