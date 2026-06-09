function MRSIStruct = op_CSIssp(MRSIStruct, minppm, maxppm, m, useEcon)
% Applies Singular-Value-based SSP lipid suppression to an MRSI structure.
%
% USAGE:
%   out = op_CSIssp(in, minppm, maxppm, m, useEcon)
%
% INPUTS
%   in         : FID-A MRSI structure (must already be spatialFT & spectralFT)
%   minppm     : lower bound of lipid ppm window
%   maxppm     : upper bound of lipid ppm window
%   m          : # SVD components to remove   (default = 6)
%   useEcon    : use economy-sized SVD        (default = false)
%
% OUTPUT
%   out        : MRSI structure with SSP applied

% -------------------------------------------------------------------------

arguments
    MRSIStruct (1,1) struct
    minppm     (1,1) double
    maxppm     (1,1) double
    m          (1,1) double  = 6
    useEcon    (1,1) logical = false
end

%% already be in spectral domain
if getFlags(MRSIStruct, "spectralFT") == 0
    error("op_CSIssp:InputError", ...
          "Data are still in the time domain; Fourier-transform first.");
end

%% 1. Identify which dimension is spectral first check f and then if not fetch t
if isfield(MRSIStruct.dims, "f") && MRSIStruct.dims.f > 0
    specLabel = "f";
elseif isfield(MRSIStruct.dims, "t") && MRSIStruct.dims.t > 0
    specLabel = "t";
else
    error("op_CSIssp:DimError", ...
          "Neither ‘f’ nor ‘t’ dimension found in MRSIStruct.dims");
end
%% same logic as prev
%% 2. Bring data to order (y, x, spec) for easy reshaping
[MRSIStruct, pastOrder, currentOrder] = reshapeDimensions(MRSIStruct, ...
                                                          {'y','x',char(specLabel)});

dori                = MRSIStruct.data;           % [ny × nx × nspec]
% Get dimensions
[ny,nx,nSpectra] = size(dori);

% Reshape into column
doricol = reshape(dori,[],nSpectra);

% Define ppm lipid bounds
lowppmlip = minppm;
highppmlip = maxppm;

%% Find the index for the ppm bounds

% Calculate absolute differences
lowdiff = abs(MRSIStruct.ppm - lowppmlip);
highdiff = abs(MRSIStruct.ppm - highppmlip);

% Find the index of the minimum difference
[~, endidx] = min(lowdiff);
[~, startidx] = min(highdiff);

%%
% Initialize dorilip matrix
numspecpoints = abs(endidx-startidx) + 1;
dorilip = zeros(ny * nx, numspecpoints);

% Extract lipid peak points at each voxel
for i = 1:ny*nx
    dorilip(i, :) = doricol(i,startidx:endidx);
end

% Perform SVD
if useEcon
    [U, ~, ~] = svd(dorilip, 'econ');
else
    [U, ~, ~] = svd(dorilip);
end

% Modify SVD matrices to only include first m spatial components
Um = U(:, 1:m);

% Calculate projection matrix
P = eye(size(Um, 1)) - Um * Um';

% Apply the projection matrix to original image
dsup = P * doricol;

% Reshape to (y,x,t)
dsup = reshape(dsup, ny, nx, nSpectra);
MRSIStruct = setData(MRSIStruct,dsup);

% Permute back into (t,y,x)
MRSIStruct = reshapeBack(MRSIStruct, pastOrder, currentOrder);
end


%learn svd 
% % op_CSIssp.m
% % Kaito Hara-Lee, SunnyBrook Hospital 2024.
% %
% % Description: Applies SSP for lipid suppression
% %
% % USAGE:
% % out = op_CSIssp(in, minppm, maxppm, m, 'econ')
% %
% % input:    in         = MRSI Structure spectral domain
% %           minppm     = minimum ppm value for lipid peak
% %           maxppm     = maximum ppm value for lipid peak
% %           m          = first m svd components (default 6)
% %           econ       = use economy-sized SVD if
% %                        MRSI data is too big (default false)
% % output:   out        = MRSI Structure with SSP applied
% 
% function MRSIStruct = op_CSIssp(MRSIStruct,minppm,maxppm,m,useEcon)
% arguments
%     MRSIStruct (1,1) struct
%     minppm (1,1) double
%     maxppm (1,1) double
%     m (1,1) double = 6;
%     useEcon (1,1) logical = false;
% end
% 
% % Check data is in spectral domain
% if (getFlags(MRSIStruct,'spectralFT') == 0)
%     disp('Data in Time domain, please Fourier Transform');
%     return;
% end
% 
% % Permute data
% [MRSIStruct,pastorder,currentorder] = reshapeDimensions(MRSIStruct ...
%     ,{'y', 'x','t'});
% 
% % Extract data
% dori = MRSIStruct.data;
% 
% % Get dimensions
% [ny,nx,nSpectra] = size(dori);
% 
% % Reshape into column
% doricol = reshape(dori,[],nSpectra);
% 
% % Define ppm lipid bounds
% lowppmlip = minppm;
% highppmlip = maxppm;
% 
% %% Find the index for the ppm bounds
% 
% % Calculate absolute differences
% lowdiff = abs(MRSIStruct.ppm - lowppmlip);
% highdiff = abs(MRSIStruct.ppm - highppmlip);
% 
% % Find the index of the minimum difference
% [~, endidx] = min(lowdiff);
% [~, startidx] = min(highdiff);
% 
% %%
% % Initialize dorilip matrix
% numspecpoints = abs(endidx-startidx) + 1;
% dorilip = zeros(ny * nx, numspecpoints);
% 
% % Extract lipid peak points at each voxel
% for i = 1:ny*nx
%     dorilip(i, :) = doricol(i,startidx:endidx);
% end
% 
% % Perform SVD
% if useEcon
%     [U, ~, ~] = svd(dorilip, 'econ');
% else
%     [U, ~, ~] = svd(dorilip);
% end
% 
% % Modify SVD matrices to only include first m spatial components
% Um = U(:, 1:m);
% 
% % Calculate projection matrix
% P = eye(size(Um, 1)) - Um * Um';
% 
% % Apply the projection matrix to original image
% dsup = P * doricol;
% 
% % Reshape to (y,x,t)
% dsup = reshape(dsup, ny, nx, nSpectra);
% MRSIStruct = setData(MRSIStruct,dsup);
% 
% % Permute back into (t,y,x)
% MRSIStruct = reshapeBack(MRSIStruct, pastorder, currentorder);
% end
