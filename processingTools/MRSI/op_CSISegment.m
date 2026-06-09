
function MRSIStruct = op_CSISegment(MRSIStruct, blob, method, threshold_mult)
% op_CSISegment - Robust segmentation for phantom and brain MRSI data
%
% USAGE:
%   For phantom (solid object):     op_CSISegment(data)
%   For brain (lipid ring):         op_CSISegment(data, 30, 'ring')
%   Stricter threshold:             op_CSISegment(data, 30, 'otsu', 1.5)
%
% INPUTS:
%   MRSIStruct      : MRSI data structure
%   blob            : Min component size to keep (default = 30)
%   method          : 'otsu' (default), 'ring', 'adaptive', 'manual'
%   threshold_mult  : Threshold multiplier (default = 1.0)
%                     >1.0 = stricter (smaller mask, brighter pixels only)
%                     <1.0 = more lenient (larger mask)
%
% OUTPUTS:
%   MRSIStruct.mask.brainmasks : Signal region mask
%   MRSIStruct.mask.lipmasks   : Lipid ring (only for 'ring' method)

arguments
    MRSIStruct (1,1) struct
    blob (1,1) double = 30
    method (1,1) string = "otsu"
    threshold_mult (1,1) double = 1.0
end

%% Extract magnitude image
ccav = abs(squeeze(MRSIStruct.data(1,:,:)));

fprintf('\n=== CSI Segmentation ===\n');
fprintf('Method: %s\n', method);
fprintf('Data size: %d x %d\n', size(ccav,1), size(ccav,2));
fprintf('Range: [%.4f, %.4f], Mean: %.4f\n', min(ccav(:)), max(ccav(:)), mean(ccav(:)));

% Smooth
I = imgaussfilt(ccav, 1.5);

%% Segmentation based on method
switch lower(method)
    case "otsu"
        [signal_mask, lipid_mask] = segment_otsu(I, blob, threshold_mult);

    case "ring"
        [signal_mask, lipid_mask] = segment_ring(I, blob, threshold_mult);

    case "adaptive"
        [signal_mask, lipid_mask] = segment_adaptive(I, blob, threshold_mult);

    case "manual"
        signal_mask = segment_manual(ccav);
        lipid_mask = [];

    otherwise
        error('Method must be: otsu, ring, adaptive, or manual');
end

%% Check if segmentation succeeded
if isempty(signal_mask) || sum(signal_mask(:)) == 0
    error('Segmentation failed! Try:\n  - Different method\n  - Lower threshold_mult (e.g., 0.8)\n  - Smaller blob size\n  - Manual mode');
end

%% Visualization
fprintf('\n=== Results ===\n');
fprintf('Signal region: %d pixels (%.1f%% of FOV)\n', ...
        sum(signal_mask(:)), 100*sum(signal_mask(:))/numel(signal_mask));
if ~isempty(lipid_mask)
    fprintf('Lipid ring: %d pixels (%.1f%% of FOV)\n', ...
            sum(lipid_mask(:)), 100*sum(lipid_mask(:))/numel(lipid_mask));
end

% Create figure
figure('Name', 'Segmentation Results', 'Position', [100 100 1200 400]);

if isempty(lipid_mask)
    % Simple segmentation (phantom mode)
    subplot(1,4,1), imagesc(ccav), axis image, colorbar, title('Original');
    subplot(1,4,2), imagesc(I), axis image, colorbar, title('Smoothed');
    subplot(1,4,3), imshow(signal_mask), axis image, title('Signal Mask');
    subplot(1,4,4), imagesc(ccav .* signal_mask), axis image, colorbar, title('Masked Signal');
else
    % Ring segmentation (brain mode)
    subplot(2,3,1), imagesc(ccav), axis image, colorbar, title('Original');
    subplot(2,3,2), imagesc(I), axis image, colorbar, title('Smoothed');
    subplot(2,3,3), imshow(lipid_mask), axis image, title('Lipid Ring');
    subplot(2,3,4), imshow(signal_mask), axis image, title('Brain (Interior)');
    subplot(2,3,5), imagesc(ccav .* signal_mask), axis image, colorbar, title('Brain Signal');
    subplot(2,3,6), imagesc(ccav .* lipid_mask), axis image, colorbar, title('Lipid Signal');
end

%% Package output
mask.brainmasks = logical(signal_mask);
mask.lipmasks = logical(lipid_mask);
mask.method = method;
mask.threshold_mult = threshold_mult;
mask.numthresh = 1;
mask.thresh = [];

MRSIStruct.mask = mask;
fprintf('\n=== Ready to apply with op_CSIapplymask() ===\n\n');
end

%% ========== HELPER FUNCTIONS ==========

function [signal_mask, lipid_mask] = segment_otsu(I, blob, threshold_mult)
% Simple Otsu thresholding - for solid objects (phantoms)
fprintf('Using Otsu thresholding...\n');

level = graythresh(I) * threshold_mult;
level = min(level, 1.0);
fprintf('  Threshold level: %.4f\n', level);

BW = imbinarize(I, level);
BW = cleanup_mask(BW, blob);

signal_mask = get_largest_component(BW);
signal_mask = imfill(signal_mask, 'holes');
lipid_mask = [];
end

function [signal_mask, lipid_mask] = segment_ring(I, blob, threshold_mult)
% Ring-based segmentation - for brain data with lipid ring
fprintf('Using ring segmentation (brain mode)...\n');

% Find the outer ring
level = graythresh(I) * threshold_mult * 0.7;  % Lower threshold for ring
level = min(level, 1.0);
fprintf('  Ring threshold: %.4f\n', level);

BW = imbinarize(I, level);
BW = cleanup_mask(BW, blob);

% Get largest component (should be ring + interior)
outer_mask = get_largest_component(BW);

% Fill to get total region
filled_mask = imfill(outer_mask, 'holes');

% Brain is the filled interior minus the outer ring
signal_mask = filled_mask & ~outer_mask;

% If brain mask is too small, use filled mask instead
if sum(signal_mask(:)) < blob * 2
    fprintf('  Warning: Ring detection failed, using filled mask\n');
    signal_mask = filled_mask;
    lipid_mask = [];
else
    lipid_mask = outer_mask;
end
end

function [signal_mask, lipid_mask] = segment_adaptive(I, blob, threshold_mult)
% Adaptive thresholding with multiple sensitivity values
fprintf('Using adaptive thresholding...\n');

best_mask = [];
best_area = 0;

sens_range = (0.3:0.1:0.7) * threshold_mult;
sens_range = max(0.1, min(sens_range, 0.9));

for sens = sens_range
    thresh = adaptthresh(I, sens);
    BW = imbinarize(I, thresh);
    BW = cleanup_mask(BW, blob);

    CC = bwconncomp(BW);
    if CC.NumObjects > 0
        stats = regionprops(CC, 'Area');
        [max_area, idx] = max([stats.Area]);

        if max_area > best_area
            best_area = max_area;
            best_mask = false(size(BW));
            best_mask(CC.PixelIdxList{idx}) = true;
        end
    end
end

if isempty(best_mask)
    warning('Adaptive failed, falling back to Otsu');
    [signal_mask, lipid_mask] = segment_otsu(I, blob, threshold_mult);
else
    signal_mask = imfill(best_mask, 'holes');
    lipid_mask = [];
    fprintf('  Found mask with area: %d\n', sum(signal_mask(:)));
end
end

function signal_mask = segment_manual(ccav)
% Manual polygon drawing
fprintf('Manual segmentation mode...\n');
figure;
imagesc(ccav), axis image, colorbar;
title('Draw polygon around signal region (double-click when done)');
roi = drawpolygon('Color', 'r', 'LineWidth', 2);
signal_mask = createMask(roi);
close(gcf);
end

function BW = cleanup_mask(BW, blob)
% Morphological cleanup
BW = imclose(BW, strel('disk', 2));
BW = imerode(BW, strel('disk', 1));
BW = imdilate(BW, strel('disk', 1));
BW = bwareaopen(BW, blob);
end

function mask = get_largest_component(BW)
% Extract largest connected component
CC = bwconncomp(BW);
if CC.NumObjects == 0
    mask = false(size(BW));
    return;
end

stats = regionprops(CC, 'Area');
[~, idx] = max([stats.Area]);

mask = false(size(BW));
mask(CC.PixelIdxList{idx}) = true;
end