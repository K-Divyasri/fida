function MRSIStruct = op_CSISegment_simple(MRSIStruct, threshold)
% OP_CSISEGMENT_SIMPLE  Simple mean-threshold segmentation for MRSI water data.
%
% USAGE:
%   ccav_w = op_CSISegment_simple(ccav_w);           % auto threshold (mean)
%   ccav_w = op_CSISegment_simple(ccav_w, 0.35);     % manual thresholded val
%
% OUTPUT:(same as b4)
%   MRSIStruct.mask.brainmask  - [32x32] logical mask

    %data3d    = squeeze(MRSIStruct.data);
    intensity = abs(squeeze(MRSIStruct.data(1,:,:)));
    %squeeze(max(abs(data3d), [], 1));   % [32x32] peak per voxel

    if nargin < 2 || isempty(threshold)
        threshold = mean(intensity(:));
        fprintf('Auto threshold (mean): %.2f\n', threshold);
    end

    mask = bwareaopen(intensity > threshold, 3); % should i make it 2/1?
    %notmask = bwareaopen(intensity < threshold, 3);
    mask = imfill(mask,"holes");
    MRSIStruct.mask.brainmasks = logical(mask);

    fprintf('%d / %d voxels selected (%.1f%%)\n', ...
        sum(mask(:)), numel(mask), 100*mean(mask(:)));

    % Visualise
    figure;
    subplot(1,3,1); imagesc(intensity);         colorbar; axis image off; title('Intensity');
    subplot(1,3,2); imagesc(mask);              colormap(gca,'gray'); axis image off; title(sprintf('Mask (thresh=%.2f)', threshold));
    subplot(1,3,3); imagesc(intensity .* mask); colorbar; axis image off; title('Masked');
    sgtitle('Segmentation thingy');
end