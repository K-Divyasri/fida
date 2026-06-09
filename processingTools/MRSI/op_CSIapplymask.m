function MRSIStruct = op_CSIapplymask(MRSIStruct)

fprintf('=== DEBUG: op_CSIapplymask ===\n');

if ~isfield(MRSIStruct,'mask')
    error("No mask found. Please segment");
end

% === DETAILED SIZE DEBUGGING ===
fprintf('MRSIStruct.data size: [%s]\n', num2str(size(MRSIStruct.data)));
fprintf('MRSIStruct.sz: [%s]\n', num2str(MRSIStruct.sz));

if isempty(MRSIStruct.mask.brainmasks)
    warning('Mask is empty - skipping masking operation');
    return;
end

fprintf('brainmasks size: [%s]\n', num2str(size(MRSIStruct.mask.brainmasks)));
fprintf('brainmasks class: %s\n', class(MRSIStruct.mask.brainmasks));
fprintf('brainmasks sum: %d voxels\n', sum(MRSIStruct.mask.brainmasks(:)));

% === ATTEMPT TO BUILD THE MASK ===
try
    % Replicate mask along spectral dimension
    mask_replicated = repmat(MRSIStruct.mask.brainmasks, 1, 1, MRSIStruct.sz(1));
    fprintf('Replicated mask size: [%s]\n', num2str(size(mask_replicated)));
    
    % Permute to match data dimensions
    mask_permuted = permute(mask_replicated, [3 1 2]);
    fprintf('Permuted mask size: [%s]\n', num2str(size(mask_permuted)));
    
    % Apply mask
    fprintf('Attempting element-wise multiplication...\n');
    MRSIStruct.data = MRSIStruct.data .* mask_permuted;
    fprintf('✓ Masking successful\n');
    
catch ME
    fprintf('✗ ERROR during masking:\n');
    fprintf('  Message: %s\n', ME.message);
    fprintf('  Expected data size: [%d, %d, %d]\n', MRSIStruct.sz);
    fprintf('  Actual data size: [%s]\n', num2str(size(MRSIStruct.data)));
    fprintf('  Mask size: [%s]\n', num2str(size(MRSIStruct.mask.brainmasks)));
    rethrow(ME);
end

end