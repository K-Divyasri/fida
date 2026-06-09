function fix_mrsi_alignment_targeted(ftSpec_smooth, t1_filename, output_filename)
% Fix MRSI alignment based on the specific affine matrix and T1 properties
% Your MRSI affine shows centering at (-120, -120, -7.5) which needs adjustment
%
% Usage: fix_mrsi_alignment_targeted(ftSpec_smooth, 'T1.nii.gz', 'mrsi_fixed.nii.gz')

    fprintf('=== TARGETED MRSI ALIGNMENT FIX ===\n');
    
    % Load T1 reference
    t1_nii = load_nii(t1_filename);
    t1_size = size(t1_nii.img);
    t1_voxel = t1_nii.hdr.dime.pixdim(2:4);
    
    fprintf('T1 image: %dx%dx%d, voxel size: %.2fx%.2fx%.2f mm\n', ...
            t1_size, t1_voxel);
    
    % Your MRSI info
    mrsi_affine = ftSpec_smooth.affineMatrix;
    mrsi_voxel = [ftSpec_smooth.voxelSize.x, ftSpec_smooth.voxelSize.y, ftSpec_smooth.voxelSize.z];
    mrsi_size = [size(ftSpec_smooth.data, 2), size(ftSpec_smooth.data, 3), 1];
    
    fprintf('MRSI: %dx%dx%d, voxel size: %.1fx%.1fx%.1f mm\n', ...
            mrsi_size, mrsi_voxel);
    fprintf('Original MRSI origin: [%.1f, %.1f, %.1f]\n', ...
            mrsi_affine(1:3, 4));
    
    % Calculate T1 image center in mm (assuming no spatial transform)
    % T1 center in voxel coordinates
    t1_center_vox = t1_size / 2;
    
    % Convert to mm (assuming image is centered at origin)
    % This is a reasonable assumption for brain images
    t1_fov = t1_size .* t1_voxel;
    t1_center_mm = -t1_fov / 2 + t1_center_vox .* t1_voxel;
    
    fprintf('T1 center (estimated): [%.1f, %.1f, %.1f] mm\n', t1_center_mm);
    
    % MRSI field of view
    mrsi_fov = mrsi_size(1:2) .* mrsi_voxel(1:2);
    fprintf('MRSI FOV: %.1f x %.1f mm\n', mrsi_fov);
    
    % Create corrected affine matrix
    % Center MRSI on T1 center
    corrected_origin = [t1_center_mm(1) - mrsi_fov(1)/2;
                        t1_center_mm(2) - mrsi_fov(2)/2;
                        t1_center_mm(3)]; % Keep original z-position or adjust as needed
    
    corrected_affine = mrsi_affine;
    corrected_affine(1:3, 4) = corrected_origin;
    
    fprintf('Corrected MRSI origin: [%.1f, %.1f, %.1f] mm\n', corrected_origin);
    
    % Create the NIfTI
    create_mrsi_with_custom_affine(ftSpec_smooth, corrected_affine, output_filename);
    
    % Also create versions with different z-positions for testing
    test_z_positions = [-20, -10, 0, 10, 20]; % mm relative to T1 center
    
    fprintf('\nCreating test versions with different z-positions:\n');
    for i = 1:length(test_z_positions)
        test_z = test_z_positions(i);
        test_affine = corrected_affine;
        test_affine(3, 4) = t1_center_mm(3) + test_z;
        
        test_filename = sprintf('mrsi_test_z%+d.nii.gz', test_z);
        create_mrsi_with_custom_affine(ftSpec_smooth, test_affine, test_filename);
        
        fprintf('  %s (z = %.1f mm)\n', test_filename, test_affine(3,4));
    end
    
    fprintf('\nTest each version:\n');
    fprintf('nii_viewer(''%s'', ''%s''); %% Main corrected version\n', t1_filename, output_filename);
    for i = 1:length(test_z_positions)
        test_z = test_z_positions(i);
        test_filename = sprintf('mrsi_test_z%+d.nii.gz', test_z);
        fprintf('nii_viewer(''%s'', ''%s'');\n', t1_filename, test_filename);
    end
end

function create_mrsi_with_custom_affine(ftSpec_smooth, affine_matrix, output_filename)
% Create MRSI NIfTI with custom affine matrix
    
    % Extract spectroscopic data
    spec_data = ftSpec_smooth.data; % [576 x 40 x 40]
    [n_freq, n_x, n_y] = size(spec_data);
    
    % Create 4D data
    nifti_data_4d = zeros(n_x, n_y, 1, n_freq, 'single');
    for freq_idx = 1:n_freq
        nifti_data_4d(:, :, 1, freq_idx) = squeeze(spec_data(freq_idx, :, :))';
    end
    
    % Create NIfTI
    voxel_x = ftSpec_smooth.voxelSize.x;
    voxel_y = ftSpec_smooth.voxelSize.y;
    voxel_z = ftSpec_smooth.voxelSize.z;
    
    nii = make_nii(nifti_data_4d, [voxel_x, voxel_y, voxel_z]);
    
    % Fix 4D header
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(2:5) = [n_x, n_y, 1, n_freq];
    nii.hdr.dime.dim(6:8) = [1, 1, 1];
    nii.hdr.dime.pixdim(5) = ftSpec_smooth.spectralDwellTime;
    nii.hdr.dime.datatype = 16;
    nii.hdr.dime.bitpix = 32;
    
    % Set custom affine transformation
    nii.hdr.hist.srow_x = affine_matrix(1,:);
    nii.hdr.hist.srow_y = affine_matrix(2,:);
    nii.hdr.hist.srow_z = affine_matrix(3,:);
    nii.hdr.hist.sform_code = 1;
    
    % Save
    save_nii(nii, output_filename);
end

function quick_manual_fix(ftSpec_smooth, t1_filename, x_shift, y_shift, z_shift, output_filename)
% Quick manual adjustment from the original affine matrix
% Positive shifts move: x=right, y=anterior, z=superior
%
% Usage: quick_manual_fix(ftSpec_smooth, 'T1.nii.gz', 120, 120, 20, 'mrsi_manual.nii.gz')

    original_affine = ftSpec_smooth.affineMatrix;
    
    % Apply shifts to the origin
    adjusted_affine = original_affine;
    adjusted_affine(1, 4) = original_affine(1, 4) + x_shift;
    adjusted_affine(2, 4) = original_affine(2, 4) + y_shift;
    adjusted_affine(3, 4) = original_affine(3, 4) + z_shift;
    
    fprintf('Original origin: [%.1f, %.1f, %.1f]\n', original_affine(1:3, 4));
    fprintf('Adjusted origin: [%.1f, %.1f, %.1f]\n', adjusted_affine(1:3, 4));
    
    create_mrsi_with_custom_affine(ftSpec_smooth, adjusted_affine, output_filename);
    
    fprintf('Manual adjustment saved: %s\n', output_filename);
    fprintf('Test with: nii_viewer(''%s'', ''%s'');\n', t1_filename, output_filename);
end

function analyze_coordinate_system(ftSpec_smooth, t1_filename)
% Analyze the coordinate systems to understand the alignment issue
    
    fprintf('=== COORDINATE SYSTEM ANALYSIS ===\n');
    
    % MRSI coordinate analysis
    mrsi_coords = ftSpec_smooth.coordinates;
    mrsi_x_range = [min(mrsi_coords.x), max(mrsi_coords.x)];
    mrsi_y_range = [min(mrsi_coords.y), max(mrsi_coords.y)];
    
    fprintf('MRSI coordinate ranges:\n');
    fprintf('  X: %.1f to %.1f mm (center: %.1f)\n', mrsi_x_range, mean(mrsi_x_range));
    fprintf('  Y: %.1f to %.1f mm (center: %.1f)\n', mrsi_y_range, mean(mrsi_y_range));
    
    % Affine matrix analysis
    affine = ftSpec_smooth.affineMatrix;
    fprintf('MRSI affine matrix origin: [%.1f, %.1f, %.1f]\n', affine(1:3, 4));
    
    % Calculate where MRSI grid actually sits
    mrsi_size = [size(ftSpec_smooth.data, 2), size(ftSpec_smooth.data, 3)];
    mrsi_voxel = [ftSpec_smooth.voxelSize.x, ftSpec_smooth.voxelSize.y];
    
    % First voxel position
    first_voxel_pos = affine(1:3, 4);
    
    % Last voxel position  
    last_voxel_pos = first_voxel_pos + [mrsi_size(1) * mrsi_voxel(1); 
                                        mrsi_size(2) * mrsi_voxel(2); 
                                        0];
    
    % Center position
    center_pos = first_voxel_pos + [mrsi_size(1) * mrsi_voxel(1) / 2;
                                    mrsi_size(2) * mrsi_voxel(2) / 2;
                                    0];
    
    fprintf('MRSI grid positions:\n');
    fprintf('  First voxel: [%.1f, %.1f, %.1f]\n', first_voxel_pos);
    fprintf('  Center: [%.1f, %.1f, %.1f]\n', center_pos);
    fprintf('  Last voxel: [%.1f, %.1f, %.1f]\n', last_voxel_pos);
    
    % T1 analysis
    t1_nii = load_nii(t1_filename);
    t1_size = size(t1_nii.img);
    t1_voxel = t1_nii.hdr.dime.pixdim(2:4);
    t1_fov = t1_size .* t1_voxel;
    
    fprintf('T1 image FOV: %.1f x %.1f x %.1f mm\n', t1_fov);
    fprintf('T1 assumed center: [%.1f, %.1f, %.1f] mm\n', -t1_fov/2 + t1_fov/2);
    
    fprintf('=============================\n');
end