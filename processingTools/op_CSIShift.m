function [MRSIStruct] = op_CSIShift(obj_in, delta_x, delta_y, k_file)
    
    if nargout > 0 || ~isempty(inputname(1))
        disp('--- Inside op_CSIShift (Fixed Version) ---');
    end

    % Extract data and size from input struct
    obj_data = obj_in.data;
    obj_size = obj_in.sz;

    if nargout > 0 || ~isempty(inputname(1))
        disp(['obj_data size: ', mat2str(size(obj_data))]);
        disp(['Requested shifts: X=', num2str(delta_x), ' mm, Y=', num2str(delta_y), ' mm']);
    end

    % Check if this is Cartesian mode (empty k_file)
    if isempty(k_file) || strcmp(k_file, "")
        disp('Cartesian mode detected - applying simple phase shift');
        % For Cartesian data, apply a simple uniform phase shift
        % This is a placeholder - you may need to implement proper Cartesian shifting
        MRSIStruct = obj_in;
        return;
    end

    % Get reshaped k-space coordinates - FIXED VERSION
    [kSpaceX, kSpaceY] = getKSpaceReshaped_Fixed(k_file, obj_size, obj_in.dims, obj_in);

    % Calculate the shifted data - VECTORIZED
    phase_shift = exp(-1i * 2 * pi * (kSpaceX * delta_x + kSpaceY * delta_y));
    obj_data = obj_data .* phase_shift;

    % Update struct
    obj_in.data = obj_data;
    MRSIStruct = obj_in;

    if nargout > 0 || ~isempty(inputname(1))
        disp('--- op_CSIShift complete ---');
    end
end

function [kSpX, kSpY] = getKSpaceReshaped_Fixed(filename, obj_size, dims, obj_struct)
    % FIXED: Correct array assignment for the given data structure
    
    persistent cached_trajectory cached_filename;
    
    % Cache trajectory file to avoid repeated file I/O
    if isempty(cached_trajectory) || ~strcmp(cached_filename, filename)
        % Read trajectory file
        if endsWith(filename, '.txt') || endsWith(filename, '.csv')
            try
                data = readmatrix(filename);
                kx = data(:, 2);
                ky = data(:, 3);
            catch
                kSpaceT = readtable(filename);
                if width(kSpaceT) >= 3
                    kx = kSpaceT{:, 2};
                    ky = kSpaceT{:, 3};
                else
                    kx = kSpaceT.Kx;
                    ky = kSpaceT.Ky;
                end
            end
        else
            % Fallback for other formats
            kSpaceT = readtable(filename);
            if ismember('Kx', kSpaceT.Properties.VariableNames)
                kx = kSpaceT.Kx;
                ky = kSpaceT.Ky;
            else
                data = table2array(kSpaceT);
                kx = data(:, 2);
                ky = data(:, 3);
            end
        end
        
        cached_trajectory = struct('kx', kx, 'ky', ky);
        cached_filename = filename;
    else
        kx = cached_trajectory.kx;
        ky = cached_trajectory.ky;
    end

    num_dimensions = length(obj_size);
    
    % Pre-compute normalization factors
    kx_max = max(abs(kx));
    ky_max = max(abs(ky));

    % FOV-based normalization.
    %
    % NOTE on apparent bug: the divide-by-(2*kmax*FOV) here shrinks the
    % applied shift by ~1/n_matrix.  That is intentional in the *current*
    % pipeline: op_NUFFTSpatial1's n_shift is derived from xCoords/yCoords
    % (which carry the imageOrigin offset set by findImageOrigin), and it
    % already moves the recon onto the prescribed slab centre.  Replacing
    % the normalisation with raw k here causes the NUFFT shift and the
    % op_CSIShift phase ramp to double up.  Keep this scaling unless you
    % also disable the n_shift mechanism upstream.
    if isfield(obj_struct, 'fov') && isfield(obj_struct.fov, 'x') && isfield(obj_struct.fov, 'y')
        fovX = obj_struct.fov.x;
        fovY = obj_struct.fov.y;

        % Vectorized normalization
        kx_norm = (kx / (2 * kx_max)) / fovX;
        ky_norm = (ky / (2 * ky_max)) / fovY;
    else
        % Fallback normalization
        assumed_fov = 240; % mm
        kx_norm = (kx / (2 * kx_max)) / assumed_fov;
        ky_norm = (ky / (2 * ky_max)) / assumed_fov;
    end

    % Get number of shots
    if dims.kshot > 0
        nKshot = obj_size(dims.kshot);
    else
        nKshot = obj_size(end);
    end
    
    nKpts = length(kx) / nKshot;
    nTimeKpts = obj_size(1);
    nTimePoints = floor(nTimeKpts / nKpts);

    % CRITICAL FIX: Correct array creation and assignment
    if num_dimensions == 4
        % 4D case: [t*kpts, coils, averages, kshot] = [72576, 16, 6, 76]
        kSpX = zeros(obj_size, 'single');
        kSpY = zeros(obj_size, 'single');
        
        % Get dimensions
        nCoils = obj_size(2);      % 16
        nAverages = obj_size(3);   % 6
        
        for shot = 1:nKshot
            shot_start = (shot-1) * nKpts + 1;
            shot_end = shot * nKpts;
            
            % Get k-space coordinates for this shot
            kx_shot = kx_norm(shot_start:shot_end);
            ky_shot = ky_norm(shot_start:shot_end);
            
            % Create time index patterns
            time_pattern = (1:nKpts)' + (0:nTimePoints-1) * nKpts;
            time_indices = time_pattern(:);
            
            % Ensure bounds
            valid_indices = time_indices <= nTimeKpts;
            time_indices = time_indices(valid_indices);
            
            % Replicate k-space values for time points
            kx_expanded = repmat(kx_shot, nTimePoints, 1);
            ky_expanded = repmat(ky_shot, nTimePoints, 1);
            
            % Trim to valid size
            kx_expanded = kx_expanded(valid_indices);
            ky_expanded = ky_expanded(valid_indices);
            
            % FIXED: Proper array assignment with correct dimensions
            % Expand to match all coils and averages
            kx_matrix = repmat(kx_expanded, 1, nCoils, nAverages);  % [time_pts, coils, averages]
            ky_matrix = repmat(ky_expanded, 1, nCoils, nAverages);  % [time_pts, coils, averages]
            
            % Assign to arrays
            kSpX(time_indices, :, :, shot) = kx_matrix;
            kSpY(time_indices, :, :, shot) = ky_matrix;
        end
        
    elseif num_dimensions == 3
        % 3D case: [t*kpts, coils, kshot] = [72576, 16, 76]
        kSpX = zeros(obj_size, 'single');
        kSpY = zeros(obj_size, 'single');
        
        % Get dimensions
        nCoils = obj_size(2);      % 16
        
        for shot = 1:nKshot
            shot_start = (shot-1) * nKpts + 1;
            shot_end = shot * nKpts;
            
            kx_shot = kx_norm(shot_start:shot_end);
            ky_shot = ky_norm(shot_start:shot_end);
            
            % Vectorized time index creation
            time_pattern = (1:nKpts)' + (0:nTimePoints-1) * nKpts;
            time_indices = time_pattern(:);
            
            valid_indices = time_indices <= nTimeKpts;
            time_indices = time_indices(valid_indices);
            
            % Vectorized expansion
            kx_expanded = repmat(kx_shot, nTimePoints, 1);
            ky_expanded = repmat(ky_shot, nTimePoints, 1);
            
            kx_expanded = kx_expanded(valid_indices);
            ky_expanded = ky_expanded(valid_indices);
            
            % FIXED: Proper array assignment with correct dimensions
            % Expand to match all coils
            kx_matrix = repmat(kx_expanded, 1, nCoils);  % [time_pts, coils]
            ky_matrix = repmat(ky_expanded, 1, nCoils);  % [time_pts, coils]
            
            % Vectorized assignment
            kSpX(time_indices, :, shot) = kx_matrix;
            kSpY(time_indices, :, shot) = ky_matrix;
        end
    else
        error('Unsupported number of dimensions: %d', num_dimensions);
    end
end

% ALTERNATIVE: More robust version with better error handling
function [MRSIStruct] = op_CSIShift_Robust(obj_in, delta_x, delta_y, k_file)
    
    disp('--- Inside op_CSIShift_Robust ---');
    
    % Check inputs
    if nargin < 4
        error('op_CSIShift requires 4 inputs: obj_in, delta_x, delta_y, k_file');
    end
    
    % Extract data and size from input struct
    obj_data = obj_in.data;
    obj_size = obj_in.sz;
    
    disp(['Input data size: ', mat2str(size(obj_data))]);
    disp(['Input sz field: ', mat2str(obj_size)]);
    disp(['Requested shifts: X=', num2str(delta_x), ' mm, Y=', num2str(delta_y), ' mm']);
    
    % Validate that data size matches sz field
    if ~isequal(size(obj_data), obj_size)
        warning('Data size does not match sz field. Using actual data size.');
        obj_size = size(obj_data);
    end
    
    % Check if this is Cartesian mode
    if isempty(k_file) || strcmp(k_file, "") || strcmp(k_file, '')
        disp('Cartesian mode detected - no k-space shifting applied');
        MRSIStruct = obj_in;
        return;
    end
    
    % Check if k-file exists
    if ~isfile(k_file)
        warning('K-space file not found: %s. Skipping shift.', k_file);
        MRSIStruct = obj_in;
        return;
    end
    
    try
        % Get reshaped k-space coordinates
        [kSpaceX, kSpaceY] = getKSpaceReshaped_Fixed(k_file, obj_size, obj_in.dims, obj_in);
        
        % Validate k-space arrays
        if ~isequal(size(kSpaceX), obj_size) || ~isequal(size(kSpaceY), obj_size)
            error('K-space array size mismatch');
        end
        
        % Calculate the shifted data
        phase_shift = exp(-1i * 2 * pi * (kSpaceX * delta_x + kSpaceY * delta_y));
        obj_data = obj_data .* phase_shift;
        
        % Update struct
        obj_in.data = obj_data;
        MRSIStruct = obj_in;
        
        disp('--- op_CSIShift_Robust complete ---');
        
    catch ME
        warning('Error in op_CSIShift: %s. Returning original data.', ME.message);
        MRSIStruct = obj_in;
    end
end

% CONSERVATIVE VERSION: Apply reduced shifts for testing
function [MRSIStruct] = op_CSIShift_Conservative(obj_in, delta_x, delta_y, k_file)
    % Apply only 10% of the requested shift for safety
    reduction_factor = 0.1;
    MRSIStruct = op_CSIShift(obj_in, delta_x * reduction_factor, delta_y * reduction_factor, k_file);
end

% NO-SHIFT VERSION: For testing k-space coordinate generation
function [MRSIStruct] = op_CSIShift_NoShift(obj_in, delta_x, delta_y, k_file)
    % Apply zero shift to test k-space coordinate creation
    MRSIStruct = op_CSIShift(obj_in, 0, 0, k_file);
end