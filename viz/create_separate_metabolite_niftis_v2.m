function out = create_separate_metabolite_niftis_v2(ftSpec, map, crlb, LW, SNR, ...
                                                 t1_path, spec2nii_empty_path, output_dir)
% CREATE_SEPARATE_METABOLITE_NIFTIS_V2 - Create metabolite NIfTIs with separate quality maps
%
% Creates:
%   - One 4D NIfTI per metabolite: [Metabolite_conc, CRLB_for_that_met, LW, SNR]
%   - One 3D CRLB NIfTI per metabolite: [CRLB_for_that_met only]
%   - One 3D LW map (shared): [LW only]
%   - One 3D SNR map (shared): [SNR only]
%
% USAGE:
%   out = create_separate_metabolite_niftis_v2(ftSpec, map, crlb, LW, SNR, ...
%                                              't1.nii', 'empty.nii', 'output_folder');

    fprintf('\n========================================================================\n');
    fprintf('=== CREATING METABOLITE NIFTIS (4D + 3D QUALITY MAPS) ===\n');
    fprintf('========================================================================\n\n');
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
        fprintf('Created output directory: %s\n\n', output_dir);
    end
    
    %% Load spec2nii empty NIfTI (with correct geometry)
    fprintf('Loading spec2nii empty NIfTI...\n');
    nii_empty = nii_tool('load', spec2nii_empty_path);
    
    X_expected = nii_empty.hdr.dim(2);
    Y_expected = nii_empty.hdr.dim(3);
    Z_expected = nii_empty.hdr.dim(4);
    
    fprintf('Empty NIfTI spatial dimensions: [%d x %d x %d]\n\n', X_expected, Y_expected, Z_expected);
    
    %% Get metabolite names
    met_names = fieldnames(map);
    n_metabolites = length(met_names);
    fprintf('Found %d metabolites: %s\n\n', n_metabolites, strjoin(met_names, ', '));
    
    %% DEBUG: Check input data types
    fprintf('=== DEBUG: Input data types ===\n');
    fprintf('map type: %s\n', class(map));
    if isstruct(map)
        fprintf('  map fields: %s\n', strjoin(fieldnames(map), ', '));
    end
    
    fprintf('\ncrlb type: %s\n', class(crlb));
    fprintf('crlb size: [%s]\n', num2str(size(crlb)));
    fprintf('crlb ndims: %d\n', ndims(crlb));
    
    if isstruct(crlb)
        crlb_fields = fieldnames(crlb);
        fprintf('  crlb is STRUCT with fields: %s\n', strjoin(crlb_fields, ', '));
    elseif iscell(crlb)
        fprintf('  crlb is CELL ARRAY with %d elements\n', length(crlb));
        if ~isempty(crlb)
            fprintf('  crlb{1} size: [%s], type: %s\n', num2str(size(crlb{1})), class(crlb{1}));
        end
    else
        fprintf('  crlb is NUMERIC ARRAY\n');
    end
    
    fprintf('\nLW size: [%s]\n', num2str(size(LW)));
    fprintf('SNR size: [%s]\n\n', num2str(size(SNR)));
    
    %% Prepare shared quality maps (LW and SNR) - 3D ONLY
    fprintf('Preparing shared quality maps (3D)...\n');
    
    % LW (shared across all metabolites) - 3D
    vol_lw = LW';
    vol_lw = flipud(vol_lw);
    vol_lw = reshape(vol_lw, [X_expected, Y_expected, Z_expected]);
    fprintf('  LW map: [%d x %d x %d], range [%.3f, %.3f] ppm\n', ...
        size(vol_lw,1), size(vol_lw,2), size(vol_lw,3), min(vol_lw(:)), max(vol_lw(:)));
    
    % SNR (shared across all metabolites) - 3D
    vol_snr = SNR';
    vol_snr = flipud(vol_snr);
    vol_snr = reshape(vol_snr, [X_expected, Y_expected, Z_expected]);
    fprintf('  SNR map: [%d x %d x %d], range [%.1f, %.1f]\n', ...
        size(vol_snr,1), size(vol_snr,2), size(vol_snr,3), min(vol_snr(:)), max(vol_snr(:)));
    
    fprintf('\n');
    
    %% Initialize output structure
    out = struct();
    out.metabolite_4d_files = {};
    out.crlb_3d_files = {};
    out.shared_lw_3d_file = '';
    out.shared_snr_3d_file = '';
    out.ftSpec_file = '';
    
    %% Create NIfTI for each metabolite
    fprintf('Creating metabolite-specific NIfTIs...\n');
    fprintf('========================================================================\n');
    
    for m = 1:n_metabolites
        met_name = met_names{m};
        fprintf('\n[%d/%d] Processing %s...\n', m, n_metabolites, met_name);
        
        try
            % Get data for this metabolite
            met_data = map.(met_name);
            [mr, mc] = size(met_data);
            fprintf('  Metabolite data size: [%d x %d]\n', mr, mc);
            
            % Get CRLB data - handle different formats
            crlb_data = [];
            
            if isstruct(crlb)
                % Case 1: crlb is a structure like map
                fprintf('  Attempting to extract from STRUCT using field: %s\n', met_name);
                if isfield(crlb, met_name)
                    crlb_data = crlb.(met_name);
                    fprintf('  ✓ Found CRLB.%s with size [%s]\n', met_name, num2str(size(crlb_data)));
                else
                    fprintf('  ✗ Field %s not found in CRLB struct\n', met_name);
                    fprintf('  Available fields: %s\n', strjoin(fieldnames(crlb), ', '));
                    error('CRLB structure missing field: %s', met_name);
                end
                
            elseif iscell(crlb)
                % Case 2: crlb is a cell array - use index
                fprintf('  Attempting to extract from CELL ARRAY at index: %d\n', m);
                if m <= length(crlb)
                    crlb_data = crlb{m};
                    fprintf('  ✓ Extracted crlb{%d} with size [%s]\n', m, num2str(size(crlb_data)));
                else
                    error('CRLB cell index %d exceeds length %d', m, length(crlb));
                end
                
            elseif isnumeric(crlb)
                % Case 3: crlb is a numeric array
                crlb_size = size(crlb);
                crlb_ndims_val = ndims(crlb);
                
                fprintf('  CRLB is NUMERIC: size=[%s], ndims=%d\n', num2str(crlb_size), crlb_ndims_val);
                
                if crlb_ndims_val == 2
                    % 2D array [Y x X] - same for all
                    fprintf('  Using 2D CRLB map (same for all metabolites)\n');
                    crlb_data = crlb;
                    
                elseif crlb_ndims_val == 3
                    % 3D array [Y x X x n_metabolites]
                    fprintf('  Extracting slice %d from 3D array (3rd dim size=%d)\n', m, crlb_size(3));
                    if m <= crlb_size(3)
                        crlb_data = crlb(:,:,m);
                        fprintf('  ✓ Extracted crlb(:,:,%d) with size [%s]\n', m, num2str(size(crlb_data)));
                    else
                        error('Index %d exceeds 3rd dimension %d', m, crlb_size(3));
                    end
                    
                else
                    error('Unexpected CRLB dimensions: %d', crlb_ndims_val);
                end
                
            else
                error('Unknown CRLB type: %s', class(crlb));
            end
            
            % Verify crlb_data
            if isempty(crlb_data)
                error('Failed to extract CRLB for %s', met_name);
            end
            
            fprintf('  Final CRLB data size: [%s]\n', num2str(size(crlb_data)));
            
            % Verify dimensions
            [crlb_r, crlb_c] = size(crlb_data);
            if crlb_r ~= mr || crlb_c ~= mc
                warning('CRLB size [%d x %d] ≠ metabolite size [%d x %d]', crlb_r, crlb_c, mr, mc);
            end
            
            if X_expected ~= mc || Y_expected ~= mr
                warning('Expected [%d x %d], got [%d x %d]. Skipping %s', ...
                    Y_expected, X_expected, mr, mc, met_name);
                continue;
            end
            
            % Transform concentration
            vol_conc = met_data';
            vol_conc = flipud(vol_conc);
            vol_conc = reshape(vol_conc, [X_expected, Y_expected, Z_expected]);
            fprintf('  Concentration range: [%.3e, %.3e]\n', min(vol_conc(:)), max(vol_conc(:)));
            
            % Transform CRLB
            vol_crlb = crlb_data';
            vol_crlb = flipud(vol_crlb);
            vol_crlb = reshape(vol_crlb, [X_expected, Y_expected, Z_expected]);
            fprintf('  CRLB range: [%.1f, %.1f]%%\n', min(vol_crlb(:)), max(vol_crlb(:)));
            
            %% Create 4D metabolite file
            data_4d = cat(4, vol_conc, vol_crlb, vol_lw, vol_snr);
            data_4d = single(data_4d);
            
            nii_4d = nii_empty;
            nii_4d.hdr.dim(1) = int16(4);
            nii_4d.hdr.dim(5) = int16(4);
            nii_4d.hdr.pixdim(5) = single(1);
            nii_4d.img = data_4d;
            nii_4d.hdr.descrip = sprintf('%s+CRLB+LW+SNR [%dx%dx%dx4]', ...
                met_name, X_expected, Y_expected, Z_expected);
            
            met_file = fullfile(output_dir, sprintf('%s_4D.nii', lower(met_name)));
            nii_tool('save', nii_4d, met_file);
            fprintf('  ✓ Saved 4D: %s\n', met_file);
            out.metabolite_4d_files{end+1} = met_file;
            
            %% Create 3D CRLB file
            data_3d = single(vol_crlb);
            
            nii_3d = nii_empty;
            nii_3d.hdr.dim(1) = int16(3);
            nii_3d.hdr.dim(5) = int16(1);
            nii_3d.hdr.pixdim(5) = single(0);
            nii_3d.img = data_3d;
            nii_3d.hdr.descrip = sprintf('%s_CRLB [%dx%dx%d]', ...
                met_name, X_expected, Y_expected, Z_expected);
            
            crlb_file = fullfile(output_dir, sprintf('%s_CRLB_3D.nii', lower(met_name)));
            nii_tool('save', nii_3d, crlb_file);
            fprintf('  ✓ Saved 3D CRLB: %s\n', crlb_file);
            out.crlb_3d_files{end+1} = crlb_file;
            
        catch ME
            fprintf('  ✗ ERROR processing %s: %s\n', met_name, ME.message);
            fprintf('     at line %d in %s\n', ME.stack(1).line, ME.stack(1).name);
            continue;
        end
    end
    
    %% Create shared LW 3D
    fprintf('\n========================================================================\n');
    fprintf('Creating shared Linewidth 3D...\n');
    
    nii_lw = nii_empty;
    nii_lw.hdr.dim(1) = int16(3);
    nii_lw.hdr.dim(5) = int16(1);
    nii_lw.hdr.pixdim(5) = single(0);
    nii_lw.img = single(vol_lw);
    nii_lw.hdr.descrip = sprintf('Linewidth [%dx%dx%d]', X_expected, Y_expected, Z_expected);
    
    lw_file = fullfile(output_dir, 'Linewidth_3D.nii');
    nii_tool('save', nii_lw, lw_file);
    fprintf('✓ Saved: %s\n', lw_file);
    out.shared_lw_3d_file = lw_file;
    
    %% Create shared SNR 3D
    fprintf('\nCreating shared SNR 3D...\n');
    
    nii_snr = nii_empty;
    nii_snr.hdr.dim(1) = int16(3);
    nii_snr.hdr.dim(5) = int16(1);
    nii_snr.hdr.pixdim(5) = single(0);
    nii_snr.img = single(vol_snr);
    nii_snr.hdr.descrip = sprintf('SNR [%dx%dx%d]', X_expected, Y_expected, Z_expected);
    
    snr_file = fullfile(output_dir, 'SNR_3D.nii');
    nii_tool('save', nii_snr, snr_file);
    fprintf('✓ Saved: %s\n', snr_file);
    out.shared_snr_3d_file = snr_file;
    
    %% Save ftSpec
    fprintf('\nSaving ftSpec...\n');
    ftSpec_path = fullfile(output_dir, 'ftSpec.mat');
    ftSpec_smooth_w = ftSpec;
    save(ftSpec_path, 'ftSpec_smooth_w', '-v7.3');
    fprintf('✓ Saved: %s\n', ftSpec_path);
    out.ftSpec_file = ftSpec_path;
    
    %% Summary
    fprintf('\n========================================================================\n');
    fprintf('=== CREATION COMPLETE ===\n');
    fprintf('========================================================================\n');
    fprintf('Output directory: %s\n\n', output_dir);
    
    fprintf('Metabolite 4D files (%d):\n', length(out.metabolite_4d_files));
    for i = 1:length(out.metabolite_4d_files)
        [~, fn, ext] = fileparts(out.metabolite_4d_files{i});
        fprintf('  %d. %s%s\n', i, fn, ext);
    end
    
    fprintf('\nCRLB 3D files (%d):\n', length(out.crlb_3d_files));
    for i = 1:length(out.crlb_3d_files)
        [~, fn, ext] = fileparts(out.crlb_3d_files{i});
        fprintf('  %d. %s%s\n', i, fn, ext);
    end
    
    fprintf('\nShared maps:\n');
    [~, fn, ext] = fileparts(out.shared_lw_3d_file);
    fprintf('  - %s%s (LW)\n', fn, ext);
    [~, fn, ext] = fileparts(out.shared_snr_3d_file);
    fprintf('  - %s%s (SNR)\n', fn, ext);
    
    fprintf('\nTotal: %d files\n', length(out.metabolite_4d_files) + length(out.crlb_3d_files) + 3);
    fprintf('========================================================================\n\n');
end