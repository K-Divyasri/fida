function tf = is_mrsi_data(p)
    % Determine if data is MRSI (signal OR metabolite - both 4D and 3D formats)
    tf = false;
    nVol = size(p.nii.img, 4);
    fprintf('=== is_mrsi_data detection ===\n');
    fprintf('Volume count: %d\n', nVol);
    
    % Get description
    desc = '';
    if isfield(p.nii.hdr, 'descrip')
        desc = strtrim(char(p.nii.hdr.descrip));
    end
    
    % METHOD 1: 4D metabolite format [Conc, CRLB, LW, SNR]
    if nVol == 4
        fprintf('4 volumes - checking for metabolite format...\n');
        
        % Check NIfTI description for metabolite indicators
        if contains(desc, 'CRLB', 'IgnoreCase', true) || ...
           contains(desc, 'LW', 'IgnoreCase', true) || ...
           contains(desc, 'SNR', 'IgnoreCase', true) || ...
           contains(desc, '+', 'IgnoreCase', true)  % Format: "NAA+CRLB+LW+SNR"
            fprintf('✓ 4D Metabolite data detected (from description)\n');
            tf = true;
            return;
        end
        
        % Fallback: Check for quality.mat or ftSpec (old format)
        [fpath, fname, ~] = fileparts(p.nii.hdr.file_name);
        if endsWith(fname, '.nii'), fname = extractBefore(fname, '.nii'); end
        if endsWith(fname, '.gz'), fname = extractBefore(fname, '.gz'); end
        
        quality_file = fullfile(fpath, [fname '_quality.mat']);
        if exist(quality_file, 'file')
            fprintf('✓ 4D Metabolite data detected (from quality.mat)\n');
            tf = true;
            return;
        end
    end
    
    % METHOD 2: 3D metabolite maps (NEW - from create_separate_metabolite_niftis_v2)
    if nVol == 1
        fprintf('Single volume - checking for 3D metabolite...\n');
        
        % Check description for 3D metabolite indicators
        if contains(desc, 'CRLB', 'IgnoreCase', true) || ...
           contains(desc, 'Linewidth', 'IgnoreCase', true) || ...
           contains(desc, 'SNR', 'IgnoreCase', true) || ...
           contains(desc, '_CRLB', 'IgnoreCase', true)  % Format: "NAA_CRLB [32x32x1]"
            fprintf('✓ 3D Metabolite map detected (CRLB/LW/SNR)\n');
            tf = true;
            return;
        end
        
        % Check filename patterns for 3D metabolite maps
        [fpath, fname, ~] = fileparts(p.nii.hdr.file_name);
        if endsWith(fname, '.nii'), fname = extractBefore(fname, '.nii'); end
        if endsWith(fname, '.gz'), fname = extractBefore(fname, '.gz'); end
        
        % Pattern: metabolite_CRLB_3D.nii or Linewidth_3D.nii or SNR_3D.nii
        if contains(fname, '_CRLB_3D', 'IgnoreCase', true) || ...
           contains(fname, 'Linewidth_3D', 'IgnoreCase', true) || ...
           contains(fname, 'SNR_3D', 'IgnoreCase', true)
            fprintf('✓ 3D Metabolite map detected (from filename)\n');
            tf = true;
            return;
        end
        
        % Check for quality.mat (old 3D metabolite format)
        quality_file = fullfile(fpath, [fname '_quality.mat']);
        if exist(quality_file, 'file')
            fprintf('✓ 3D Metabolite data detected (from quality.mat)\n');
            tf = true;
            return;
        end
        
        % Check for ftSpec file (signal data)
        ftSpec_file = fullfile(fpath, [fname '_ftSpec.mat']);
        if exist(ftSpec_file, 'file')
            fprintf('✓ MRSI signal data detected (from ftSpec)\n');
            tf = true;
            return;
        end
    end
    
    % METHOD 3: Multi-volume signal data (original behavior)
    if nVol >= 2
        fprintf('Multi-volume signal data\n');
        
        % Check intent code
        if any(p.nii.hdr.intent_code == [2003, 2001, 2002])
            fprintf('✓ MRSI from intent code\n');
            tf = true;
            return;
        end
        
        % Check for ppm field
        if isfield(p, 'ppm') || isfield(p.nii, 'ppm')
            fprintf('✓ MRSI from ppm field\n');
            tf = true;
            return;
        end
        
        % Check dimensions
        dim = p.nii.hdr.dim(2:4);
        if nVol > 64 && all(dim <= 128)
            fprintf('✓ MRSI from dimensions\n');
            tf = true;
            return;
        end
        
        % Check temporal unit
        temporal_unit = bitand(p.nii.hdr.xyzt_units, 56);
        if any(temporal_unit == [32, 40])
            fprintf('✓ MRSI from temporal unit\n');
            tf = true;
            return;
        end
    end
    
    fprintf('Not MRSI data\n');
end