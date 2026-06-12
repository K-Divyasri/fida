function show_mrsi_spectrum_onclick(hs, p, c, ax)
    % Display metabolite info OR spectrum based on data type
    % Works with both 4D and 3D metabolite formats
    
    fprintf('\n=== Click Display ===\n');
    
    % Get coordinates safely
    ijk_bg_raw = get(hs.ijk, 'Value');
    if iscell(ijk_bg_raw)
        ijk_bg = cell2mat(ijk_bg_raw);
    else
        ijk_bg = ijk_bg_raw;
    end
    
    if numel(ijk_bg) == 2
        ijk_bg(3) = 1;
    elseif numel(ijk_bg) < 2
        return;
    end
    
    fprintf('Background voxel: [%d, %d, %d]\n', ijk_bg(1), ijk_bg(2), ijk_bg(3));
    
    try
        img_size_raw = size(p.nii.img);
        
        % Handle different dimension cases
        if length(img_size_raw) == 2
            % 2D image [X Y] - add Z=1 and T=1
            img_size = [img_size_raw(1), img_size_raw(2), 1, 1];
        elseif length(img_size_raw) == 3
            % 3D image [X Y Z] - add T=1
            img_size = [img_size_raw(1), img_size_raw(2), img_size_raw(3), 1];
        else
            % 4D or higher
            img_size = img_size_raw;
        end
        
        fprintf('Image size: [%s] (normalized to 4D: [%s])\n', ...
            num2str(img_size_raw), num2str(img_size));
        
        % Transform coordinates
        if isfield(p, 'R0')
            xyz = hs.bg.R * [ijk_bg(1:3)-1; 1];
            xyz = xyz(1:3);
            
            if img_size(3) == 1  % 2D
                try
                    ijk_mrsi_temp = p.Ri * [xyz; 1];
                    if numel(ijk_mrsi_temp) >= 2
                        ijk_mrsi = zeros(3, 1);
                        ijk_mrsi(1) = round(ijk_mrsi_temp(1)) + 1;
                        ijk_mrsi(2) = round(ijk_mrsi_temp(2)) + 1;
                        ijk_mrsi(3) = 1;
                    else
                        return;
                    end
                catch
                    return;
                end
            else  % 3D
                ijk_mrsi_temp = round(p.Ri * [xyz; 1]) + 1;
                if numel(ijk_mrsi_temp) >= 3
                    ijk_mrsi = ijk_mrsi_temp(1:3);
                else
                    return;
                end
            end
        else
            ijk_mrsi = ijk_bg(1:3);
        end
        
        % Panel index of the pulled ftSpec voxel (ix=a, iy=Ny-b+1) is (a, b),
        % matching the integration panel for the same spectrum.
        px_panel = ijk_mrsi(1);
        py_panel = ijk_mrsi(2);
        fprintf('MRSI voxel: NIfTI [%d, %d, %d]   panel (%d,%d)\n', ...
            ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), px_panel, py_panel);

        % Bounds check
        if ijk_mrsi(1) < 1 || ijk_mrsi(1) > img_size(1) || ...
           ijk_mrsi(2) < 1 || ijk_mrsi(2) > img_size(2) || ...
           ijk_mrsi(3) < 1 || ijk_mrsi(3) > img_size(3)
            fprintf('Out of bounds\n');
            return;
        end

        % ---- World-coordinate report (mm, RAS) for T1 and MRSI ----
        % Purely informational. Wrapped ENTIRELY in try/catch (and indices
        % forced to columns) so a coordinate-reporting hiccup can NEVER abort
        % the FID/spectrum display further down.
        try
            w_bg = [];
            if isfield(hs,'bg') && isfield(hs.bg,'R') && ~isempty(hs.bg.R)
                ib = double(ijk_bg(:));                  % force column
                w_bg = hs.bg.R * [ib(1:3)-1; 1];
                fprintf('  T1   voxel [%3d %3d %3d] -> world [%+8.2f %+8.2f %+8.2f] mm (RAS)\n', ...
                    ib(1), ib(2), ib(3), w_bg(1), w_bg(2), w_bg(3));
            else
                fprintf('  T1 world unavailable (no hs.bg.R)\n');
            end
            S_mrsi = get_sform_3x4(p.nii.hdr);
            im = double(ijk_mrsi(:));                    % force column
            w_mrsi = S_mrsi * [im(1:3)-1; 1];
            vsz = [norm(S_mrsi(:,1)), norm(S_mrsi(:,2)), norm(S_mrsi(:,3))];
            fprintf('  MRSI voxel [%3d %3d %3d] -> world [%+8.2f %+8.2f %+8.2f] mm (RAS)\n', ...
                im(1), im(2), im(3), w_mrsi(1), w_mrsi(2), w_mrsi(3));
            fprintf('  MRSI voxel size: [%.2f %.2f %.2f] mm\n', vsz(1), vsz(2), vsz(3));
            if ~isempty(w_bg)
                d = norm(w_bg(1:3) - w_mrsi(1:3));
                if d <= max(vsz(1:2)), verdict = 'OK (<= one MRSI voxel)';
                else,                  verdict = 'MISALIGNED (> one MRSI voxel)'; end
                fprintf('  T1<->MRSI world mismatch: %.2f mm  %s\n', d, verdict);
            end
        catch ME
            fprintf('  [world-coord report skipped: %s]\n', ME.message);
        end

        % Get description
        descrip = '';
        if isfield(p.nii.hdr, 'descrip')
            descrip = strtrim(char(p.nii.hdr.descrip));
            fprintf('NIfTI description: "%s"\n', descrip);
        end
        
        % DETECT DATA TYPE
        nVol = img_size(4);
        is_metabolite = false;
        met_name = 'Unknown';
        data_type = '';  % '4D_metabolite', '3D_CRLB', '3D_LW', '3D_SNR', 'signal'
        
        % Case 1: 4D Metabolite [Conc, CRLB, LW, SNR]
        if nVol == 4
            if contains(descrip, 'CRLB', 'IgnoreCase', true) || ...
               contains(descrip, 'LW', 'IgnoreCase', true) || ...
               contains(descrip, 'SNR', 'IgnoreCase', true) || ...
               contains(descrip, '+', 'IgnoreCase', true)
                
                fprintf('✓ 4D metabolite format detected\n');
                is_metabolite = true;
                data_type = '4D_metabolite';
                
                % Extract metabolite name
                met_name = extract_metabolite_name(descrip);
            else
                fprintf('→ 4D signal data (will show spectrum)\n');
                is_metabolite = false;
                data_type = 'signal';
            end
            
        % Case 2: 3D Metabolite maps (NEW) or single volume
        elseif nVol == 1
            % Check for 3D CRLB map
            if contains(descrip, 'CRLB', 'IgnoreCase', true) || ...
               contains(descrip, '_CRLB', 'IgnoreCase', true)
                fprintf('✓ 3D CRLB map detected\n');
                is_metabolite = true;
                data_type = '3D_CRLB';
                met_name = extract_metabolite_name(descrip);
                
            % Check for 3D Linewidth map
            elseif contains(descrip, 'Linewidth', 'IgnoreCase', true)
                fprintf('✓ 3D Linewidth map detected\n');
                is_metabolite = true;
                data_type = '3D_LW';
                met_name = 'Linewidth';
                
            % Check for 3D SNR map
            elseif contains(descrip, 'SNR', 'IgnoreCase', true)
                fprintf('✓ 3D SNR map detected\n');
                is_metabolite = true;
                data_type = '3D_SNR';
                met_name = 'SNR';
                
            else
                fprintf('→ Single volume data (signal or other)\n');
                is_metabolite = false;
                data_type = 'signal';
            end
        else
            fprintf('→ Multi-volume signal data\n');
            is_metabolite = false;
            data_type = 'signal';
        end
        
        % DISPLAY
        if is_metabolite
            fprintf('Displaying metabolite info (type: %s)...\n', data_type);
            show_metabolite_info(hs, p, ijk_mrsi, met_name, data_type, img_size(3));
        else
            fprintf('Displaying spectrum...\n');
            show_spectrum(hs, p, ijk_mrsi, img_size);
        end
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        
        if ishandle(hs.ax(4))
            cla(hs.ax(4));
            axis(hs.ax(4), 'on');
            set(hs.ax(4), 'Color', 'k');
            text(0.5, 0.5, {sprintf('Error: %s', ME.message)}, ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', 10, 'Color', 'r', 'Parent', hs.ax(4));
        end
    end
end

%% Extract metabolite name from description
function met_name = extract_metabolite_name(descrip)
    met_name = 'Unknown';
    
    % Common metabolites
    metabolites = {'CrPCr', 'NAANAAG', 'GPCPCh', 'GluGln', 'Ins', ...
                   'NAA', 'Cho', 'Cr', 'Glu', 'Gln', 'mI'};
    
    for i = 1:length(metabolites)
        if contains(descrip, metabolites{i}, 'IgnoreCase', true)
            met_name = metabolites{i};
            return;
        end
    end
    
    % Try to extract from description format "MetName+CRLB+LW+SNR" or "MetName_CRLB"
    tokens = split(descrip, {'+', '[', ' ', '_'});
    if ~isempty(tokens) && ~isempty(tokens{1})
        met_name = strtrim(tokens{1});
    end
end

%% Show metabolite info - handles both 4D and 3D formats
function show_metabolite_info(hs, p, ijk_mrsi, met_name, data_type, nz)
    % Extract values based on data type
    
    if strcmp(data_type, '4D_metabolite')
        % 4D format: [Conc, CRLB, LW, SNR]
        conc = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), 1);
        crlb_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), 2);
        lw_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), 3);
        snr_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), 4);
        
        fprintf('4D Metabolite - Voxel [%d, %d, %d]:\n', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        fprintf('  %s = %.3f\n', met_name, conc);
        fprintf('  CRLB = %.1f%%\n', crlb_val);
        fprintf('  LW = %.3f ppm\n', lw_val);
        fprintf('  SNR = %.1f\n', snr_val);
        
        display_metabolite_4d(hs, ijk_mrsi, met_name, conc, crlb_val, lw_val, snr_val, nz);
        
    elseif strcmp(data_type, '3D_CRLB')
        % 3D CRLB map only - handle 2D/3D indexing
        img_dims = size(p.nii.img);
        if length(img_dims) == 2
            crlb_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2));
        else
            crlb_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        end
        
        fprintf('3D CRLB - Voxel [%d, %d, %d]:\n', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        fprintf('  %s CRLB = %.1f%%\n', met_name, crlb_val);
        
        display_single_quality_metric(hs, ijk_mrsi, met_name, 'CRLB', crlb_val, '%', nz);
        
    elseif strcmp(data_type, '3D_LW')
        % 3D Linewidth map only - handle 2D/3D indexing
        img_dims = size(p.nii.img);
        if length(img_dims) == 2
            lw_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2));
        else
            lw_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        end
        
        fprintf('3D Linewidth - Voxel [%d, %d, %d]:\n', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        fprintf('  LW = %.3f ppm\n', lw_val);
        
        display_single_quality_metric(hs, ijk_mrsi, 'All Metabolites', 'Linewidth', lw_val, 'ppm', nz);
        
    elseif strcmp(data_type, '3D_SNR')
        % 3D SNR map only - handle 2D/3D indexing
        img_dims = size(p.nii.img);
        if length(img_dims) == 2
            snr_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2));
        else
            snr_val = p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        end
        
        fprintf('3D SNR - Voxel [%d, %d, %d]:\n', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3));
        fprintf('  SNR = %.1f\n', snr_val);
        
        display_single_quality_metric(hs, ijk_mrsi, 'All Metabolites', 'SNR', snr_val, '', nz);
    end
end

%% Display 4D metabolite (original)
function display_metabolite_4d(hs, ijk_mrsi, met_name, conc, crlb_val, lw_val, snr_val, nz)
    if ~ishandle(hs.ax(4)), return; end
    
    % Store/restore position
    if ~isappdata(hs.ax(4), 'OriginalPosition')
        setappdata(hs.ax(4), 'OriginalPosition', get(hs.ax(4), 'Position'));
    end
    
    % Clear old content
    clear_axis_content(hs.ax(4));
    
    % Restore settings
    set(hs.ax(4), 'Position', getappdata(hs.ax(4), 'OriginalPosition'));
    axis(hs.ax(4), 'off');
    set(hs.ax(4), 'Color', 'k');
    
    % Get quality colors and status
    [crlb_color, crlb_status] = get_crlb_quality(crlb_val);
    [lw_color, lw_status] = get_lw_quality(lw_val);
    [snr_color, snr_status] = get_snr_quality(snr_val);
    
    % Display
    y_pos = 0.95;
    line_height = 0.055;
    
    % Voxel coordinates (Cyan)
    if nz == 1
        text(0.1, y_pos, sprintf('Voxel [%d, %d]', ijk_mrsi(1), ijk_mrsi(2)), ...
            'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', ...
            'Color', [0 1 1], 'Parent', hs.ax(4));
    else
        text(0.1, y_pos, sprintf('Voxel [%d, %d, %d]', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3)), ...
            'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', ...
            'Color', [0 1 1], 'Parent', hs.ax(4));
    end
    y_pos = y_pos - line_height * 1.5;
    
    % Metabolite name (Yellow)
    text(0.1, y_pos, [met_name ':'], ...
        'Units', 'normalized', 'FontSize', 12, ...
        'Color', [1 1 0], 'Parent', hs.ax(4));
    y_pos = y_pos - line_height;
    
    % Concentration (White)
    text(0.1, y_pos, sprintf('%.3f', conc), ...
        'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
        'Color', [1 1 1], 'Parent', hs.ax(4));
    y_pos = y_pos - line_height * 1.5;
    
    % Quality header (Yellow)
    text(0.1, y_pos, 'Quality Metrics:', ...
        'Units', 'normalized', 'FontSize', 11, ...
        'Color', [1 1 0], 'Parent', hs.ax(4));
    y_pos = y_pos - line_height * 1.2;
    
    % CRLB
    text(0.1, y_pos, sprintf('CRLB: %.1f%% (%s)', crlb_val, crlb_status), ...
        'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', crlb_color, 'Parent', hs.ax(4));
    y_pos = y_pos - line_height;
    
    % Linewidth
    text(0.1, y_pos, sprintf('LW: %.3f ppm (%s)', lw_val, lw_status), ...
        'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', lw_color, 'Parent', hs.ax(4));
    y_pos = y_pos - line_height;
    
    % SNR
    text(0.1, y_pos, sprintf('SNR: %.1f (%s)', snr_val, snr_status), ...
        'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', snr_color, 'Parent', hs.ax(4));
    
    fprintf('✓ 4D Metabolite info displayed\n');
end

%% Display single quality metric (NEW - for 3D maps)
function display_single_quality_metric(hs, ijk_mrsi, met_name, metric_name, value, unit, nz)
    if ~ishandle(hs.ax(4)), return; end
    
    % Clear old content
    clear_axis_content(hs.ax(4));
    
    % Settings
    if ~isappdata(hs.ax(4), 'OriginalPosition')
        setappdata(hs.ax(4), 'OriginalPosition', get(hs.ax(4), 'Position'));
    end
    set(hs.ax(4), 'Position', getappdata(hs.ax(4), 'OriginalPosition'));
    axis(hs.ax(4), 'off');
    set(hs.ax(4), 'Color', 'k');
    
    % Get quality color and status
    if strcmp(metric_name, 'CRLB')
        [qual_color, qual_status] = get_crlb_quality(value);
        value_str = sprintf('%.1f%s', value, unit);
    elseif strcmp(metric_name, 'Linewidth')
        [qual_color, qual_status] = get_lw_quality(value);
        value_str = sprintf('%.3f %s', value, unit);
    elseif strcmp(metric_name, 'SNR')
        [qual_color, qual_status] = get_snr_quality(value);
        value_str = sprintf('%.1f%s', value, unit);
    else
        qual_color = [1 1 1];
        qual_status = '';
        value_str = sprintf('%.3f %s', value, unit);
    end
    
    % Display
    y_pos = 0.90;
    line_height = 0.08;
    
    % Voxel coordinates (Cyan)
    if nz == 1
        text(0.1, y_pos, sprintf('Voxel [%d, %d]', ijk_mrsi(1), ijk_mrsi(2)), ...
            'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', ...
            'Color', [0 1 1], 'Parent', hs.ax(4));
    else
        text(0.1, y_pos, sprintf('Voxel [%d, %d, %d]', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3)), ...
            'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', ...
            'Color', [0 1 1], 'Parent', hs.ax(4));
    end
    y_pos = y_pos - line_height * 1.5;
    
    % Metabolite name (Yellow)
    text(0.1, y_pos, [met_name ':'], ...
        'Units', 'normalized', 'FontSize', 13, ...
        'Color', [1 1 0], 'Parent', hs.ax(4));
    y_pos = y_pos - line_height * 1.5;
    
    % Metric name (White)
    text(0.1, y_pos, [metric_name ':'], ...
        'Units', 'normalized', 'FontSize', 13, ...
        'Color', [1 1 1], 'Parent', hs.ax(4));
    y_pos = y_pos - line_height * 1.2;
    
    % Value with quality color (LARGE)
    text(0.1, y_pos, value_str, ...
        'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', ...
        'Color', qual_color, 'Parent', hs.ax(4));
    y_pos = y_pos - line_height * 1.2;
    
    % Status
    if ~isempty(qual_status)
        text(0.1, y_pos, sprintf('(%s)', qual_status), ...
            'Units', 'normalized', 'FontSize', 12, ...
            'Color', qual_color, 'Parent', hs.ax(4));
        y_pos = y_pos - line_height * 1.5;
    end
    
    % Add threshold info
    y_pos = y_pos - line_height * 0.5;
    text(0.1, y_pos, 'Quality Thresholds:', ...
        'Units', 'normalized', 'FontSize', 9, ...
        'Color', [0.7 0.7 0.7], 'Parent', hs.ax(4));
    y_pos = y_pos - line_height * 0.7;
    
    if strcmp(metric_name, 'CRLB')
        text(0.1, y_pos, '<10%=Excellent, 10-20%=Good, 20-50%=Acceptable', ...
            'Units', 'normalized', 'FontSize', 8, ...
            'Color', [0.7 0.7 0.7], 'Parent', hs.ax(4));
    elseif strcmp(metric_name, 'Linewidth')
        text(0.1, y_pos, '<0.05=Excellent, 0.05-0.1=Good, 0.1-0.15=Acceptable', ...
            'Units', 'normalized', 'FontSize', 8, ...
            'Color', [0.7 0.7 0.7], 'Parent', hs.ax(4));
    elseif strcmp(metric_name, 'SNR')
        text(0.1, y_pos, '>30=Excellent, 10-30=Good, 3-10=Acceptable', ...
            'Units', 'normalized', 'FontSize', 8, ...
            'Color', [0.7 0.7 0.7], 'Parent', hs.ax(4));
    end
    
    fprintf('✓ 3D quality metric displayed\n');
end

%% Helper: Clear axis content
function clear_axis_content(ax)
    cla(ax);
    if isappdata(ax, 'TimeAxes')
        old_ax = getappdata(ax, 'TimeAxes');
        if ishandle(old_ax), delete(old_ax); end
        rmappdata(ax, 'TimeAxes');
    end
    if isappdata(ax, 'FreqAxes')
        old_ax = getappdata(ax, 'FreqAxes');
        if ishandle(old_ax), delete(old_ax); end
        rmappdata(ax, 'FreqAxes');
    end
end

%% Helper: Get CRLB quality
function [color, status] = get_crlb_quality(val)
    if val < 10
        color = [0.2 1.0 0.2]; status = 'Excellent';
    elseif val < 20
        color = [0.5 1.0 0.3]; status = 'Good';
    elseif val < 50
        color = [1.0 0.8 0.0]; status = 'Acceptable';
    else
        color = [1.0 0.2 0.0]; status = 'Poor';
    end
end

%% Helper: Get LW quality
function [color, status] = get_lw_quality(val)
    if val < 0.05
        color = [0.2 1.0 0.2]; status = 'Excellent';
    elseif val < 0.1
        color = [0.5 1.0 0.3]; status = 'Good';
    elseif val < 0.15
        color = [1.0 0.8 0.0]; status = 'Acceptable';
    else
        color = [1.0 0.2 0.0]; status = 'Poor';
    end
end

%% Helper: Get SNR quality
function [color, status] = get_snr_quality(val)
    if val > 30
        color = [0.2 1.0 0.2]; status = 'Excellent';
    elseif val >= 10
        color = [0.5 1.0 0.3]; status = 'Good';
    elseif val >= 3
        color = [1.0 0.8 0.0]; status = 'Acceptable';
    else
        color = [1.0 0.2 0.0]; status = 'Poor';
    end
end

%% Show spectrum (original - unchanged)
function show_spectrum(hs, p, ijk_mrsi, img_size)
    show_freq = isfield(p, 'ftSpec_smooth_w') && ~isempty(p.ftSpec_smooth_w);
    
    fid = squeeze(p.nii.img(ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3), :));
    fid = fid(:);
    
    if all(fid == 0) || all(isnan(fid))
        if ishandle(hs.ax(4))
            cla(hs.ax(4));
            axis(hs.ax(4), 'off');
            set(hs.ax(4), 'Color', 'k');
            text(0.5, 0.5, 'No signal', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 12, ...
                'Color', 'y', 'Parent', hs.ax(4));
        end
        return;
    end
    
    if isfield(p, 'scl_slope') && p.scl_slope ~= 0
        fid = single(fid) * p.scl_slope + p.scl_inter;
    end
    
    dwelltime = 5e-6;
    if isfield(p.nii.hdr, 'pixdim') && p.nii.hdr.pixdim(5) > 0
        dwelltime = double(p.nii.hdr.pixdim(5));
    end
    
    t_ms = (0:length(fid)-1)' * dwelltime * 1000;
    
    if show_freq
        try
            % Extract the spectrum EXACTLY like mrsi_integration_panel:
            % raw real(ftSpec.data) at the voxel, with ftSpec.ppm.  We do NOT
            % use op_CSItoMRS here because it conditionally ifft's the data
            % when the 'spectralft' flag is not set, producing a different
            % curve.  Indexing the data directly guarantees the nii_viewer
            % spectrum is identical to the integration panel.
            %
            % The ccav_w overlay (written by mrsi_on_t1_map) is correct, but
            % ftSpec_smooth is stored 180 deg rotated (flipped in BOTH x and y)
            % relative to ccav_w.  The overlay voxel shown at NIfTI (a,b) is
            % ccav_w(Nx-a+1, b); the ftSpec voxel that PHYSICALLY matches it is
            % (ix = a, iy = Ny - b + 1).  Pull that so the spectrum belongs to
            % the clicked overlay voxel.
            ft   = p.ftSpec_smooth_w;
            fDim = ft.dims.f;  if fDim == 0, fDim = ft.dims.t; end
            Ny_d = ft.sz(ft.dims.y);
            ix_d = ijk_mrsi(1);
            iy_d = Ny_d - ijk_mrsi(2) + 1;
            specCube = permute(ft.data, [fDim, ft.dims.y, ft.dims.x]);  % (Nf, Ny, Nx)
            vox_MRS  = struct('specs', squeeze(specCube(:, iy_d, ix_d)), ...
                              'ppm',   ft.ppm(:));
        catch
            show_freq = false;
        end
    end
    
    if ~ishandle(hs.ax(4)), return; end
    
    if ~isappdata(hs.ax(4), 'OriginalPosition')
        setappdata(hs.ax(4), 'OriginalPosition', get(hs.ax(4), 'Position'));
    end
    origPos = getappdata(hs.ax(4), 'OriginalPosition');
    
    clear_axis_content(hs.ax(4));
    
    if show_freq
        % Panel index of the pulled ftSpec voxel (ix=a, iy=Ny-b+1) is (a, b),
        % which matches the integration panel's index for the same spectrum.
        px_p = ijk_mrsi(1);
        py_p = ijk_mrsi(2);

        width = (origPos(3) - 0.08) / 2;
        height = origPos(4) - 0.15;
        bottom = origPos(2) + 0.08;

        ax_t = axes('Parent', get(hs.ax(4), 'Parent'), 'Units', 'normalized');
        set(ax_t, 'Position', [origPos(1)+0.03 bottom width height]);
        plot(ax_t, t_ms, real(fid), 'b-', 'LineWidth', 2);
        set(ax_t, 'XColor', 'w', 'YColor', 'w', 'Color', 'k', 'FontSize', 12, 'Box', 'off');
        xlabel(ax_t, 'Time (ms)', 'FontSize', 14, 'Color', 'w');
        title(ax_t, sprintf('Time (Real)  panel (%d,%d)', px_p, py_p), ...
            'FontSize', 12, 'Color', 'w');
        xlim(ax_t, [0 max(t_ms)]);

        ax_f = axes('Parent', get(hs.ax(4), 'Parent'), 'Units', 'normalized');
        set(ax_f, 'Position', [origPos(1)+0.06+width bottom width height]);

        if isfield(vox_MRS, 'ppm') && isfield(vox_MRS, 'specs')
            plot(ax_f, vox_MRS.ppm, real(vox_MRS.specs), 'b-', 'LineWidth', 2);
            set(ax_f, 'XDir', 'reverse', 'XColor', 'w', 'YColor', 'w', ...
                'Color', 'k', 'FontSize', 12, 'Box', 'off');
            xlabel(ax_f, 'ppm', 'FontSize', 14, 'Color', 'w');
            title(ax_f, sprintf('Spectrum (Real)  panel (%d,%d)', px_p, py_p), ...
                'FontSize', 12, 'Color', 'w');
            xlim(ax_f, [0 6]);
        end
        
        setappdata(hs.ax(4), 'TimeAxes', ax_t);
        setappdata(hs.ax(4), 'FreqAxes', ax_f);
    else
        cla(hs.ax(4));
        set(hs.ax(4), 'Position', [origPos(1)+0.03 origPos(2)+0.08 origPos(3)-0.06 origPos(4)-0.15]);
        plot(hs.ax(4), t_ms, real(fid), 'b-', 'LineWidth', 2);
        set(hs.ax(4), 'XColor', 'w', 'YColor', 'w', 'Color', 'k', 'FontSize', 12, 'Box', 'off');
        xlabel(hs.ax(4), 'Time (ms)', 'FontSize', 14, 'Color', 'w');
        title(hs.ax(4), sprintf('Time (Real) [%d,%d,%d]', ijk_mrsi(1), ijk_mrsi(2), ijk_mrsi(3)), ...
            'FontSize', 12, 'Color', 'w');
        xlim(hs.ax(4), [0 max(t_ms)]);
    end
    
    fprintf('✓ Spectrum displayed\n');
end