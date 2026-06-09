function nii_viewer_with_spectra(t1_filename, mrsi_filename, ftSpec_smooth)
% NII viewer style interface with click-to-show-spectrum functionality
% Mimics nii_viewer but shows spectra instead of time courses
%
% Usage: nii_viewer_with_spectra('T1.nii.gz', 'mrsi.nii.gz', ftSpec_smooth)

    % Load data
    t1_nii = nii_tool('load', t1_filename);
    mrsi_nii = nii_tool('load', mrsi_filename);
    
    % Create main figure (similar to nii_viewer layout)
    fig = figure('Name', sprintf('MRSI Viewer - %s', mrsi_filename), ...
                 'Position', [100 100 1200 800], ...
                 'Color', 'k', ...
                 'KeyPressFcn', @key_callback);
    
    % Create panels similar to nii_viewer
    % Main image panel
    img_panel = uipanel(fig, 'Position', [0.02 0.15 0.65 0.83], ...
                        'BackgroundColor', 'k', 'BorderType', 'none');
    
    % Control panel (top)
    ctrl_panel = uipanel(fig, 'Position', [0.02 0.02 0.65 0.12], ...
                         'BackgroundColor', [0.2 0.2 0.2]);
    
    % Spectrum panel (right)
    spec_panel = uipanel(fig, 'Position', [0.68 0.02 0.30 0.96], ...
                         'BackgroundColor', 'k', 'BorderType', 'none');
    
    % Get image data
    t1_img = double(t1_nii.img);
    mrsi_img = mrsi_nii.img;
    
    % Normalize T1 for display
    t1_img = (t1_img - min(t1_img(:))) / (max(t1_img(:)) - min(t1_img(:)));
    
    % Initialize current slice indices
    current_slices = [round(size(t1_img,1)/2), round(size(t1_img,2)/2), round(size(t1_img,3)/2)];
    
    % Create axes for the three views
    ax_cor = axes('Parent', img_panel, 'Position', [0.02 0.52 0.46 0.46]);
    ax_sag = axes('Parent', img_panel, 'Position', [0.52 0.52 0.46 0.46]);
    ax_ax = axes('Parent', img_panel, 'Position', [0.02 0.02 0.46 0.46]);
    
    % Create spectrum axes
    ax_spec = axes('Parent', spec_panel, 'Position', [0.1 0.55 0.85 0.4]);
    
    % Info text area in spectrum panel
    info_text = annotation(spec_panel, 'textbox', [0.05 0.05 0.9 0.45], ...
                          'String', 'Click on MRSI voxel to show spectrum', ...
                          'Color', 'w', 'FontSize', 10, 'EdgeColor', 'none', ...
                          'BackgroundColor', 'k', 'FitBoxToText', 'off');
    
    % Create control elements
    create_controls(ctrl_panel, current_slices, size(t1_img));
    
    % Initial display
    update_display();
    
    % Store data for callbacks
    setappdata(fig, 't1_img', t1_img);
    setappdata(fig, 'mrsi_img', mrsi_img);
    setappdata(fig, 'ftSpec_smooth', ftSpec_smooth);
    setappdata(fig, 'current_slices', current_slices);
    setappdata(fig, 'axes', struct('cor', ax_cor, 'sag', ax_sag, 'ax', ax_ax, 'spec', ax_spec));
    setappdata(fig, 'info_text', info_text);
    
    function update_display()
        % Update all three orthogonal views
        
        % Coronal view (Y slice)
        y_slice = current_slices(2);
        cor_img = squeeze(t1_img(:, y_slice, :))';
        cor_mrsi = squeeze(abs(mrsi_img(:, y_slice, :, :)));
        
        imagesc(ax_cor, cor_img);
        colormap(ax_cor, 'gray');
        axis(ax_cor, 'equal', 'tight', 'xy');
        title(ax_cor, sprintf('Coronal Y=%d', y_slice), 'Color', 'w');
        hold(ax_cor, 'on');
        overlay_mrsi(ax_cor, cor_mrsi, size(cor_img), 'coronal');
        hold(ax_cor, 'off');
        
        % Sagittal view (X slice)
        x_slice = current_slices(1);
        sag_img = squeeze(t1_img(x_slice, :, :))';
        sag_mrsi = squeeze(abs(mrsi_img(x_slice, :, :, :)));
        
        imagesc(ax_sag, sag_img);
        colormap(ax_sag, 'gray');
        axis(ax_sag, 'equal', 'tight', 'xy');
        title(ax_sag, sprintf('Sagittal X=%d', x_slice), 'Color', 'w');
        hold(ax_sag, 'on');
        overlay_mrsi(ax_sag, sag_mrsi, size(sag_img), 'sagittal');
        hold(ax_sag, 'off');
        
        % Axial view (Z slice)
        z_slice = current_slices(3);
        ax_img = t1_img(:, :, z_slice);
        ax_mrsi = abs(mrsi_img(:, :, z_slice, :));
        
        imagesc(ax_ax, ax_img);
        colormap(ax_ax, 'gray');
        axis(ax_ax, 'equal', 'tight', 'xy');
        title(ax_ax, sprintf('Axial Z=%d', z_slice), 'Color', 'w');
        hold(ax_ax, 'on');
        overlay_mrsi(ax_ax, ax_mrsi, size(ax_img), 'axial');
        hold(ax_ax, 'off');
        
        % Set click callbacks
        set(ax_cor, 'ButtonDownFcn', @(src,evt) click_callback(src, evt, 'coronal'));
        set(ax_sag, 'ButtonDownFcn', @(src,evt) click_callback(src, evt, 'sagittal'));
        set(ax_ax, 'ButtonDownFcn', @(src,evt) click_callback(src, evt, 'axial'));
        
        % Update stored current slices
        setappdata(fig, 'current_slices', current_slices);
    end
    
    function overlay_mrsi(ax, mrsi_data, img_size, view_type)
        % Overlay MRSI voxels as colored rectangles
        
        if isempty(mrsi_data) || all(mrsi_data(:) == 0)
            return;
        end
        
        % Get MRSI dimensions
        mrsi_size = size(mrsi_data);
        if length(mrsi_size) < 2
            return;
        end
        
        % Calculate scaling
        x_scale = img_size(2) / mrsi_size(2);
        y_scale = img_size(1) / mrsi_size(1);
        
        % Draw MRSI voxels
        for i = 1:mrsi_size(1)
            for j = 1:mrsi_size(2)
                if length(mrsi_size) > 2 && any(mrsi_data(i, j, :) > 0)
                    x_pos = (j - 0.5) * x_scale;
                    y_pos = (mrsi_size(1) - i + 0.5) * y_scale; % Flip Y
                    
                    rectangle(ax, 'Position', [x_pos-x_scale/2, y_pos-y_scale/2, x_scale, y_scale], ...
                             'FaceColor', [1 0 0 0.3], 'EdgeColor', 'r', 'LineWidth', 1);
                elseif length(mrsi_size) == 2 && mrsi_data(i, j) > 0
                    x_pos = (j - 0.5) * x_scale;
                    y_pos = (mrsi_size(1) - i + 0.5) * y_scale; % Flip Y
                    
                    rectangle(ax, 'Position', [x_pos-x_scale/2, y_pos-y_scale/2, x_scale, y_scale], ...
                             'FaceColor', [1 0 0 0.3], 'EdgeColor', 'r', 'LineWidth', 1);
                end
            end
        end
    end
    
    function click_callback(src, evt, view_type)
        % Handle clicks on MRSI voxels
        
        pos = get(src, 'CurrentPoint');
        x_click = pos(1, 1);
        y_click = pos(1, 2);
        
        % Convert to MRSI coordinates based on current view and slice
        [vox_x, vox_y, vox_z] = convert_click_to_voxel(x_click, y_click, view_type);
        
        if ~isempty(vox_x)
            display_spectrum(vox_x, vox_y, vox_z);
        end
    end
    
    function [vox_x, vox_y, vox_z] = convert_click_to_voxel(x_click, y_click, view_type)
        % Convert click coordinates to MRSI voxel indices
        
        ftSpec = getappdata(fig, 'ftSpec_smooth');
        [~, n_x, n_y] = size(ftSpec.data);
        current_slices = getappdata(fig, 'current_slices');
        t1_img = getappdata(fig, 't1_img');
        
        vox_x = []; vox_y = []; vox_z = [];
        
        switch view_type
            case 'coronal' % YZ plane, X fixed
                img_size = [size(t1_img, 3), size(t1_img, 1)]; % [Z, X]
                x_scale = img_size(2) / n_x;
                z_scale = img_size(1) / 1; % Assuming single Z slice
                
                vox_x = max(1, min(n_x, round(y_click / x_scale)));
                vox_y = current_slices(2); % Fixed Y slice
                vox_z = 1;
                
            case 'sagittal' % XZ plane, Y fixed  
                img_size = [size(t1_img, 3), size(t1_img, 2)]; % [Z, Y]
                y_scale = img_size(2) / n_y;
                z_scale = img_size(1) / 1;
                
                vox_x = current_slices(1); % Fixed X slice
                vox_y = max(1, min(n_y, round(y_click / y_scale)));
                vox_z = 1;
                
            case 'axial' % XY plane, Z fixed
                img_size = [size(t1_img, 2), size(t1_img, 1)]; % [Y, X]
                x_scale = img_size(2) / n_x;
                y_scale = img_size(1) / n_y;
                
                vox_x = max(1, min(n_x, round(y_click / x_scale)));
                vox_y = max(1, min(n_y, round(x_click / y_scale)));
                vox_z = 1;
        end
        
        % Validate voxel indices
        if vox_x < 1 || vox_x > n_x || vox_y < 1 || vox_y > n_y
            vox_x = []; vox_y = []; vox_z = [];
        end
    end
    
    function display_spectrum(vox_x, vox_y, vox_z)
        % Display spectrum for selected voxel
        
        ftSpec = getappdata(fig, 'ftSpec_smooth');
        axes_struct = getappdata(fig, 'axes');
        info_text = getappdata(fig, 'info_text');
        
        % Get spectrum
        spectrum = ftSpec.data(:, vox_x, vox_y);
        ppm = ftSpec.ppm;
        
        % Clear and plot spectrum
        cla(axes_struct.spec);
        
        plot(axes_struct.spec, ppm, abs(spectrum), 'b-', 'LineWidth', 2);
        hold(axes_struct.spec, 'on');
        plot(axes_struct.spec, ppm, real(spectrum), 'r--', 'LineWidth', 1);
        hold(axes_struct.spec, 'off');
        
        % Format spectrum plot
        set(axes_struct.spec, 'XDir', 'reverse', 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
        xlabel(axes_struct.spec, 'Chemical Shift (ppm)', 'Color', 'w');
        ylabel(axes_struct.spec, 'Signal', 'Color', 'w');
        title(axes_struct.spec, sprintf('Voxel (%d,%d)', vox_x, vox_y), 'Color', 'w');
        legend(axes_struct.spec, {'Magnitude', 'Real'}, 'TextColor', 'w', 'Location', 'best');
        grid(axes_struct.spec, 'on');
        xlim(axes_struct.spec, [0 4.5]);
        
        % Add metabolite annotations
        add_metabolite_lines(axes_struct.spec);
        
        % Update info text
        snr = max(abs(spectrum)) / std(abs(spectrum(1:50))); % Rough SNR estimate
        info_str = sprintf(['Voxel: (%d, %d)\n' ...
                           'Peak intensity: %.2e\n' ...
                           'Estimated SNR: %.1f\n\n' ...
                           'Metabolite regions:\n' ...
                           'NAA: ~2.0 ppm\n' ...
                           'Cr: ~3.0 ppm\n' ...
                           'Cho: ~3.2 ppm\n' ...
                           'Lac: ~1.3 ppm'], ...
                          vox_x, vox_y, max(abs(spectrum)), snr);
        
        set(info_text, 'String', info_str);
        
        fprintf('Displayed spectrum for voxel (%d, %d)\n', vox_x, vox_y);
    end
    
    function add_metabolite_lines(ax)
        % Add metabolite reference lines
        metabolites = struct('NAA', 2.0, 'Cr', 3.0, 'Cho', 3.2, 'Lac', 1.3);
        y_lim = get(ax, 'YLim');
        
        hold(ax, 'on');
        met_names = fieldnames(metabolites);
        for i = 1:length(met_names)
            ppm_val = metabolites.(met_names{i});
            line(ax, [ppm_val ppm_val], y_lim, 'Color', [0.7 0.7 0.7], ...
                 'LineStyle', ':', 'LineWidth', 1);
            text(ax, ppm_val, y_lim(2)*0.9, met_names{i}, 'Color', 'w', ...
                 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
        hold(ax, 'off');
    end
    
    function create_controls(parent, slices, img_size)
        % Create slice navigation controls
        
        % X slice control
        uicontrol(parent, 'Style', 'text', 'Position', [10 60 30 20], ...
                  'String', 'X:', 'ForegroundColor', 'w', 'BackgroundColor', [0.2 0.2 0.2]);
        x_slider = uicontrol(parent, 'Style', 'slider', 'Position', [45 60 150 20], ...
                            'Min', 1, 'Max', img_size(1), 'Value', slices(1), ...
                            'SliderStep', [1/(img_size(1)-1) 10/(img_size(1)-1)], ...
                            'Callback', @(src,evt) slider_callback(1, src));
        
        % Y slice control  
        uicontrol(parent, 'Style', 'text', 'Position', [10 35 30 20], ...
                  'String', 'Y:', 'ForegroundColor', 'w', 'BackgroundColor', [0.2 0.2 0.2]);
        y_slider = uicontrol(parent, 'Style', 'slider', 'Position', [45 35 150 20], ...
                            'Min', 1, 'Max', img_size(2), 'Value', slices(2), ...
                            'SliderStep', [1/(img_size(2)-1) 10/(img_size(2)-1)], ...
                            'Callback', @(src,evt) slider_callback(2, src));
        
        % Z slice control
        uicontrol(parent, 'Style', 'text', 'Position', [10 10 30 20], ...
                  'String', 'Z:', 'ForegroundColor', 'w', 'BackgroundColor', [0.2 0.2 0.2]);
        z_slider = uicontrol(parent, 'Style', 'slider', 'Position', [45 10 150 20], ...
                            'Min', 1, 'Max', img_size(3), 'Value', slices(3), ...
                            'SliderStep', [1/(img_size(3)-1) 10/(img_size(3)-1)], ...
                            'Callback', @(src,evt) slider_callback(3, src));
        
        % Instructions
        uicontrol(parent, 'Style', 'text', 'Position', [220 10 300 60], ...
                  'String', 'Use sliders to navigate slices. Click red MRSI voxels to view spectra.', ...
                  'ForegroundColor', 'w', 'BackgroundColor', [0.2 0.2 0.2], ...
                  'HorizontalAlignment', 'left');
    end
    
    function slider_callback(dim, src)
        % Handle slice slider changes
        current_slices = getappdata(fig, 'current_slices');
        current_slices(dim) = round(get(src, 'Value'));
        setappdata(fig, 'current_slices', current_slices);
        update_display();
    end
    
    function key_callback(src, evt)
        % Handle keyboard shortcuts
        current_slices = getappdata(fig, 'current_slices');
        t1_img = getappdata(fig, 't1_img');
        
        switch evt.Key
            case 'uparrow'
                current_slices(3) = min(size(t1_img,3), current_slices(3) + 1);
            case 'downarrow'
                current_slices(3) = max(1, current_slices(3) - 1);
            case 'leftarrow'
                current_slices(1) = max(1, current_slices(1) - 1);
            case 'rightarrow'
                current_slices(1) = min(size(t1_img,1), current_slices(1) + 1);
            case 'pageup'
                current_slices(2) = min(size(t1_img,2), current_slices(2) + 1);
            case 'pagedown'
                current_slices(2) = max(1, current_slices(2) - 1);
            otherwise
                return;
        end
        
        setappdata(fig, 'current_slices', current_slices);
        update_display();
    end
end