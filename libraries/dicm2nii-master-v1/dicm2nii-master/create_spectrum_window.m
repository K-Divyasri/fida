function fh = create_spectrum_window(spectrum_data, x_axis, ijk, p)
    % Create spectrum display window
    [~, fname] = fileparts(p.nii.hdr.file_name);
    if length(fname) > 40
        fname = [fname(1:37) '...'];
    end
    
    fh = figure('Name', sprintf('MRSI Spectrum - %s - Voxel [%d,%d,%d]', ...
        fname, ijk(1), ijk(2), ijk(3)), ...
        'NumberTitle', 'off', 'Position', [100 100 800 400], ...
        'MenuBar', 'none', 'Toolbar', 'figure');
    
    if ~isreal(spectrum_data)
        % Complex data: show real and imaginary parts
        subplot(2,1,1);
        plot(x_axis, real(spectrum_data), 'b-', 'LineWidth', 1.5);
        title('Real Part'); grid on;
        xlabel(get_axis_label(p)); ylabel('Signal Intensity');
        
        subplot(2,1,2);
        plot(x_axis, imag(spectrum_data), 'r-', 'LineWidth', 1.5);
        title('Imaginary Part'); grid on;
        xlabel(get_axis_label(p)); ylabel('Signal Intensity');
    else
        % Real data: show magnitude
        plot(x_axis, spectrum_data, 'b-', 'LineWidth', 1.5);
        title(sprintf('MRSI Spectrum - Voxel [%d,%d,%d]', ijk(1), ijk(2), ijk(3)));
        grid on; xlabel(get_axis_label(p)); ylabel('Signal Intensity');
        
        % Reverse x-axis for ppm
        if contains(get_axis_label(p), 'ppm') || (max(x_axis) < 15 && min(x_axis) >= 0)
            set(gca, 'XDir', 'reverse');
        end
    end
    
    % Add context menu
    cmenu = uicontextmenu('Parent', fh);
    uimenu(cmenu, 'Label', 'Export Data to Workspace', ...
        'Callback', @(~,~) export_spectrum_data(spectrum_data, x_axis, ijk));
    uimenu(cmenu, 'Label', 'Copy Figure', 'Callback', @(~,~) print('-dbitmap', '-noui'));
    uimenu(cmenu, 'Label', 'Save as PNG', 'Callback', @(~,~) save_spectrum_figure(fh, ijk));
    set(gca, 'UIContextMenu', cmenu);
    
    drawnow;
end