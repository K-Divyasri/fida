function label = get_axis_label(p)
    if isfield(p, 'ppm') || isfield(p.nii, 'ppm')
        label = 'Chemical Shift (ppm)';
    elseif bitand(p.nii.hdr.xyzt_units, 56) == 40
        label = 'Chemical Shift (ppm)';
    elseif bitand(p.nii.hdr.xyzt_units, 56) == 32
        label = 'Frequency (Hz)';
    elseif bitand(p.nii.hdr.xyzt_units, 56) == 8
        label = 'Time (seconds)';
    else
        label = 'Spectral Point';
    end
end