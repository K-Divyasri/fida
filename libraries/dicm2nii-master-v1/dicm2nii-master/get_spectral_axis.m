function x_axis = get_spectral_axis(p)
    % Generate spectral axis
    nPoints = size(p.nii.img, 4);
    fprintf('Creating spectral axis for %d points\n', nPoints);
    
    if isfield(p, 'ppm')
        fprintf('Using p.ppm field\n');
        x_axis = p.ppm; return;
    elseif isfield(p.nii, 'ppm')
        fprintf('Using p.nii.ppm field\n');
        x_axis = p.nii.ppm; return;
    end
    
    if isfield(p, 'spectralTime')
        fprintf('Using p.spectralTime field\n');
        x_axis = p.spectralTime; return;
    end
    
    TR = p.nii.hdr.pixdim(5);
    fprintf('TR value: %.6f\n', TR);
    if TR > 0
        fprintf('Using TR-based axis\n');
        x_axis = (0:nPoints-1) * TR; return;
    end
    
    temporal_unit = bitand(p.nii.hdr.xyzt_units, 56);
    fprintf('Temporal unit: %d\n', temporal_unit);
    if temporal_unit == 40 % ppm
        fprintf('Using ppm axis (0-10)\n');
        x_axis = linspace(0, 10, nPoints);
    elseif temporal_unit == 32 % Hz
        fprintf('Using Hz axis (0-2000)\n');
        x_axis = linspace(0, 2000, nPoints);
    else
        fprintf('Using default point-based axis\n');
        x_axis = 1:nPoints;
    end
end