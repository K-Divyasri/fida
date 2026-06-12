function out = mrsi_on_t1_map(S_time, S_freq, outPath, emptyNiiPath, options)
% MRSI_ON_T1_MAP  Pour a FID-A time-domain CSI struct into a spec2nii empty
% NIfTI template so nii_viewer can overlay it on a T1.
%
% INPUTS
%   S_time       - FID-A time-domain CSI struct (e.g. ccav_w).  Must have
%                  dims.t, dims.x, dims.y set, complex data of shape
%                  [t, ..., y, x, ...].
%   S_freq       - FID-A frequency-domain CSI struct (saved as companion
%                  .mat so the boxed/integrating viewers can find it).
%                  nii_viewer itself does NOT read this file.
%   outPath      - Output path for the filled NIfTI (.nii or .nii.gz).
%   emptyNiiPath - spec2nii empty NIfTI template (sets the sform).
%   options      - struct, optional fields:
%                    .mode  'signal' (default) | 'metabolite'
%                    .map, .crlb, .LW, .SNR, .metabolite_name (for metabolite mode)
%
% ORIENTATION NOTE
%   The spatial layout is built with the SAME transform that
%   mrsi_integration_panel uses for its (verified-perfect) T1-aligned NIfTI:
%       img = flip(mp.', 1)        where mp is (Ny, Nx) = (iy, ix)
%   i.e. arrange the data as (iy, ix), transpose to (ix, iy), then flip the
%   ix axis.  Signal mode applies this per time point; metabolite mode applies
%   it to the 2-D map.  Because the writer and the integration panel use the
%   identical transform, nii_viewer and the panel overlay agree exactly on
%   the T1.  No additional flipud/fliplr/rot90 is needed.

    if nargin < 5, options = struct(); end
    if ~isfield(options, 'mode'), options.mode = 'signal'; end

    fprintf('\n=== mrsi_on_t1_map (mode: %s) ===\n', upper(options.mode));

    dSpec_time = S_time.dims.t;
    nii_empty  = nii_tool('load', emptyNiiPath);

    switch lower(options.mode)
    case 'signal'
        %% ------------------ SIGNAL MODE: 4D time-domain --------------------
        % Lay out each spatial slice with the EXACT transform that
        % mrsi_integration_panel uses for its (verified-perfect) T1-aligned
        % NIfTI:  img = flip(mp.', 1)  where mp is (Ny, Nx) = (iy, ix).
        % We apply that same transform per time point so the 4D overlay lands
        % on the T1 in exactly the panel's orientation.
        %   1. arrange data as (iy, ix, t)
        %   2. transpose iy<->ix, then flip dim 1 (ix)  ==  flip(slice.',1) per t
        arr = permute(S_time.data, [S_time.dims.y, S_time.dims.x, dSpec_time]); % (Ny, Nx, T)
        [Y, X, T] = size(arr);                    % Y = Ny (iy), X = Nx (ix)
        vol = flip(permute(arr, [2 1 3]), 1);     % (X, Y, T): vol(a,b,t)=data(ix=X-a+1, iy=b)
        vol = single(reshape(vol, [X, Y, 1, T])); % nii_tool keeps complex64

        % Sync the empty NIfTI's time dimension to whatever the data actually has
        T_expected = nii_empty.hdr.dim(5);
        if T ~= T_expected
            fprintf('  Adjusting NIfTI T: empty=%d -> data=%d\n', T_expected, T);
            nii_empty.hdr.dim(5)    = int16(T);
            nii_empty.hdr.pixdim(5) = single(S_time.adcDwellTime);
        end

        % Sanity check: spatial dims should match.  If they don't, something
        % upstream is wrong -- bail loudly rather than silently transposing.
        if nii_empty.hdr.dim(2) ~= X || nii_empty.hdr.dim(3) ~= Y
            error(['Spatial dim mismatch.\n  Empty NIfTI: [%d %d]\n  ' ...
                   'Data (after permute [dims.x dims.y dims.t]): [%d %d].\n' ...
                   'Check that the spec2nii empty was generated from this ' ...
                   'same .dat file and that dims.x/dims.y are set on S_time.'], ...
                nii_empty.hdr.dim(2), nii_empty.hdr.dim(3), X, Y);
        end

        % Fill and save
        nii_filled         = nii_empty;
        nii_filled.img     = vol;
        nii_filled.hdr.descrip = sprintf('MRSI time-domain [%dx%dx1x%d]', X, Y, T);
        nii_tool('save', nii_filled, outPath);
        fprintf('  Wrote 4D NIfTI: %s\n', outPath);

        % Companion .mat (untouched, byte-identical to caller's S_freq).
        % Only the boxed viewers use this; nii_viewer ignores it.
        [p, n, ~] = fileparts(outPath);
        if endsWith(n, '.nii'), n = extractBefore(n, '.nii'); end
        ftSpecPath = fullfile(p, [n '_ftSpec.mat']);
        ftSpec_smooth_w = S_freq;
        save(ftSpecPath, 'ftSpec_smooth_w', '-v7.3');
        fprintf('  Wrote companion ftSpec: %s\n', ftSpecPath);

        out = struct('mrsi4D_time', outPath, ...
                     'ftSpec_mat',  ftSpecPath, ...
                     'mode',        'signal');

    case 'metabolite'
        %% ------------------ METABOLITE MODE: 3D conc map -------------------
        % map.(metName) is [Y, X] (PE x RO).  Transpose to [RO, PE] so it
        % lands in the NIfTI with the same orientation as the signal-mode 4D.
        assert(isfield(options, 'map') && isfield(options, 'metabolite_name'), ...
               'metabolite mode requires options.map and options.metabolite_name');
        metName = options.metabolite_name;
        assert(isfield(options.map, metName), 'metabolite %s not in options.map', metName);

        met2d = options.map.(metName);              % [Y, X] = (iy, ix), same layout as panel's mp
        img2d = flip(met2d.', 1);                    % == integration panel saveNifti: flip(mp.',1)
        [X, Y] = size(img2d);                        % X = Nx (ix), Y = Ny (iy)
        vol    = single(reshape(img2d, [X, Y, 1, 1]));

        nii_filled = nii_empty;
        nii_filled.hdr.dim(1)    = int16(3);
        nii_filled.hdr.dim(2:4)  = int16([X, Y, 1]);
        nii_filled.hdr.dim(5)    = int16(1);
        nii_filled.img           = vol;
        nii_filled.hdr.descrip   = sprintf('MRSI metabolite: %s [%dx%dx1]', metName, X, Y);

        nii_tool('save', nii_filled, outPath);
        fprintf('  Wrote 3D metabolite NIfTI: %s\n', outPath);

        % Quality companions
        [p, n, ~] = fileparts(outPath);
        if endsWith(n, '.nii'), n = extractBefore(n, '.nii'); end
        qualityPath = fullfile(p, [n '_quality.mat']);
        quality_data = struct( ...
            'metabolite_name', metName, ...
            'crlb',  options.crlb.(metName).', ...   % match the .' on the map
            'LW',    options.LW.', ...
            'SNR',   options.SNR.', ...
            'map',   options.map, ...
            'crlb_full', options.crlb);
        save(qualityPath, 'quality_data', '-v7.3');

        ftSpecPath      = fullfile(p, [n '_ftSpec.mat']);
        ftSpec_smooth_w = S_freq;
        save(ftSpecPath, 'ftSpec_smooth_w', '-v7.3');

        out = struct('metabolite_nii',  outPath, ...
                     'quality_mat',     qualityPath, ...
                     'ftSpec_mat',      ftSpecPath, ...
                     'mode',            'metabolite', ...
                     'metabolite_name', metName);

    otherwise
        error('Unknown mode "%s". Use ''signal'' or ''metabolite''.', options.mode);
    end

    fprintf('=== mrsi_on_t1_map: done ===\n\n');
end
