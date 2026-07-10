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
%   The overlay is laid out [RO, PE, slice, T] (permute [dims.x, dims.y,
%   dims.t]) and flipped along dim 1 so it aligns with the spec2nii sform on
%   the T1 -- the same overlay orientation the integration panel uses, so the
%   two agree on where each voxel sits.
%   NOTE: the companion ftSpec is stored 180 deg rotated relative to this
%   overlay (S_time = ccav_w and S_freq = ftSpec_smooth come from different
%   recon branches), so nii_viewer's click handler compensates when it pulls a
%   voxel spectrum.  See show_mrsi_spectrum_onclick.m.

    if nargin < 5, options = struct(); end
    if ~isfield(options, 'mode'), options.mode = 'signal'; end

    fprintf('\n=== mrsi_on_t1_map (mode: %s) ===\n', upper(options.mode));

    dSpec_time = S_time.dims.t;
    nii_empty  = nii_tool('load', emptyNiiPath);

    switch lower(options.mode)
    case 'signal'
        %% ------------------ SIGNAL MODE: 4D time-domain --------------------
        % Lay the data out as [RO, PE, 1, T] (permute [dims.x, dims.y, dims.t])
        % then flip dim 1 so it aligns with the spec2nii sform on the T1.  This
        % is the same overlay orientation the integration panel uses, so the
        % two agree on where each voxel sits.
        vol = permute(S_time.data, [S_time.dims.x, S_time.dims.y, dSpec_time]);
        [X, Y, T] = size(vol);
        vol = single(reshape(vol, [X, Y, 1, T]));  % nii_tool keeps complex64
        vol = flip(vol, 1);

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

        met2d = options.map.(metName);              % [Y, X]
        vol   = single(reshape(met2d.', [size(met2d,2), size(met2d,1), 1, 1]));
        vol   = flip(vol, 1);                        % same dim-1 flip as signal mode
        [X, Y, ~, ~] = size(vol);

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
