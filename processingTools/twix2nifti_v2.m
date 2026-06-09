function twix2nifti_v2(twix_file, output_file)
% TWIX2NIFTI_V2  Convert a Siemens TWIX (.dat) MRSI scan to an empty
%                NIfTI-MRS file matching `spec2nii twix -e image`.
%   Usage
%       twix2nifti_v2('meas_MID00123.dat')
%       twix2nifti_v2('meas_MID00123.dat', 'output.nii')

    if nargin < 1 || isempty(twix_file)
        error('twix2nifti_v2: twix file path is required.');
    end
    if nargin < 2 || isempty(output_file)
        [~, fname, ~] = fileparts(twix_file);
        output_file = [fname '_empty.nii'];
    end
    if ~endsWith(output_file, '.nii', 'IgnoreCase', true)
        output_file = [output_file '.nii'];
    end

    %Read TWIX header
    fprintf('Loading twix file: %s\n', twix_file);
    twixObj = mapVBVD(twix_file);
    if iscell(twixObj)
        twixObj = twixObj{end};
    end
    hdr  = twixObj.hdr;
    yaps = hdr.MeasYaps;

    %Matrix size  (lBaseResolution priority)
    fprintf('\n=== Resolution fields ===\n');
    n_pe = pick_matrix_size(hdr.Meas, 'lBaseResolution', 'lFinalMatrixSizePhase', 'PE   ');
    n_ro = pick_matrix_size(hdr.Meas, 'lBaseResolution', 'lFinalMatrixSizeRead',  'RO   ');
    n_sl = pick_matrix_size(hdr.Meas, '',                'lFinalMatrixSizeSlice', 'SLICE');

    if isfield(hdr.Meas, 'lVectorSize') && ~isempty(hdr.Meas.lVectorSize)
        n_t = double(hdr.Meas.lVectorSize);
    else
        n_t = 1;
    end
    fprintf('  TIME : lVectorSize = %d\n', n_t);

    data_size = [n_pe, n_ro, n_sl];
    fprintf('  Spatial matrix [PE RO SL] = [%d %d %d]\n', data_size);

    %Dwell time   ()
    alDwellNs = yaps_get(yaps, {{'sRXSPEC','alDwellTime',1}}, 0.0);
    dwellTime = alDwellNs / 1e9 * 2;
    fprintf('  alDwellTime = %.0f ns -> pixdim[4] = %.6f s (remove_os doubled)\n', ...
            alDwellNs, dwellTime);

    %Slice / VOI orientation parameters from MeasYaps
    NormaldSag = yaps_get(yaps, {{'sSliceArray','asSlice',1,'sNormal','dSag'}, ...
                                  {'sSpecPara','sVoI','sNormal','dSag'}}, 0.0);
    NormaldCor = yaps_get(yaps, {{'sSliceArray','asSlice',1,'sNormal','dCor'}, ...
                                  {'sSpecPara','sVoI','sNormal','dCor'}}, 0.0);
    NormaldTra = yaps_get(yaps, {{'sSliceArray','asSlice',1,'sNormal','dTra'}, ...
                                  {'sSpecPara','sVoI','sNormal','dTra'}}, 0.0);
    inplaneRot = yaps_get(yaps, {{'sSliceArray','asSlice',1,'dInPlaneRot'}, ...
                                  {'sSpecPara','sVoI','dInPlaneRot'}}, 0.0);
    RoFoV = yaps_get(yaps, {{'sSliceArray','asSlice',1,'dReadoutFOV'}, ...
                             {'sSpecPara','sVoI','dReadoutFOV'}}, 10000.0);
    PeFoV = yaps_get(yaps, {{'sSliceArray','asSlice',1,'dPhaseFOV'}, ...
                             {'sSpecPara','sVoI','dPhaseFOV'}}, 10000.0);
    sliceThickness = yaps_get(yaps, {{'sSliceArray','asSlice',1,'dThickness'}, ...
                                      {'sSpecPara','sVoI','dThickness'}}, 10000.0);
    PosdSag = yaps_get(yaps, {{'sSliceArray','asSlice',1,'sPosition','dSag'}, ...
                               {'sSpecPara','sVoI','sPosition','dSag'}}, 0.0);
    PosdCor = yaps_get(yaps, {{'sSliceArray','asSlice',1,'sPosition','dCor'}, ...
                               {'sSpecPara','sVoI','sPosition','dCor'}}, 0.0);
    PosdTra = yaps_get(yaps, {{'sSliceArray','asSlice',1,'sPosition','dTra'}, ...
                               {'sSpecPara','sVoI','sPosition','dTra'}}, 0.0);

    PosdSag = PosdSag + yaps_get(yaps, {{'lScanRegionPosSag'}}, 0.0);
    PosdCor = PosdCor + yaps_get(yaps, {{'lScanRegionPosCor'}}, 0.0);
    PosdTra = PosdTra + yaps_get(yaps, {{'lScanRegionPosTra'}}, 0.0);

    sliceNormal = [NormaldSag, NormaldCor, NormaldTra];
    if ~any(sliceNormal); sliceNormal(1) = 1.0; end

    fprintf('\n=== Orientation parameters ===\n');
    fprintf('  Slice normal     : [% .4f % .4f % .4f]\n', sliceNormal);
    fprintf('  In-plane rot     :  % .4f rad\n',          inplaneRot);
    fprintf('  FOV  RO / PE     :  %.2f / %.2f mm\n',     RoFoV, PeFoV);
    fprintf('  Slice thickness  :  %.2f mm\n',            sliceThickness);
    fprintf('  Position (S/C/T) : [% .2f % .2f % .2f] mm\n', PosdSag, PosdCor, PosdTra);

    %CSI orientation
    base_pos = [PosdSag, PosdCor, PosdTra];
    [iop, ipp, pixSpc, sliceThickness, dim_swapped] = csi_orientations( ...
        sliceNormal, inplaneRot, PeFoV, RoFoV, sliceThickness, ...
        n_pe, n_ro, n_sl, base_pos);

    fprintf('\n=== Orientation results ===\n');
    fprintf('  imagePositionPatient : [% .4f % .4f % .4f]\n', ipp);
    fprintf('  imageOrientationPatient (row, col):\n');
    fprintf('    row = [% .4f % .4f % .4f]\n', iop(1, :));
    fprintf('    col = [% .4f % .4f % .4f]\n', iop(2, :));
    fprintf('  pixelSpacing      : [%.4f %.4f] mm\n',  pixSpc);
    fprintf('  sliceThickness    : %.4f mm\n',         sliceThickness);
    fprintf('  dim_swapped       : %d\n',              dim_swapped);

    %dcm2niix-style 4x4 affine + half-voxel shift
    xyzMM = [pixSpc(:); sliceThickness];
    Q44   = nifti_dicom2mat(iop, ipp, xyzMM);
    Q44   = verify_slice_dir(Q44, [n_pe, n_ro, n_sl]);
    Q44(1:2, :) = -Q44(1:2, :);
    half = [0.5, 0.5, 0] * Q44(1:3, 1:3).';
    Q44(1:3, 4) = Q44(1:3, 4) + half(:);

    fprintf('\n=== Final affine (Q44) ===\n');
    disp(Q44);

    %Apply transverse PE/RO data-shape swap
    if dim_swapped
        data_size = [data_size(2), data_size(1), data_size(3)];
        fprintf('  Transverse: PE/RO swapped -> [%d %d %d]\n', data_size);
    end

    %Quaternion (qform) from Q44
    [qb, qc, qd, qfac, qoff] = mat44_to_quatern(Q44);

    %Voxel pixel sizes
    pixdim_xyz = zeros(1, 3);
    for c = 1:3
        v = norm(Q44(1:3, c));
        if v == 0; v = 1; end
        pixdim_xyz(c) = v;
    end

    %Write 4-D complex128 NIfTI-MRS-shaped empty file
    full_size = [data_size, n_t];
    write_empty_nifti_mrs(output_file, full_size, Q44, pixdim_xyz, ...
                          dwellTime, qb, qc, qd, qfac, qoff);

    fprintf('\nWrote empty NIfTI-MRS: %s\n', output_file);
    fprintf('  shape    : [%d %d %d %d]\n', full_size);
    fprintf('  datatype : 1792 (complex128, bitpix 128)\n');
    fprintf('  pixdim   : [%.4f %.4f %.4f %.6f] (mm,mm,mm,s)\n', ...
            pixdim_xyz, dwellTime);
    fprintf('  sform_code = 2 / qform_code = 2\n');
    fprintf('Done.\n');
end


% =========================================================================
% Helpers
% =========================================================================
function n = pick_matrix_size(meas, primary_field, fallback_field, label) %%this function is similar to the thingy i modified for base resolution or fallback to interpolated res
    n   = 1;
    src = 'default 1';
    if ~isempty(primary_field) && isfield(meas, primary_field) ...
            && ~isempty(meas.(primary_field)) ...
            && double(meas.(primary_field)) > 0
        n   = double(meas.(primary_field));
        src = sprintf('%s = %d', primary_field, n);
    elseif isfield(meas, fallback_field) ...
            && ~isempty(meas.(fallback_field)) ...
            && double(meas.(fallback_field)) > 0
        n   = double(meas.(fallback_field));
        src = sprintf('%s = %d', fallback_field, n);
    end
    fprintf('  %s: %s\n', label, src);
end


function val = yaps_get(yaps, paths, default_val) %%the long pymapVBVD if structure to find
    for ii = 1:numel(paths)
        [ok, v] = walk(yaps, paths{ii});
        if ok
            val = double(v);
            return;
        end
    end
    val = default_val;
end


function [ok, v] = walk(s, path) %%this is that dictionary reference table 
    ok = true;  v = s;
    for k = 1:numel(path)
        p = path{k};
        if isnumeric(p)
            if iscell(v) && numel(v) >= p && ~isempty(v{p})
                v = v{p};
            elseif isstruct(v) && numel(v) >= p
                v = v(p);
            else
                ok = false; v = []; return;
            end
        else
            if isstruct(v) && isfield(v, p)
                v = v.(p);
            else
                ok = false; v = []; return;
            end
        end
        if isempty(v); ok = false; v = []; return; end
    end
end


function ori = class_ori(sag, cor, tra) %%the same thing i did in io_load_twix2 to find out which plane dominates sag corr or trans
    as = abs(sag); ac = abs(cor); at = abs(tra);
    eq_sc = abs(as - ac) < 1e-12;
    eq_st = abs(as - at) < 1e-12;
    eq_ct = abs(ac - at) < 1e-12;
    if (eq_sc && eq_st) || (eq_sc && as < at) || (eq_st && as > ac) ...
            || (eq_ct && ac > as) || (as > ac && as < at) ...
            || (as < ac && ac < at) || (as < at && at > ac) ...
            || (ac < at && at > as)
        ori = 2;
    elseif (eq_sc && as > at) || (eq_st && as < ac) ...
            || (as < ac && ac > at) || (as > at && as < ac) ...
            || (as < at && at < ac)
        ori = 1;
    elseif (eq_ct && ac < as) || (as > ac && as > at) ...
            || (ac > at && ac < as) || (ac < at && at < as)
        ori = 0;
    else
        error('class_ori: invalid slice orientation');
    end
end

% given the slice-normal gs and an in-plane rotation phi, return:
%   - gp — phase-encode direction unit vector
%   - gr — readout direction unit vector
% 
%   Algorithm:
%   1. Determine orientation case via class_ori.
%   2. Pick a default gp perpendicular to gs, choosing the formula by case:
%     - Transverse: gp = (0, gs_z, -gs_y) / sqrt(gs_y² + gs_z²)
%     - Coronal: gp = (gs_y, -gs_x, 0) / sqrt(gs_x² + gs_y²)
%     - Sagittal: gp = (-gs_y, gs_x, 0) / sqrt(gs_x² + gs_y²)
%   3. gr = gs × gp (right-handed orthogonal frame).
%   4. Rotate gp by angle phi around gs: gp' = cos(phi)·gp − sin(phi)·gr.
%   5. Recompute gr = gs × gp' so the frame stays orthogonal.
% step 5: rotating only gp would leave gr un-rotated, breaking orthogonality of the (gs, gp, gr) triad. The Siemens code re-derives gr from the cross product
%   after rotation

function [gp, gr] = calc_prs(gs, phi)
    SAGITTAL = 0; CORONAL = 1; TRANSVERSE = 2;
    gs = gs(:);
    ori = class_ori(gs(1), gs(2), gs(3));
    gp = zeros(3, 1);
    switch ori
        case TRANSVERSE
            d = sqrt(1.0 / (gs(2)^2 + gs(3)^2));
            gp(1) = 0; gp(2) = gs(3)*d; gp(3) = -gs(2)*d;
        case CORONAL
            d = sqrt(1.0 / (gs(1)^2 + gs(2)^2));
            gp(1) = gs(2)*d; gp(2) = -gs(1)*d; gp(3) = 0;
        case SAGITTAL
            d = sqrt(1.0 / (gs(1)^2 + gs(2)^2));
            gp(1) = -gs(2)*d; gp(2) = gs(1)*d; gp(3) = 0;
    end
    gr = [gs(2)*gp(3) - gs(3)*gp(2);
          gs(3)*gp(1) - gs(1)*gp(3);
          gs(1)*gp(2) - gs(2)*gp(1)];
    gp_rot = cos(phi)*gp - sin(phi)*gr;
    gp = gp_rot;
    gr = [gs(2)*gp(3) - gs(3)*gp(2);
          gs(3)*gp(1) - gs(1)*gp(3);
          gs(1)*gp(2) - gs(2)*gp(1)];
end


% convert Siemens orientation/FOV/position into DICOM-style:
%   - iop — imageOrientationPatient (2×3: row direction; column direction)
%   - ipp — imagePositionPatient (1×3: position of the first voxel)
%   - pixSpc — pixel spacing ([dy, dx] in DICOM order)
%   - fov_sl — slice thickness (per-partition for 3D MRSI)
%   - dim_swapped — boolean: do PE/RO need to be swapped in the data block?
% 
%   Walkthrough:
%   1. Determine orientation case (mo_case).
%   2. Get default gp, gr from calc_prs. dColVec = gp (column = phase), dRowVec = gr (row = readout).
%   3. 3D MRSI adjustment (n_sl > 1): Siemens encodes position as the centre of the slab; per-partition we want the centre of the first partition. So shift base_pos along
%   −sliceNormal by (fov_sl/2 − fov_sl/n_sl/2), and divide fov_sl by n_sl to get per-slice thickness.
%   4. Per-orientation handling:
%     - Sagittal: flip dRowVec (image-row mirror); pixel spacing = [fov_ro/n_ro, fov_pe/n_pe]; swap col/row vectors so the resulting IOP matches the DICOM convention; compute ipp as
%    the slab-centre minus half-extents.
%     - Coronal: same as sagittal but no row flip.
%     - Transverse: flip dRowVec; pixel spacing = [fov_pe/n_pe, fov_ro/n_ro] (note PE/RO order swapped here — that's why we set dim_swapped=true); compute ipp differently because
%   the ro/pe axes are swapped relative to the image grid.
%   5. Return iop = [dRowVec; dColVec].
% 
%   Why the "swap col/row" lines (tmp = dColVec; dColVec = dRowVec; dRowVec = tmp;): Siemens stores the slice frame with phase as the first in-plane axis, but DICOM's
%   imageOrientationPatient expects (row, column) = (readout, phase). For sagittal/coronal the swap converts between the two conventions; for transverse the swap is achieved instead
%    by setting dim_swapped=true and permuting the data array later in step 7 of the main routine.
function [iop, ipp, pixSpc, fov_sl, dim_swapped] = csi_orientations(...
        sliceNormal, ip_rot, fov_pe, fov_ro, fov_sl, n_pe, n_ro, n_sl, base_pos)
    mo_case  = class_ori(sliceNormal(1), sliceNormal(2), sliceNormal(3));
    [gp, gr] = calc_prs(sliceNormal, ip_rot);
    dColVec = gp.';
    dRowVec = gr.';
    if n_sl > 1
        base_pos = base_pos - sliceNormal(:).' * (fov_sl/2 - fov_sl/n_sl/2);
        fov_sl   = fov_sl / n_sl;
    end
    dim_swapped = false;
    switch mo_case
        case 0
            dRowVec = -dRowVec;
            pixSpc  = [fov_ro / n_ro, fov_pe / n_pe];
            tmp = dColVec; dColVec = dRowVec; dRowVec = tmp;
            ipp = base_pos - dRowVec*fov_pe/2 - dColVec*fov_ro/2;
        case 1
            pixSpc  = [fov_ro / n_ro, fov_pe / n_pe];
            tmp = dColVec; dColVec = dRowVec; dRowVec = tmp;
            ipp = base_pos - dRowVec*fov_pe/2 - dColVec*fov_ro/2;
        case 2
            dRowVec = -dRowVec;
            pixSpc  = [fov_pe / n_pe, fov_ro / n_ro];
            dim_swapped = true;
            ipp = base_pos - dRowVec*fov_ro/2 - dColVec*fov_pe/2;
    end
    iop = [dRowVec; dColVec];
end

% nifti_dicom2mat in spec2nii/dcm2niiOrientation/orientationFuncs.py (which ports dcm2niix/console/nifti1_io_core.cpp)
function Q44 = nifti_dicom2mat(orient, patientPosition, xyzMM)
    Q = zeros(3, 3);
    Q(1:2, :) = orient;
    for r = 1:3
        nrm = norm(Q(r, :)); if nrm == 0; nrm = 1; end
        Q(r, :) = Q(r, :) / nrm;
    end
    Q(3, :) = cross(Q(1, :), Q(2, :));
    Q = Q.';
    if det(Q) < 0; Q(:, 3) = -Q(:, 3); end
    xyzMM = xyzMM(:);
    xyzMM([1 2]) = xyzMM([2 1]);
    Q = Q * diag(xyzMM);
    Q44 = eye(4);
    Q44(1:3, 1:3) = Q;
    Q44(1:3, 4)   = patientPosition(:);
end


function R = verify_slice_dir(R, dim)
    if numel(dim) < 3 || dim(3) < 2; return; end
    iSL = 1;
    if abs(R(2,3)) >= abs(R(1,3)) && abs(R(2,3)) >= abs(R(3,3)); iSL = 2; end
    if abs(R(3,3)) >= abs(R(1,3)) && abs(R(3,3)) >= abs(R(2,3)); iSL = 3; end
    x = [0, 0, dim(3) - 1, 1];
    pos1v = x * R.';
    if pos1v(iSL) < R(iSL, 4); R(:, 3) = -R(:, 3); end
end


% --- Build a NIfTI-1 quaternion from a 4x4 affine ------------------------
function [qb, qc, qd, qfac, qoff] = mat44_to_quatern(M)
    qoff = M(1:3, 4);
    R = M(1:3, 1:3);

    sx = norm(R(:,1)); if sx > 0; R(:,1) = R(:,1)/sx; end
    sy = norm(R(:,2)); if sy > 0; R(:,2) = R(:,2)/sy; end
    sz = norm(R(:,3)); if sz > 0; R(:,3) = R(:,3)/sz; end

    z = cross(R(:,1), R(:,2));
    if dot(z, R(:,3)) < 0
        R(:,3) = -R(:,3);
        qfac = -1;
    else
        qfac = 1;
    end

    a = R(1,1) + R(2,2) + R(3,3) + 1;
    if a > 0.5
        a = 0.5 * sqrt(a);
        qb = 0.25 * (R(3,2) - R(2,3)) / a;
        qc = 0.25 * (R(1,3) - R(3,1)) / a;
        qd = 0.25 * (R(2,1) - R(1,2)) / a;
    else
        xd = 1 + R(1,1) - (R(2,2) + R(3,3));
        yd = 1 + R(2,2) - (R(1,1) + R(3,3));
        zd = 1 + R(3,3) - (R(1,1) + R(2,2));
        if xd > 1
            qb = 0.5 * sqrt(xd);
            qc = 0.25 * (R(1,2) + R(2,1)) / qb;
            qd = 0.25 * (R(1,3) + R(3,1)) / qb;
            a  = 0.25 * (R(3,2) - R(2,3)) / qb;
        elseif yd > 1
            qc = 0.5 * sqrt(yd);
            qb = 0.25 * (R(1,2) + R(2,1)) / qc;
            qd = 0.25 * (R(2,3) + R(3,2)) / qc;
            a  = 0.25 * (R(1,3) - R(3,1)) / qc;
        else
            qd = 0.5 * sqrt(zd);
            qb = 0.25 * (R(1,3) + R(3,1)) / qd;
            qc = 0.25 * (R(2,3) + R(3,2)) / qd;
            a  = 0.25 * (R(2,1) - R(1,2)) / qd;
        end
        if a < 0
            qb = -qb; qc = -qc; qd = -qd;
        end
    end
end


% --- Minimal NIfTI-1 (.nii) writer for 4-D complex128 empty MRS ---------
function write_empty_nifti_mrs(filename, full_size, Q44, pixdim_xyz, ...
                                dwellTime, qb, qc, qd, qfac, qoff)
    if numel(full_size) < 4
        full_size(end+1:4) = 1;
    end

    fid = fopen(filename, 'w', 'l');
    if fid < 0
        error('write_empty_nifti_mrs: cannot open %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    fwrite(fid, 348,                'int32');        % sizeof_hdr
    fwrite(fid, zeros(1, 10),       'uint8');        % data_type[10]
    fwrite(fid, zeros(1, 18),       'uint8');        % db_name[18]
    fwrite(fid, 0,                  'int32');        % extents
    fwrite(fid, 0,                  'int16');        % session_error
    fwrite(fid, uint8('r'),         'uint8');        % regular
    fwrite(fid, 0,                  'uint8');        % dim_info

    % dim[8]
    dim    = ones(1, 8);
    dim(1) = 4;
    dim(2:5) = full_size(1:4);
    fwrite(fid, dim,                'int16');

    fwrite(fid, [0 0 0],            'float32');      % intent_p1,2,3
    fwrite(fid, 0,                  'int16');        % intent_code
    fwrite(fid, 1792,               'int16');        % datatype = COMPLEX128
    fwrite(fid, 128,                'int16');        % bitpix
    fwrite(fid, 0,                  'int16');        % slice_start

    % pixdim[8]: [qfac vox_x vox_y vox_z dwellTime 1 1 1]
    pix    = ones(1, 8);
    pix(1) = qfac;
    pix(2:4) = pixdim_xyz(1:3);
    pix(5)   = dwellTime;
    fwrite(fid, pix,                'float32');

    fwrite(fid, 352,                'float32');      % vox_offset
    fwrite(fid, 1,                  'float32');      % scl_slope
    fwrite(fid, 0,                  'float32');      % scl_inter
    fwrite(fid, 0,                  'int16');        % slice_end
    fwrite(fid, 0,                  'uint8');        % slice_code
    fwrite(fid, 2 + 8,              'uint8');        % xyzt_units (mm + sec)
    fwrite(fid, [0 0],              'float32');      % cal_max, cal_min
    fwrite(fid, 0,                  'float32');      % slice_duration
    fwrite(fid, 0,                  'float32');      % toffset
    fwrite(fid, [0 0],              'int32');        % glmax, glmin

    descrip     = 'spec2nii twix MATLAB port (NIfTI-MRS)';
    descrip_buf = zeros(1, 80, 'uint8');
    descrip_buf(1:numel(descrip)) = uint8(descrip);
    fwrite(fid, descrip_buf,        'uint8');
    fwrite(fid, zeros(1, 24),       'uint8');        % aux_file[24]

    fwrite(fid, 2,                  'int16');        % qform_code = 2 (Aligned)
    fwrite(fid, 2,                  'int16');        % sform_code = 2 (Aligned)

    fwrite(fid, qb,                 'float32');      % quatern_b
    fwrite(fid, qc,                 'float32');      % quatern_c
    fwrite(fid, qd,                 'float32');      % quatern_d
    fwrite(fid, qoff(1),            'float32');      % qoffset_x
    fwrite(fid, qoff(2),            'float32');      % qoffset_y
    fwrite(fid, qoff(3),            'float32');      % qoffset_z

    fwrite(fid, Q44(1, :),          'float32');      % srow_x[4]
    fwrite(fid, Q44(2, :),          'float32');      % srow_y[4]
    fwrite(fid, Q44(3, :),          'float32');      % srow_z[4]

    fwrite(fid, zeros(1, 16),       'uint8');        % intent_name[16]
    fwrite(fid, uint8('n+1'),       'uint8');        % magic[3]
    fwrite(fid, 0,                  'uint8');        % null
    fwrite(fid, [0 0 0 0],          'uint8');        % pad to 352

    % Voxel data: complex128 = 16 bytes per voxel, all zeros.
    % Write in chunks to avoid allocating the full array.
    n_vox    = prod(double(full_size));
    n_doubles = n_vox * 2;                           % real + imag per voxel
    chunk    = 1e6;                                  % ~8 MB per chunk
    written  = 0;
    while written < n_doubles
        nw = min(chunk, n_doubles - written);
        fwrite(fid, zeros(nw, 1, 'double'), 'double');
        written = written + nw;
    end

    delete(cleaner);
end
