function A = affine_from_hdr(H)
% 4x4 RAS affine from dicm2nii header (sform preferred, qform fallback)
if isfield(H,'sform_code') && H.sform_code>0
    if isfield(H,'sform_mat') && ~isempty(H.sform_mat)
        A = [double(H.sform_mat); 0 0 0 1];
    else
        A = [double(H.srow_x); double(H.srow_y); double(H.srow_z); 0 0 0 1];
    end
elseif isfield(H,'qform_code') && H.qform_code>0
    pix = double(H.pixdim(2:4));
    bcd = double(H.quatern_bcd(:).');
    qx=bcd(1); qy=bcd(2); qz=bcd(3); qw = sqrt(max(0,1-(qx*qx+qy*qy+qz*qz)));
    R = [1-2*(qy^2+qz^2), 2*(qx*qy - qz*qw), 2*(qx*qz + qy*qw);
         2*(qx*qy + qz*qw), 1-2*(qx^2+qz^2), 2*(qy*qz - qx*qw);
         2*(qx*qz - qy*qw), 2*(qy*qz + qx*qw), 1-2*(qx^2+qy^2)];
    A = eye(4); A(1:3,1:3) = R .* pix; A(1:3,4) = double(H.qoffset_xyz(:));
else
    error('No sform/qform in header.');
end
end