function S = get_sform_3x4(H)
if isfield(H,'sform_mat') && ~isempty(H.sform_mat)
    S = double(H.sform_mat);
else
    S = [double(H.srow_x); double(H.srow_y); double(H.srow_z)];
end
end