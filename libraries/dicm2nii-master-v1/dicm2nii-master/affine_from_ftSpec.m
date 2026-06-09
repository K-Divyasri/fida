function A = affine_from_ftSpec(S, AssumeLPS, UseCoordsForOrigin, SwapXY, FlipX, FlipY)
% Build candidate MRSI RAS affine from ftSpec fields and options
A = S.affineMatrix;
if AssumeLPS, A = diag([-1 -1 1 1]) * A; end % LPS→RAS (Siemens)
% use coordinates for first-voxel center if provided
if UseCoordsForOrigin && isfield(S,'coordinates') && isfield(S.coordinates,'x') && isfield(S.coordinates,'y')
    ox = S.coordinates.x(1); oy = S.coordinates.y(1);
    if AssumeLPS, ox = -ox; oy = -oy; end
    A(1:3,4) = [ox; oy; A(3,4)];
end
% in-plane ops
if SwapXY, A(:,[1 2]) = A(:,[2 1]); end
if FlipX,  A(:,1) = -A(:,1); end
if FlipY,  A(:,2) = -A(:,2); end
end