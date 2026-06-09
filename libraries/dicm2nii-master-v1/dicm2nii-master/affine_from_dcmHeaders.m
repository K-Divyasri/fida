function A = affine_from_dcmHeaders(dcmMat)
% Build T1 RAS affine from dicm2nii's dcmHeaders.mat (IOP/IPP/PixelSpacing)
X = load(dcmMat); fn = fieldnames(X); H0 = X.(fn{1});
if isstruct(H0) && numel(fieldnames(H0))==1, H0 = H0.(fieldnames(H0){1}); end
files = fieldnames(H0);
IOP=[]; pxsp=[]; IPP=[]; for k=1:numel(files)
    D = H0.(files{k});
    if isempty(IOP)  && isfield(D,'ImageOrientationPatient'), IOP=double(D.ImageOrientationPatient(:)); end
    if isempty(pxsp) && isfield(D,'PixelSpacing'),           pxsp=double(D.PixelSpacing(:));           end
    if isempty(IPP)  && isfield(D,'ImagePositionPatient'),   IPP=double(D.ImagePositionPatient(:));    end
end
assert(~isempty(IOP)&&~isempty(pxsp)&&~isempty(IPP), 'Missing IOP/IPP/PixelSpacing in dcmHeaders.');
% DICOM LPS basis → RAS
M = diag([-1 -1 1]);
ex_lps = IOP(1:3);  ey_lps = IOP(4:6);  ez_lps = cross(ex_lps, ey_lps);
ex = M*ex_lps; ey = M*ey_lps; ez = M*ez_lps;  o = M*IPP;
dx = pxsp(1); dy = pxsp(2);
A      = eye(4);
A(1:3,1) = ex*dx; A(1:3,2) = ey*dy; A(1:3,3) = ez; A(1:3,4) = o;
end