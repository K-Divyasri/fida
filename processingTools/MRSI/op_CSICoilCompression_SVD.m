function [out, W] = op_CSICoilCompression_SVD(in, nComp, useVoxelWise)
% op_CSICoilCompression_SVD  Reduce coil count by SVD compression
%
% [out, W] = op_CSICoilCompression_SVD(in, nComp)
% [out, W] = op_CSICoilCompression_SVD(in, nComp, useVoxelWise)
%
% INPUTS
%   in           – FID-A CSI struct (has .data or .fids, .sz, .dims.*)
%   nComp        – number of virtual channels to keep (≤ physical coils)
%   useVoxelWise – optional (false).  true → SVD per voxel (slow)
%
% OUTPUTS
%   out – same struct but with coil axis length = nComp
%   W   – Nc×nComp weight matrix (global) or  Nc×nComp×Nvox  (voxel-wise)
%
% ----------------------------------------------------------------------
% 2025-06-19  DK  •   fixed inverse-permutation bug (ORDER length mismatch)
% ----------------------------------------------------------------------

if nargin < 3, useVoxelWise = false; end
assert(isstruct(in), 'Input must be a struct.');

% -------- pick raw array field -----------------------------------------
rawField = 'data';
if ~isfield(in,'data') && isfield(in,'fids'), rawField = 'fids'; end
A   = in.(rawField);                 % complex array
sz  = in.sz;
dim = in.dims;                       % dims struct

coilDim = dim.coils;
assert(coilDim > 0, 'dims.coils = 0 – no coil dimension present.');

Nc_phys = sz(coilDim);
assert(nComp <= Nc_phys, 'nComp (%d) exceeds # physical coils (%d).', ...
                          nComp, Nc_phys);

% -------- make coils last, time first -----------------------------------
allDims = 1:numel(sz);
perm    = [dim.t, setdiff(allDims,[dim.t coilDim],'stable'), coilDim];
A       = permute(A, perm);
spatSz  = size(A);                   % [Nt ... Nc]
Nt      = spatSz(1);
Nc      = spatSz(end);
Nvox    = prod(spatSz(2:end-1));

A       = reshape(A, Nt, Nvox, Nc);  % Nt × Nvox × Nc

% -------- compute SVD weights ------------------------------------------
if useVoxelWise
    W = zeros(Nc, nComp, Nvox, 'like', A);
    for v = 1:Nvox
        [U,~,~]   = svd(squeeze(A(:,v,:)).', 'econ');  % Nc×Nt → U Nc×Nc
        W(:,:,v)  = U(:,1:nComp);
    end
    B = zeros(Nt, Nvox, nComp, 'like', A);
    for v = 1:Nvox
        B(:,v,:) = squeeze(A(:,v,:)) * conj(W(:,:,v)); % (Nt×Nc)*(Nc×nComp)
    end
else
    D = reshape(A, Nt*Nvox, Nc).';   % Nc×(Nt·Nvox)
    [U,~,~] = svd(D, 'econ');
    W       = U(:,1:nComp);          % Nc×nComp
    % combine: Nt×Nc×Nvox  *  Nc×nComp  →  Nt×Nvox×nComp
    B       = pagemtimes(permute(A,[1 3 2]), conj(W));
    B       = permute(B, [1 3 2]);
end

% -------- reshape back to original layout ------------------------------
spatSz(end) = nComp;                 % replace Nc by nComp
B = reshape(B, spatSz);              % Nt × ... × nComp

% inverse permutation (no element dropped now!)
invPerm          = zeros(1, numel(perm));
invPerm(perm)    = 1:numel(perm);
B                = permute(B, invPerm);

% -------- pack output ---------------------------------------------------
out            = in;
out.(rawField) = B;
out.sz(coilDim)= nComp;
out.coilCompression.W     = W;
out.coilCompression.nComp = nComp;
out.flags.coilCompressed  = true;
end
