function [out, stats] = op_CSISubtract(A, B, varargin)
% op_CSISubtract
%   Subtract two CSI/MRSI structs (A - B) and keep A's metadata.
%   Works for 3D [t,y,x] (e.g., ftSpec) and 5D [t,coils,averages,x,y] (e.g., ftSpatial).
%
% Usage:
%   out = op_CSISubtract(ftSpec, ftSpec_1);
%   out = op_CSISubtract(ftSpatialA, ftSpatialB, 'coilIndex',1, 'averageIndex',1);
%
% Name-Value:
%   'coilIndex'     : scalar or [] (default [] -> use all if present)
%   'averageIndex'  : scalar or [] (default [] -> use all if present)
%   'precision'     : 'likeA' (default), 'double', or 'single'
%   'allowSqueeze'  : true (default). If true, drop trailing singleton dims before shape match.
%
% Output:
%   out   : struct like A, with .data = A - B (aligned)
%   stats : struct with simple RMSE metrics over full data
%
% Notes:
%   - If coils/averages exist and you pass an index, both A and B are sliced there
%     before subtraction (sizes must then match).
%   - If sizes still don’t match after optional slicing/squeezing, an error is thrown.

    ip = inputParser;
    ip.addParameter('coilIndex',    [], @(x) isempty(x) || (isscalar(x) && x>=1));
    ip.addParameter('averageIndex', [], @(x) isempty(x) || (isscalar(x) && x>=1));
    ip.addParameter('precision', 'likeA', @(s) ischar(s) || isstring(s));
    ip.addParameter('allowSqueeze', true, @islogical);
    ip.parse(varargin{:});
    opt = ip.Results;

    % --- 0) Basic checks
    assert(isstruct(A) && isfield(A,'data'), 'A must be a struct with .data');
    assert(isstruct(B) && isfield(B,'data'), 'B must be a struct with .data');

    % --- 1) Optional slicing for coils/averages (if those dims exist)
    Adata = A.data;
    Bdata = B.data;

    % helper to read dim index safely
    getdim = @(S,lab) (isfield(S,'dims') && isfield(S.dims,lab) && S.dims.(lab)>0) * S.dims.(lab);

    % slice coil
    cA = getdim(A,'coils'); cB = getdim(B,'coils');
    if ~isempty(opt.coilIndex) && (cA || cB)
        if cA, Adata = indexAlong(Adata, cA, opt.coilIndex); end
        if cB, Bdata = indexAlong(Bdata, cB, opt.coilIndex); end
    end

    % slice averages
    aA = getdim(A,'averages'); aB = getdim(B,'averages');
    if ~isempty(opt.averageIndex) && (aA || aB)
        if aA, Adata = indexAlong(Adata, aA, opt.averageIndex); end
        if aB, Bdata = indexAlong(Bdata, aB, opt.averageIndex); end
    end

    % --- 2) Optionally squeeze trailing singleton dims
    if opt.allowSqueeze
        Adata = squeezeSafe(Adata);
        Bdata = squeezeSafe(Bdata);
    end

    % --- 3) Final shape check (must match)
    if ~isequal(size(Adata), size(Bdata))
        error(['op_CSISubtract: size mismatch after optional slicing/squeezing.\n' ...
               'size(A) = %s, size(B) = %s'], mat2str(size(Adata)), mat2str(size(Bdata)));
    end

    % --- 4) Subtract with chosen precision
    clsA = class(Adata);
    switch lower(string(opt.precision))
        case "double"
            outData = double(Adata) - double(Bdata);
        case "single"
            outData = single(Adata) - single(Bdata);
        otherwise % 'likeA'
            outData = cast(Adata, clsA) - cast(Bdata, clsA);
    end

    % --- 5) Build output struct (inherit A’s metadata)
    out = A;
    out.data = outData;
    if isfield(out,'sz'), out.sz = int32(size(outData)); end

    % --- 6) Simple global RMSE metrics (over complex values)
    d = double(outData(:));
    stats.RMSE_global_complex = sqrt( mean( abs(d).^2 ) );
    stats.RMSE_global_real    = sqrt( mean( real(d).^2 ) );
end

% ---------- helpers ----------
function X = indexAlong(X, dim, idx)
    % Index along dimension 'dim' with 'idx' (keeps dimensionality where possible)
    sz = size(X);
    subs = repmat({':'}, 1, max(dim, ndims(X)));
    subs{dim} = idx;
    X = X(subs{:});
end

function Y = squeezeSafe(X)
    % Squeeze but keep at least 2D in spatial/CSI contexts
    Y = squeeze(X);
    if isvector(Y) % keep column shape (Nt x 1) instead of row if needed
        Y = Y(:);
    end
end
