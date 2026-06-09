function MRSIStruct = op_NUFFTSpatial(MRSIStruct, kFile, varargin)
% Minimal & readable NUFFT recon that matches your slow SFT path.
% Works for both dComp (5D) and dComp_w (4D) once reshapeDimensions is used.
%
% Example:
%   ftFast  = op_NUFFTSpatial(dComp,   kFile, 'Calib', true, 'nCalT', 12, 'Verbose', true);
%   ftFastW = op_NUFFTSpatial(dComp_w, kFile, 'Calib', true, 'nCalT', 12, 'Verbose', true);

% -------- options --------
p = inputParser;
p.addParameter('Calib',   true,  @(x)islogical(x)||isscalar(x));
p.addParameter('nCalT',   12,    @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Verbose', true,  @(x)islogical(x)||isscalar(x));
p.parse(varargin{:});
opt = p.Results;
vprintf = @(varargin) ifvprintf(opt.Verbose, varargin{:});
%% --------- Handle 4D data structure (no averages) -----------
if ndims(MRSIStruct.data) == 5
    %B = MRSIStruct;
    Nt = MRSIStruct.sz(MRSIStruct.dims.t);
    kPtsPerCycle = MRSIStruct.sz(MRSIStruct.dims.kpts);
    nCoil = MRSIStruct.sz(MRSIStruct.dims.coils);
    nAvg = MRSIStruct.sz(MRSIStruct.dims.averages);
    nKy = MRSIStruct.sz(MRSIStruct.dims.kshot);
    fprintf('Original size: %s\n', mat2str(size(MRSIStruct.data))); % e.g. [576 126 16 4 63]

    % Permute to [kpts t coils avg shot]
    A = permute(MRSIStruct.data, [4 1 2 3 5]);

    % Reshape to [kpts*t, coils, avg, shot]
    A = reshape(A, [kPtsPerCycle*Nt, nCoil, nAvg, nKy]);
    fprintf('Reshaped size: %s\n', mat2str(size(A)));

    MRSIStruct.data = A;
    MRSIStruct.sz   = double([size(A,1), size(A,2), size(A,3), size(A,4)]);
    %MRSIStruct = B;
    MRSIStruct.dims.kshot=0;
    MRSIStruct.dims.kpts=0;
    MRSIStruct.dims.ky=4;
    MRSIStruct.dims.t=1;

elseif ndims(MRSIStruct.data) == 4
    %B = MRSIStruct;
    Nt = MRSIStruct.sz(MRSIStruct.dims.t);
    kPtsPerCycle = MRSIStruct.sz(MRSIStruct.dims.kpts);
    nCoil = MRSIStruct.sz(MRSIStruct.dims.coils);
    %nAvg = MRSIStruct.sz(MRSIStruct.dims.averages);
    nKy = MRSIStruct.sz(MRSIStruct.dims.kshot);
    fprintf('Original size: %s\n', mat2str(size(MRSIStruct.data))); % e.g. [576 126 16 63]

    % Permute to [kpts t coils shot]
    A = permute(MRSIStruct.data, [3 1 2 4]); % [16 576 126 63]

    % Reshape to [kpts*t, coils, shot]
    A = reshape(A, [kPtsPerCycle*Nt, nCoil, nKy]);
    fprintf('Reshaped size: %s\n', mat2str(size(A)));

    MRSIStruct.data = A;
    MRSIStruct.sz   = double([size(A,1), size(A,2), size(A,3)]);
    MRSIStruct.dims.kshot=0;
    MRSIStruct.dims.kpts=0;
    MRSIStruct.dims.ky=3;
    MRSIStruct.dims.t=1;
    %MRSIStruct = B;
end

vprintf('=== START: op_NUFFTSpatial ===\n');

% -------- read k (1/mm) --------
[kTab, kXY] = readKFile_simple(kFile);       % kXY = [Kx Ky] in 1/mm
if isempty(kXY), error('Could not read Kx,Ky from "%s".', kFile); end
kXY = double(kXY(:,1:2));                    % [Nk x 2]
Nk  = size(kXY,1);

% -------- sizes like slow path --------
sz   = MRSIStruct.sz;
dims = MRSIStruct.dims;

Nt_total = sz(dims.t);        % e.g., 72576
nKy      = sz(dims.ky);       % e.g., 63
if mod(Nk, nKy)~=0, error('Nk (%d) not divisible by nKy (%d).', Nk, nKy); end
kPtsPerCycle = Nk / nKy;      % e.g., 126
if mod(Nt_total, kPtsPerCycle)~=0
    error('Nt_total (%d) not multiple of kPtsPerCycle (%d).', Nt_total, kPtsPerCycle);
end
NPtemporal = Nt_total / kPtsPerCycle;

xCoords = getCoordinates(MRSIStruct, 'x');   % mm
yCoords = getCoordinates(MRSIStruct, 'y');   % mm
Nx = numel(xCoords); Ny = numel(yCoords);

% reshape to [t ky extras] like SFT
[MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t','ky'});
X = getData(MRSIStruct);                      % [Nt_total  nKy  Nextra]
Nextra = size(X,3);                           % coils*averages, etc.

vprintf('Nk=%d, nKy=%d, kPtsPerCycle=%d, NPtemporal=%d, Nx=%d, Ny=%d, Nextra=%d\n', ...
        Nk, nKy, kPtsPerCycle, NPtemporal, Nx, Ny, Nextra);

% -------- NUFFT grid spacing & shift (from coords) --------
dx = median(abs(diff(xCoords(:))));
dy = median(abs(diff(yCoords(:))));
x_idx = (0:Nx-1).';  y_idx = (0:Ny-1).';
x_shift = median( x_idx - xCoords(:)/dx );    % ~Nx/2 if centered at 0
y_shift = median( y_idx - yCoords(:)/dy );    % ~Ny/2
n_shift = [y_shift, x_shift];                 % [Ny-shift, Nx-shift] order
Nd      = [Ny, Nx];
vprintf('Inferred n_shift = [%.3f, %.3f]\n', y_shift, x_shift);

% -------- SLOW SFT operator for calibration --------
% (x,y) in mm; phase = + i*2π*(x*Kx + y*Ky) / Nk  (adjoint-equivalent scale)
[xx, yy] = meshgrid(xCoords, yCoords);
XY  = [xx(:), yy(:)];                         % [Ny*Nx x 2]
SFT = exp(1i*2*pi * (XY * kXY.'));            % [Ny*Nx x Nk]
SFT = SFT / Nk;

% -------- NUFFT candidates (axis/sign) --------
% NUFFT needs radians-per-index; convert 1/mm -> rad using 2π*k*Δ
om_yx = [2*pi*kXY(:,2)*dy, 2*pi*kXY(:,1)*dx]; % [Ky,Kx] -> [y,x]
cands = {
  'om=[Ky,Kx] -> [y,x]',            om_yx
  'om=[Kx,Ky] -> [y,x]',            [2*pi*kXY(:,1)*dy, 2*pi*kXY(:,2)*dx]
  'om=[-Ky,+Kx]',                   [-om_yx(:,1),  +om_yx(:,2)]
  'om=[+Ky,-Kx]',                   [ +om_yx(:,1), -om_yx(:,2)]
  'om=[-Kx,+Ky]',                   [-2*pi*kXY(:,1)*dy, +2*pi*kXY(:,2)*dx]
  'om=[+Kx,-Ky]',                   [ +2*pi*kXY(:,1)*dy, -2*pi*kXY(:,2)*dx]
  'om=[-Ky,-Kx]',                   [-om_yx(:,1),  -om_yx(:,2)]
  'om=[-Kx,-Ky]',                   [-2*pi*kXY(:,1)*dy, -2*pi*kXY(:,2)*dx]
};

% calibration frames
cal_idx = 1:NPtemporal;
if opt.Calib
    nCal = min(opt.nCalT, NPtemporal);
    cal_idx = unique(max(1, round(linspace(1, NPtemporal, nCal))));
end

% -------- pick best candidate + complex scale alpha --------
Jd = [6,6];  Kd = 2*Nd;
best.resid = inf;
for ic = 1:size(cands,1)
    name = cands{ic,1};
    om   = cands{ic,2};
    st   = nufft_init(om, Nd, Jd, Kd, n_shift);

    % Least-squares alpha using a few t-frames
    num = 0+0i; den = 0; resid = 0;
    for it = cal_idx
        i0 = (it-1)*kPtsPerCycle + 1;     i1 = it*kPtsPerCycle;
        Y   = double(reshape(X(i0:i1, :, :), [], Nextra));     % [Nk x Nextra]
        Zs  = SFT * Y;                                        % [Ny*Nx x Nextra]

        Zn  = nufft_adj(Y, st);                                % [Ny x Nx x Nextra] or [Ny*Nx x Nextra]
        ZnM = as_mat(Zn, Ny, Nx, Nextra);                      % -> [Ny*Nx x Nextra]

        a = ZnM(:);  b = Zs(:);
        num = num + (a' * b);
        den = den + (a' * a);
    end
    alpha = (den~=0) * (num/den) + (den==0) * 1;

    % residual
    for it = cal_idx
        i0 = (it-1)*kPtsPerCycle + 1;     i1 = it*kPtsPerCycle;
        Y   = double(reshape(X(i0:i1, :, :), [], Nextra));
        Zs  = SFT * Y;                                        % [Ny*Nx x Nextra]

        Zn  = nufft_adj(Y, st);
        ZnM = as_mat(Zn, Ny, Nx, Nextra);                     % [Ny*Nx x Nextra]

        r   = alpha * ZnM - Zs;
        resid = resid + norm(r(:))^2;
    end

    vprintf('Cand %d: %s  | alpha=%+.6e%+.6ei  | resid=%.3e\n', ...
            ic, name, real(alpha), imag(alpha), resid);

    if resid < best.resid
        best.name  = name;
        best.st    = st;
        best.alpha = alpha;
        best.resid = resid;
    end
end
vprintf('Selected: %s  | alpha=%+.6e%+.6ei  | resid=%.3e\n', ...
        best.name, real(best.alpha), imag(best.alpha), best.resid);

% -------- full recon --------
img = zeros(NPtemporal, Ny, Nx, Nextra, 'like', double(X));
for it = 1:NPtemporal
    i0 = (it-1)*kPtsPerCycle + 1;     i1 = it*kPtsPerCycle;
    Y  = double(reshape(X(i0:i1, :, :), [], Nextra));          % [Nk x Nextra]

    Z  = best.alpha * nufft_adj(Y, best.st);                   % [Ny x Nx x Nextra] or [Ny*Nx x Nextra]
    Z3 = as_cube(Z, Ny, Nx, Nextra);                           % -> [Ny x Nx x Nextra]

    img(it,:,:,:) = Z3;
end

% -------- reshape back like slow SFT --------
MRSIStruct = setData(MRSIStruct, double(img));
kyDim = getDimension(MRSIStruct,'ky');
prevPermute = removeDimPrevPermute(prevPermute, kyDim);
prevPermute = addDimPrevPermute(prevPermute, 'y', kyDim);
prevPermute = addDimPrevPermute(prevPermute, 'x', kyDim+1);
prevSize(1) = NPtemporal; prevSize(2) = Ny;
prevSize    = [prevSize(1:2), Nx, prevSize(3:end)];
MRSIStruct  = reshapeBack(MRSIStruct, prevPermute, prevSize);

% spectral values identical to slow path
adcDT  = getAdcDwellTime(MRSIStruct);
specDT = adcDT * kPtsPerCycle;
specSW = 1 / specDT;
specT  = 0:specDT:specDT*(NPtemporal-1);
MRSIStruct = setSpectralWidth(MRSIStruct,  specSW);
MRSIStruct = setSpectralDwellTime(MRSIStruct, specDT);
MRSIStruct = setSpectralTime(MRSIStruct, specT);

MRSIStruct = setFlags(MRSIStruct,'spatialFT',true);
vprintf('=== END: op_NUFFTSpatial ===\n');
end

% ================= helpers =================
function M = as_mat(Z, Ny, Nx, Nextra)
% Return [Ny*Nx x Nextra]
    if ndims(Z)==3 && all(size(Z)==[Ny Nx Nextra])
        M = reshape(Z, Ny*Nx, Nextra);
    elseif ismatrix(Z)
        if size(Z,1)==Ny*Nx
            M = Z;
        elseif size(Z,1)==Ny && size(Z,2)==Nx*Nextra
            M = reshape(Z, Ny*Nx, Nextra);
        else
            M = reshape(Z, Ny*Nx, []);  % best-effort fallback
        end
    else
        M = reshape(Z, Ny*Nx, []);      % best-effort fallback
    end
end

function C = as_cube(Z, Ny, Nx, Nextra)
% Return [Ny x Nx x Nextra]
    if ndims(Z)==3 && all(size(Z)==[Ny Nx Nextra])
        C = Z;
    elseif ismatrix(Z) && size(Z,1)==Ny*Nx
        C = reshape(Z, Ny, Nx, []);
    elseif ismatrix(Z) && size(Z,1)==Ny
        C = reshape(Z, Ny, Nx, []);
    else
        C = reshape(Z, Ny, Nx, []);     % best-effort fallback
    end
end

function [kTable, kArray] = readKFile_simple(kFile)
    kTable = []; kArray = [];
    if isempty(kFile) || ~isfile(kFile), return; end
    try
        T = readtable(kFile);
        kTable = T;
        vn = lower(string(T.Properties.VariableNames));
        ix = find(vn=="kx",1); iy = find(vn=="ky",1);
        if ~isempty(ix) && ~isempty(iy)
            kArray = [T{:,ix}, T{:,iy}];
        elseif width(T)>=3
            kArray = [T{:,2}, T{:,3}];
        end
    catch
        A = readmatrix(kFile);
        if size(A,2)>=3, kArray = A(:,2:3); kTable = A; end
    end
end

function ifvprintf(flag, varargin)
    if flag, fprintf(varargin{:}); end
end


% function MRSIStruct = op_NUFFTSpatial(MRSIStruct, kFile, varargin)
% % NUFFT-based spatial recon that matches the slow SFT path on your data.
% % Usage:
% %   ftFast  = op_NUFFTSpatial(dComp,   kFile, 'Calib', true, 'nCalT', 12, 'Verbose', true);
% %   ftFastW = op_NUFFTSpatial(dComp_w, kFile, 'Calib', true, 'nCalT', 12, 'Verbose', true);
% 
% % -------- options --------
% p = inputParser;
% p.addParameter('Calib',   true,  @(x)islogical(x)||isscalar(x));
% p.addParameter('nCalT',   12,    @(x)isnumeric(x)&&isscalar(x)&&x>=1);
% p.addParameter('Verbose', true,  @(x)islogical(x)||isscalar(x));
% p.parse(varargin{:});
% opt = p.Results;
% 
% vprint = @(varargin) ifvprintf(opt.Verbose, varargin{:});
% vprint('=== START: op_NUFFTSpatial ===\n');
% 
% % -------- read k-trajectory (1/mm) --------
% [kTab, kXY] = readKFile_local(kFile);   % kXY: [Nk x 2] = [Kx, Ky] in 1/mm
% if isempty(kXY) || size(kXY,2) < 2
%     error('Failed to read Kx,Ky from "%s".', kFile);
% end
% kXY = double(kXY(:,1:2));
% Nk  = size(kXY,1);
% 
% % -------- sizes & reshape like SFT path --------
% sz   = MRSIStruct.sz;
% dims = MRSIStruct.dims;
% 
% it_t  = dims.t;   Nt_total = sz(it_t);     % 72576
% it_ky = dims.ky;  nKy      = sz(it_ky);    % 63
% 
% kPtsPerCycle = getKPtsPerCycle_local(kTab); % 126
% if mod(Nt_total, kPtsPerCycle) ~= 0
%     error('Nt_total (%d) not multiple of kPtsPerCycle (%d).', Nt_total, kPtsPerCycle);
% end
% NPtemporal = Nt_total / kPtsPerCycle;       % 576
% 
% % Target image grid (40x40 from coordinates)
% xCoords = getCoordinates(MRSIStruct, 'x');  % mm
% yCoords = getCoordinates(MRSIStruct, 'y');  % mm
% Nx = length(xCoords);
% Ny = length(yCoords);
% if Nx<=0 || Ny<=0
%     error('Empty coordinates.x / coordinates.y; cannot determine image size.');
% end
% 
% % Reshape to [t_total, ky, extras] (same as slow path)
% [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t','ky'});
% X = getData(MRSIStruct);                    % [Nt_total  nKy  Nextra]
% if ndims(X)=3 || size(X,1)=Nt_total || size(X,2)~=nKy
%     error('Unexpected shape after reshapeDimensions: %s', mat2str(size(X)));
% end
% Nextra = size(X,3);
% 
% vprint('Nk=%d, nKy=%d, kPtsPerCycle=%d, NPtemporal=%d, Nx=%d, Ny=%d, Nextra=%d\n', ...
%        Nk, nKy, kPtsPerCycle, NPtemporal, Nx, Ny, Nextra);
% 
% if Nk ~= kPtsPerCycle * nKy
%     error('Trajectory rows (%d) != kPtsPerCycle (%d) * nKy (%d).', Nk, kPtsPerCycle, nKy);
% end
% 
% % -------- voxel pitch & NUFFT grid shift from coordinates --------
% dx = median(abs(diff(xCoords(:))));
% dy = median(abs(diff(yCoords(:))));
% if ~(isfinite(dx)&&dx>0 && isfinite(dy)&&dy>0)
%     error('Bad voxel pitch: dx=%.3g, dy=%.3g', dx, dy);
% end
% 
% % Find exact n_shift so that (n - n_shift)*d? ~ coords (least squares)
% x_idx = (0:Nx-1).';  x_shift = median( x_idx - xCoords(:)/dx );
% y_idx = (0:Ny-1).';  y_shift = median( y_idx - yCoords(:)/dy );
% n_shift = [y_shift, x_shift];     % [Ny-shift, Nx-shift]
% Nd      = [Ny, Nx];
% vprint('Inferred n_shift = [%.3f, %.3f] (ideal would be [%g, %g])\n', ...
%        y_shift, x_shift, Ny/2, Nx/2);
% 
% % -------- SFT operator (for calibration only; uses physical coords) --------
% [xg, yg]  = meshgrid(xCoords, yCoords);     % Ny x Nx
% imgTrajXY = [xg(:), yg(:)];                 % [Ny*Nx x 2]
% sftOp = sft2_Operator_local(kXY, imgTrajXY, 1);  % [Ny*Nx x Nk], normalized by Nk
% 
% % -------- NUFFT candidates (axis/sign) --------
% % NUFFT expects om (radians per index), so om = 2π * k(1/mm) * Δ(mm)
% om_xy = [ 2*pi*kXY(:,1)*dx , 2*pi*kXY(:,2)*dy ];  % [ωy? ωx? careful]
% om_yx = [ 2*pi*kXY(:,2)*dy , 2*pi*kXY(:,1)*dx ];
% cands = {};
% cands = pushcand(cands, 'om=[Ky,Kx] -> [y,x]', om_yx, Nd);
% cands = pushcand(cands, 'om=[Kx,Ky] -> [y,x]', om_xy, Nd);
% % sign flips
% cands = pushcand(cands, '-Ky,+Kx', [-om_yx(:,1), +om_yx(:,2)], Nd);
% cands = pushcand(cands, '+Ky,-Kx', [+om_yx(:,1), -om_yx(:,2)], Nd);
% cands = pushcand(cands, '-Kx,+Ky', [-om_xy(:,1), +om_xy(:,2)], Nd);
% cands = pushcand(cands, '+Kx,-Ky', [+om_xy(:,1), -om_xy(:,2)], Nd);
% cands = pushcand(cands, '-Ky,-Kx', [-om_yx(:,1), -om_yx(:,2)], Nd);
% cands = pushcand(cands, '-Kx,-Ky', [-om_xy(:,1), -om_xy(:,2)], Nd);
% 
% % -------- quick calibration (pick axis/sign and complex scale α) --------
% cal_idx = 1:NPtemporal;
% if opt.Calib
%     nCal = min(opt.nCalT, NPtemporal);
%     cal_idx = unique(max(1, round(linspace(1, NPtemporal, nCal))));
% end
% 
% best.st = [];
% best.alpha = 1;
% best.resid = inf;
% best.name = '';
% 
% Jd = [6,6];               % kernel size
% Kd = 2*Nd;                % oversamp
% 
% for ic = 1:size(cands,1)
%     name = cands{ic,1};
%     om   = cands{ic,2};
%     Nd_  = cands{ic,3};
% 
%     st = nufft_init(om, Nd_, Jd, Kd, n_shift);
% 
%     if isempty(cal_idx)
%         alpha = 1; resid = 0;
%     else
%         num = 0+0i; den = 0; rsum = 0;
%         for it = cal_idx
%             i0 = (it-1)*kPtsPerCycle + 1;
%             i1 = it*kPtsPerCycle;
% 
%             Y = double(reshape(X(i0:i1, :, :), [], Nextra));   % [Nk x Nextra]
%             Zsft = sftOp * Y;                                  % [Ny*Nx x Nextra]
%             Znu  = nufft_adj(Y, st);                           % [Ny*Nx x Nextra]
% 
%             % LS alpha: minimize ||alpha*Znu - Zsft||
%             a = Znu(:); b = Zsft(:);
%             num = num + (a' * b);
%             den = den + (a' * a);
%         end
%         alpha = (den~=0) * (num/den) + (den==0) * 1;
% 
%         for it = cal_idx
%             i0 = (it-1)*kPtsPerCycle + 1;
%             i1 = it*kPtsPerCycle;
%             Y   = double(reshape(X(i0:i1, :, :), [], Nextra));
%             Zsft= sftOp * Y;
%             Znu = alpha * nufft_adj(Y, st);
%             rsum= rsum + norm(Znu(:) - Zsft(:))^2;
%         end
%         resid = rsum;
%     end
% 
%     vprint('Cand %d: %s  | alpha=%+.6e%+.6ei  | resid=%.3e\n', ic, name, real(alpha), imag(alpha), resid);
% 
%     if resid < best.resid
%         best.st    = st;
%         best.alpha = alpha;
%         best.resid = resid;
%         best.name  = name;
%     end
% end
% 
% vprint('Selected: %s  | alpha=%+.6e%+.6ei  | resid=%.3e\n', best.name, real(best.alpha), imag(best.alpha), best.resid);
% 
% % -------- full reconstruction --------
% image = zeros(NPtemporal, Ny, Nx, Nextra, 'like', double(X));
% 
% for it = 1:NPtemporal
%     i0 = (it-1)*kPtsPerCycle + 1;
%     i1 = it*kPtsPerCycle;
% 
%     Y = double(reshape(X(i0:i1, :, :), [], Nextra));      % [Nk x Nextra]
%     Z = best.alpha * nufft_adj(Y, best.st);               % [Ny*Nx x Nextra]
%     Z = reshape(Z, [Ny, Nx, Nextra]);
%     image(it, :, :, :) = Z;
% end
% 
% % -------- restore dims like slow path --------
% MRSIStruct = setData(MRSIStruct, double(image));
% 
% kyDim = getDimension(MRSIStruct, 'ky');
% prevPermute = removeDimPrevPermute(prevPermute, kyDim);
% prevPermute = addDimPrevPermute(prevPermute, 'y', kyDim);
% prevPermute = addDimPrevPermute(prevPermute, 'x', kyDim + 1);
% 
% prevSize(1) = NPtemporal;
% prevSize(2) = Ny;
% prevSize    = [prevSize(1:2), Nx, prevSize(3:end)];
% 
% MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
% 
% % spectral values identical to slow path
% MRSIStruct = calculateSpectralValues_local(MRSIStruct, kPtsPerCycle, NPtemporal);
% 
% % flag
% MRSIStruct = setFlags(MRSIStruct, 'spatialFT', true);
% 
% vprint('=== END: op_NUFFTSpatial ===\n');
% end
% 
% % ================= helpers =================
% function v = ifvprintf(flag, varargin)
% if flag, fprintf(varargin{:}); end
% v = [];
% end
% 
% function [kTable, kArray] = readKFile_local(kFileName)
% kTable = [];
% kArray = [];
% if isempty(kFileName) || ~isfile(kFileName), return; end
% try
%     T = readtable(kFileName);
%     kTable = T;
%     vn = lower(string(T.Properties.VariableNames));
%     if any(vn=="kx") && any(vn=="ky")
%         kArray = [T{:, find(vn=="kx",1)}, T{:, find(vn=="ky",1)}];
%     elseif width(T) >= 3
%         kArray = [T{:,2}, T{:,3}];
%     end
% catch
%     try
%         A = readmatrix(kFileName);
%         if size(A,2) >= 3
%             kArray = A(:,2:3);
%             kTable = A;
%         end
%     catch
%         kArray = [];
%         kTable = [];
%     end
% end
% end
% 
% function kPtsPerCycle = getKPtsPerCycle_local(kTable)
% kPtsPerCycle = 126;
% try
%     if istable(kTable)
%         vn = lower(string(kTable.Properties.VariableNames));
%         if any(vn=="tr")
%             nTR = max(kTable{:, find(vn=="tr",1)}, [], 'all');
%             if nTR > 0
%                 kPtsPerCycle = round(height(kTable)/nTR);
%             end
%         end
%     elseif isnumeric(kTable) && size(kTable,2) >= 5
%         nTR = max(kTable(:,5), [], 'all');
%         if nTR > 0, kPtsPerCycle = round(size(kTable,1)/nTR); end
%     end
% catch
%     kPtsPerCycle = 126;
% end
% end
% 
% function S = sft2_Operator_local(InTraj, OutTraj, Ift_flag)
% % InTraj: [Nk x 2] = [Kx, Ky] in 1/mm
% % OutTraj: [Ny*Nx x 2] = [x(mm), y(mm)]
% % Ift_flag=1 => adjoint (divide by Nk)
% if Ift_flag, Expy = 2*pi*1i; else, Expy = -2*pi*1i; end
% NOut = size(OutTraj,1);
% Nk   = size(InTraj,1);
% S = zeros(NOut, Nk);
% % exp(i 2π (x*Kx + y*Ky))
% for j = 1:NOut
%     S(j,:) = exp( Expy*( OutTraj(j,1)*InTraj(:,1).' + OutTraj(j,2)*InTraj(:,2).' ) );
% end
% if Ift_flag, S = S / Nk; end
% end
% 
% function MRSIStruct = calculateSpectralValues_local(MRSIStruct, kPtsPerCycle, NPtemporal)
% adcDwellTime       = getAdcDwellTime(MRSIStruct);
% spectralDwellTime  = kPtsPerCycle * adcDwellTime;
% spectralWidth      = 1 / spectralDwellTime;
% spectralTime       = 0:spectralDwellTime:spectralDwellTime*(NPtemporal-1);
% MRSIStruct = setSpectralWidth(MRSIStruct, spectralWidth);
% MRSIStruct = setSpectralDwellTime(MRSIStruct, spectralDwellTime);
% MRSIStruct = setSpectralTime(MRSIStruct, spectralTime);
% end
% 
% function c = pushcand(c, name, om, Nd)
% % append one candidate row
% c(end+1,:) = {name, om, Nd};
% end