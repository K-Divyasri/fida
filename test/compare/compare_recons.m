%% compare_recons.m — diff Jason's recon vs the user's, stage by stage.
% Loads the saved intermediates (no toolbox functions needed). Run AFTER both
% rosette4divya_save.m (Jason) and run_mine_save.m (user) have been run.
% Handles the two dims conventions: Jason labels the spectral axis dims.t (no .f);
% the user uses dims.f (dims.t=0).

jDir = 'F:\fida\divya\20260831_ulinePhantom_test\subject04\compare\jason';
mDir = 'F:\fida\divya\20260831_ulinePhantom_test\subject04\compare\mine';

fprintf('\n=================== RECON COMPARISON ===================\n');

%% --- 1. RECON OUTPUT (ccav) : the decisive spatial image ---------------
J = loadsnap(jDir,'ccav');   M = loadsnap(mDir,'ccav');
Ji = img_t1(J);  Mi = img_t1(M);
showpair('ccav (recon output) |img(t=1)|', Ji, Mi);
centerReport('ccav metab', Ji, Mi);

Jw = loadsnap(jDir,'ccav_w'); Mw = loadsnap(mDir,'ccav_w');
Jwi = img_t1(Jw); Mwi = img_t1(Mw);
showpair('ccav_w (water ref) |img(t=1)|', Jwi, Mwi);
centerReport('ccav_w water', Jwi, Mwi);

%% --- 2. FINAL SPECTRUM at center + a good (edge) voxel -----------------
Js = loadsnap(jDir,'ftSpec_smooth'); Ms = loadsnap(mDir,'ftSpec_smooth');
[ny,nx] = size(Mi);
cy = round(ny/2); cx = round(nx/2);          % center voxel
ey = round(ny/2); ex = max(round(nx/2)-8,1); % an off-center voxel

specOverlay(Js, Ms, cy, cx, 'CENTER voxel');
specOverlay(Js, Ms, ey, ex, 'off-center voxel');

fprintf('\nRead: center/median >1 = filled (no dropout). If Jason and Mine\n');
fprintf('differ in the center ratio, that is the exact-DFT vs NUFFT + DCF\n');
fprintf('difference. Transposed/flipped images = display convention, not a bug.\n');

%% ---- locals ----------------------------------------------------------
function S = loadsnap(d,name)
    L = load(fullfile(d,[name '.mat']));  S = L.snapStruct;
end

function v = fld(dimsStruct, name)
    if isfield(dimsStruct,name), v = dimsStruct.(name); else, v = 0; end
end

function im = img_t1(S)
    d = S.data; td = fld(S.dims,'t'); if td==0, td = 1; end
    nd = ndims(d); d = permute(d, [td, setdiff(1:nd, td)]);
    sub = repmat({':'},1,ndims(d)); sub{1} = 1; d = d(sub{:});
    im = squeeze(abs(d));
    while ndims(im) > 2, im = squeeze(sum(im, ndims(im))); end
end

function [sp,ppm] = spec_at(S,iy,ix)
    d = S.data; ppm = S.ppm(:);
    yd = fld(S.dims,'y'); xd = fld(S.dims,'x');
    % spectral axis: prefer .f, else .t, else the dim whose size == numel(ppm),
    % else whatever remains after y,x.
    fd = fld(S.dims,'f');
    if fd==0, fd = fld(S.dims,'t'); end
    if fd==0 && ~isempty(ppm)
        cand = find(size(d)==numel(ppm),1); if ~isempty(cand), fd = cand; end
    end
    if fd==0
        rem = setdiff(1:ndims(d), [yd xd]); fd = rem(1);
    end
    if yd==0 || xd==0
        rem = setdiff(1:ndims(d), fd);
        if yd==0, yd = rem(1); end
        if xd==0, xd = rem(min(2,numel(rem))); end
    end
    D = real(permute(d, [fd yd xd]));
    iy = min(iy,size(D,2)); ix = min(ix,size(D,3));
    sp = squeeze(D(:,iy,ix));
    if isempty(ppm) || numel(ppm)~=numel(sp), ppm = (1:numel(sp))'; end
end

function showpair(ttl, A, B)
    figure('Name',ttl,'Color','w');
    subplot(1,3,1); imagesc(A); axis image; colorbar; title('Jason');
    subplot(1,3,2); imagesc(B); axis image; colorbar; title('Mine');
    subplot(1,3,3);
    if isequal(size(A),size(B))
        imagesc(abs(A./max(A(:)) - B./max(B(:)))); title('|norm diff|');
    else
        title(sprintf('size mismatch J%s M%s', mat2str(size(A)), mat2str(size(B))));
    end
    axis image; colorbar; try, sgtitle(ttl,'Interpreter','none'); catch, end
end

function centerReport(tag, A, B)
    ca = A(round(end/2),round(end/2));  cb = B(round(end/2),round(end/2));
    fprintf('%-14s center/median :  Jason=%.2f   Mine=%.2f   (dropout if <<1)\n', ...
        tag, ca/median(A(:)), cb/median(B(:)));
end

function specOverlay(Js, Ms, iy, ix, tag)
    [jsp,jppm] = spec_at(Js,iy,ix);  [msp,mppm] = spec_at(Ms,iy,ix);
    figure('Name',[tag ' spectrum'],'Color','w');
    plot(jppm, jsp./max(abs(jsp)), 'b', mppm, msp./max(abs(msp)), 'r', 'LineWidth',1);
    set(gca,'XDir','reverse'); grid on; legend('Jason','Mine');
    xlabel('ppm'); title(sprintf('%s (%d,%d) — normalized', tag, iy, ix));
end
