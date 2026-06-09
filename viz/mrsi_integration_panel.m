function mrsi_integration_panel(t1Path, mrsi4DPath, brainmask)
%MRSI_INTEGRATION_PANEL  Brain-masked ppm-integration heatmap + click-to-spectrum.
%
%   mrsi_integration_panel(t1Path, mrsi4DPath, brainmask)
%
%   - t1Path:     T1 NIfTI (background for nii_viewer)
%   - mrsi4DPath: 4-D MRSI NIfTI written by mrsi_on_t1_map.m.
%                 Its companion <base>_ftSpec.mat is loaded for the spectra.
%   - brainmask:  2-D logical [Ny x Nx] = ccav_w.mask.brainmasks.
%
%   What the function does, in order:
%     1. Load ftSpec from the .mat companion of mrsi4DPath
%     2. Build a control window (PPM min/max, preset buttons, heatmap, spectrum)
%     3. On Update / preset:  integrate |Re(spec)| over [lo, hi] inside the
%        brain mask, normalise by the 99th percentile, redraw the heatmap,
%        save it as a 3-D NIfTI, and overlay it on the T1 in nii_viewer with
%        the "hot" colormap.
%     4. Click a heatmap voxel -> plot its spectrum with the integration
%        window highlighted in green.

arguments
    t1Path     (1,:) char {mustBeFile}
    mrsi4DPath (1,:) char {mustBeFile}
    brainmask         {mustBeNumericOrLogical}
end

%% ----- 1.  load ftSpec ---------------------------------------------------
[outDir, base] = fileparts(mrsi4DPath);
if endsWith(base, '.nii'), base = extractBefore(base, '.nii'); end
matFile = fullfile(outDir, [base '_ftSpec.mat']);
S = load(matFile);
fn = fieldnames(S);
ftSpec = S.(fn{1});                          % first struct in the .mat

ppm  = ftSpec.ppm(:);
spec = real(permute(ftSpec.data, ...
              [ftSpec.dims.f, ftSpec.dims.y, ftSpec.dims.x]));   % (Nf, Ny, Nx)
[~, Ny, Nx] = size(spec);

bm = logical(brainmask);
assert(isequal(size(bm), [Ny, Nx]), 'brainmask must be %dx%d', Ny, Nx);
fprintf('Brain mask: %d / %d voxels (%.1f%%).\n', sum(bm(:)), numel(bm), 100*mean(bm(:)));

heatPath = fullfile(outDir, [base '_intHeatmap.nii.gz']);

%% ----- 2.  state shared across callbacks --------------------------------
mp    = zeros(Ny, Nx);                  % current normalised heatmap
sel   = [round(Ny/2), round(Nx/2)];     % selected voxel  [iy, ix]
lo    = min(ppm);   hi  = max(ppm);     % current PPM window
n_pts = numel(ppm);                     % # spectral points used

%% ----- 3.  build the UI -------------------------------------------------
ctrl = figure('Name','MRSI Integration Panel', 'NumberTitle','off', ...
              'MenuBar','none', 'ToolBar','none', ...
              'Position',[100 100 920 540], 'Color',[0.94 0.94 0.94]);

label('PPM min:', [20 490 60 22]);
e_lo = uicontrol(ctrl, 'Style','edit', 'String',num2str(lo,'%.2f'), ...
                 'Position',[85 490 70 25]);
label('PPM max:', [165 490 60 22]);
e_hi = uicontrol(ctrl, 'Style','edit', 'String',num2str(hi,'%.2f'), ...
                 'Position',[230 490 70 25]);

% metabolite preset buttons
presets = {'NAA',1.95,2.10; 'Cr',2.95,3.05; 'Cho',3.15,3.30; ...
           'Lac',1.25,1.40; 'mI',3.50,3.65; 'Lip',0.90,1.45;
           'Full', min(ppm), max(ppm)};
for k = 1:size(presets,1)
    uicontrol(ctrl, 'Style','pushbutton', 'String',presets{k,1}, ...
              'Position',[310 + (k-1)*62, 490, 58, 25], ...
              'Callback', @(~,~) preset(presets{k,2}, presets{k,3}));
end
uicontrol(ctrl, 'Style','pushbutton', 'String','Update', ...
          'FontWeight','bold', 'BackgroundColor',[0.65 0.85 0.65], 'FontSize',11, ...
          'Position',[800 490 100 30], 'Callback', @update);

% info banner: shows the actual ppm window + how many points were summed
info = uicontrol(ctrl, 'Style','text', 'String','', ...
                 'Position',[20 455 880 24], 'HorizontalAlignment','left', ...
                 'BackgroundColor',[0.94 0.94 0.94], ...
                 'FontWeight','bold', 'FontSize',11, ...
                 'ForegroundColor',[0.1 0.35 0.6]);

% heatmap axes  (black background + hot colormap)
heat = axes('Parent',ctrl, 'Units','pixels', 'Position',[50 50 380 380], ...
            'Color','k', 'XColor',[0.6 0.6 0.6], 'YColor',[0.6 0.6 0.6]);
% spectrum axes
spec_ax = axes('Parent',ctrl, 'Units','pixels', 'Position',[490 50 400 380]);

compute();      % do the initial integration over the full ppm range


%% ----- 4.  callbacks ----------------------------------------------------
    function preset(a, b)
        set(e_lo, 'String', num2str(a,'%.2f'));
        set(e_hi, 'String', num2str(b,'%.2f'));
        update();
    end

    function update(~,~)
        a = str2double(get(e_lo,'String'));
        b = str2double(get(e_hi,'String'));
        if ~(isfinite(a) && isfinite(b) && a < b)
            set(info, 'String','Invalid PPM range', 'ForegroundColor',[0.8 0 0]);
            return
        end
        lo = a;  hi = b;
        compute();
    end

    function clickHeat(~,~)
        cp = get(heat, 'CurrentPoint');
        % Display = fliplr(rot90(mp, 2)).  Undo so the click maps back to
        % the original (iy, ix) voxel:
        %   display (disp_x, disp_y)  ->  (ix = disp_x, iy = Ny - disp_y + 1)
        disp_x = round(cp(1,1));   disp_y = round(cp(1,2));
        ix = disp_x;
        iy = Ny - disp_y + 1;
        if iy>=1 && iy<=Ny && ix>=1 && ix<=Nx
            sel = [iy, ix];  drawHeat();  drawSpec();
        end
    end


%% ----- 5.  integrate + normalise + render -----------------------------
    function compute()
        % integrate Re(spec) over [lo, hi]
        m     = (ppm >= lo) & (ppm <= hi);
        n_pts = sum(m);
        if n_pts > 0
            raw = abs(squeeze(sum(spec(m, :, :), 1)));   % (Ny, Nx)
        else
            raw = zeros(Ny, Nx);
        end

        % brain-only, normalise by 99th percentile inside brain
        raw(~bm) = 0;
        if any(raw(bm) > 0)
            raw = raw / max(prctile(raw(bm), 99), eps);
        end
        mp = min(raw, 1);

        % info banner
        if n_pts == 0
            set(info, 'ForegroundColor',[0.8 0 0], ...
                'String', sprintf('Integrating [%.2f, %.2f] ppm  -  NO points in window', lo, hi));
        else
            p = ppm(m);
            set(info, 'ForegroundColor',[0.1 0.35 0.6], ...
                'String', sprintf( ...
                  'Integrating Re(ftSpec) over [%.2f, %.2f] ppm  |  %d points summed  |  spans %.3f to %.3f ppm', ...
                  lo, hi, n_pts, min(p), max(p)));
        end

        saveNifti();
        drawHeat();
        drawSpec();
        refreshViewer();
    end

    function drawHeat()
        cla(heat);
        % Heatmap display transforms (display-only -- the underlying mp /
        % bm arrays and the NIfTI on disk are NOT touched):
        %   1) rot90(.,2)  -- 180 deg rotation
        %   2) fliplr(.)   -- horizontal mirror (left <-> right)
        % The two together are equivalent to a single flipud(.).
        mp_disp = fliplr(rot90(mp, 2));
        bm_disp = fliplr(rot90(bm, 2));
        h = imagesc(heat, mp_disp, 'AlphaData', double(bm_disp));
        set(h, 'ButtonDownFcn', @clickHeat);
        set(heat, 'Color','k', 'CLim',[0 1], 'YDir','normal');
        axis(heat, 'image');  axis(heat, 'on');
        colormap(heat, hot(256));
        c = colorbar(heat, 'Color',[0.7 0.7 0.7]);
        c.Label.String = 'normalised |integral|';
        title(heat, sprintf('Heatmap  [%.2f, %.2f] ppm  (%d pts)', lo, hi, n_pts), ...
              'Color',[0.2 0.2 0.2]);
        hold(heat,'on');
        % Cyan marker follows the same composed transform:
        %   voxel (iy, ix) ends up at display coords (ix, Ny - iy + 1)
        plot(heat, sel(2), Ny - sel(1) + 1, 'co', ...
             'MarkerSize',10, 'LineWidth',1.5, ...
             'HitTest','off', 'PickableParts','none');
        hold(heat,'off');
    end

    function drawSpec()
        cla(spec_ax);
        iy = sel(1);  ix = sel(2);
        if ~bm(iy, ix)
            text(spec_ax, 0.5, 0.5, sprintf('(x=%d, y=%d) outside brain mask', ix, iy), ...
                'Units','normalized', 'HorizontalAlignment','center', 'Color',[0.55 0.55 0.55]);
            set(spec_ax,'XTick',[],'YTick',[]);  title(spec_ax,'');
            return
        end
        plot(spec_ax, ppm, spec(:, iy, ix), 'b-', 'LineWidth',1.2);
        set(spec_ax, 'XDir','reverse');  grid(spec_ax,'on');
        xlim(spec_ax, [min(ppm) max(ppm)]);
        xlabel(spec_ax, 'Chemical shift (ppm)');  ylabel(spec_ax, 'Re(ftSpec)');
        title(spec_ax, sprintf('Spectrum @ voxel (x=%d, y=%d)  -  green band = %d summed points', ...
              ix, iy, n_pts));
        if n_pts > 0
            yl = ylim(spec_ax);
            hold(spec_ax,'on');
            patch(spec_ax, [lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.2 0.7 0.2], 'FaceAlpha',0.18, 'EdgeColor',[0.1 0.5 0.1]);
            hold(spec_ax,'off');
        end
    end


%% ----- 6.  save NIfTI + open viewer ------------------------------------
    function saveNifti()
        % Match mrsi_on_t1_map's orientation: transpose then flip dim 1
        img = flip(mp.', 1);
        nii = nii_tool('load', mrsi4DPath);
        nii.img            = single(img);
        nii.hdr.dim(1)     = int16(3);
        nii.hdr.dim(2:4)   = int16([size(img,1), size(img,2), 1]);
        nii.hdr.dim(5)     = int16(1);
        nii.hdr.pixdim(5)  = single(0);
        nii.hdr.datatype   = int16(16);
        nii.hdr.bitpix     = int16(32);
        nii.hdr.cal_min    = single(0.02);
        nii.hdr.cal_max    = single(1.0);
        nii.hdr.descrip    = sprintf('|Re ftSpec| [%.2f, %.2f] ppm', lo, hi);
        nii_tool('save', nii, heatPath);
    end

    function refreshViewer()
        fh = findobj('Type','figure', '-regexp','Name','^nii_viewer');
        if ~isempty(fh), try, close(fh); catch, end, end
        nii_viewer(t1Path, heatPath);
        drawnow;
        styleOverlay();
    end

    function styleOverlay()
        % After nii_viewer opens, switch the heatmap overlay to the "hot"
        % colormap at full alpha (matches the in-panel heatmap).
        fh = findobj('Type','figure', '-regexp','Name','^nii_viewer');
        if isempty(fh), return; end
        hs = guidata(fh(1));
        if ~isstruct(hs), return; end

        % nii_viewer's LUT list (see nii_viewer.m line 355):
        %   1 gray  2 red  ...  15 HOT  16 cool ...
        names = get(hs.files,'String');
        idx = find(iscell(names) & contains(string(names), '_intHeatmap'), 1);
        if ~isempty(idx)
            set(hs.files, 'Value', idx);  fire(hs.files);
        end
        set(hs.lut, 'Value', 15);            fire(hs.lut);
        try, hs.alpha.setValue(1.0);  fire(hs.alpha);  catch, end
        try, hs.lb.setValue(0.02);    fire(hs.lb);     catch, end
        try, hs.ub.setValue(1.0);     fire(hs.ub);     catch, end
    end


%% ----- 7.  tiny helpers ------------------------------------------------
    function label(str, pos)
        uicontrol(ctrl, 'Style','text', 'String',str, 'Position',pos, ...
            'BackgroundColor',[0.94 0.94 0.94], 'HorizontalAlignment','right');
    end
end


%% =========================================================================
%   top-level helpers
%% =========================================================================
function fire(h)
    cb = get(h, 'Callback');
    if isempty(cb), return; end
    try
        if iscell(cb),                  feval(cb{1}, h, [], cb{2:end});
        elseif isa(cb,'function_handle'), cb(h, []);
        end
    catch
    end
end

function mustBeNumericOrLogical(x)
    if ~(isnumeric(x) || islogical(x)) || isempty(x)
        error('brainmask must be a non-empty numeric/logical array.');
    end
end
