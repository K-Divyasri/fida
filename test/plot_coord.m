function plot_coord(coordFile)
% plot_coord  Plot an LCModel .coord: data, fit, residual, baseline.
%   Residual = data - fit. FLAT random noise = clean fit.
%   Structured oscillation under peaks = truncation / lineshape mismatch.
%
%   plot_coord('F:\...\20x20_ftSpec_smooth_lcm.coord')

    txt   = fileread(coordFile);
    lines = regexp(txt, '\r?\n', 'split');

    ppm = grab(lines, 'points on ppm-axis');
    dat = grab(lines, 'phased data points');
    fit = grab(lines, 'points of the fit');
    bg  = grab(lines, 'background');           % baseline (may be empty)

    n = min([numel(ppm) numel(dat) numel(fit)]);
    assert(n > 0, 'Could not parse data/fit/ppm blocks from %s', coordFile);
    ppm = ppm(1:n); dat = dat(1:n); fit = fit(1:n);
    res = dat - fit;
    if numel(bg) >= n, bg = bg(1:n); else, bg = []; end

    figure('Color','w','Name',coordFile);
    subplot(3,1,1);
        plot(ppm, dat, 'k', 'LineWidth', 1); hold on;
        plot(ppm, fit, 'r', 'LineWidth', 1);
        if ~isempty(bg), plot(ppm, bg, 'Color',[.6 .6 .6]); end
        set(gca,'XDir','reverse'); grid on; xlim([min(ppm) max(ppm)]);
        legend('data','fit','baseline','Location','best');
        title('data vs fit'); ylabel('a.u.');
    subplot(3,1,2);
        plot(ppm, res, 'b', 'LineWidth', 1);
        set(gca,'XDir','reverse'); grid on; xlim([min(ppm) max(ppm)]);
        title('RESIDUAL = data - fit   (flat noise = clean, wiggles = structure)');
        ylabel('a.u.');
    subplot(3,1,3);
        plot(ppm, res, 'b'); set(gca,'XDir','reverse'); grid on;
        xlim([1.0 3.6]);                       % zoom NAA/Cr/Cho/Glu region
        title('residual zoom 1.0-3.6 ppm'); xlabel('ppm'); ylabel('a.u.');

    fprintf('residual std / max|data| = %.4f\n', std(res)/max(abs(dat)));
end

function v = grab(lines, key)
% Return the numeric block that FOLLOWS the header line containing `key`,
% reading floats until the next non-numeric header line.
    v = [];
    started = false;
    for i = 1:numel(lines)
        L = lines{i};
        if ~started
            if ~isempty(regexpi(L, key, 'once')), started = true; end
            continue;
        end
        nums = sscanf(L, '%f');
        if isempty(nums)
            if ~isempty(strtrim(L)), break; end   % hit next header -> stop
        else
            v = [v; nums]; %#ok<AGROW>
        end
    end
end
