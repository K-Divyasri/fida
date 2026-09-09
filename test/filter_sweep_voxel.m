%% filter_sweep_voxel.m
% Look at ONE voxel and sweep the time-domain line-broadening filter, plotting
% the spectrum after each Hz value so you can pick the lineshape you want.
%
% PREREQ: run run_uline_full FIRST with lb = 0 (so the sweep filter is the
% ABSOLUTE broadening, not added on top of a pipeline filter). ftSpec_smooth
% must be in the workspace.
%
% Usage:  filter_sweep_voxel
% ---------------------------------------------------------------------------

vX = 14;  vY = 18;                       % op_CSItoMRS(struct, x, y) -> your (14,18)
lbs = 1 : 0.25 : 5;                      % 1, 1.25, ... 5 Hz (exponential LB)

assert(exist('ftSpec_smooth','var')==1, ...
    'ftSpec_smooth not found. Run run_uline_full first (recommend lb=0).');

vox0 = op_CSItoMRS(ftSpec_smooth, vX, vY);   % single-voxel FID-A struct
n = numel(lbs);

% --- (1) tiled: one panel per Hz ---
figure('Color','w','Name',sprintf('voxel (%d,%d): filter sweep (panels)', vX, vY));
nc = ceil(sqrt(n));  nr = ceil(n/nc);
for i = 1:n
    voxf = op_filter(vox0, lbs(i));          % exponential LB, recomputes specs
    subplot(nr, nc, i);
    plot(voxf.ppm, real(voxf.specs), 'LineWidth', 1);
    set(gca,'XDir','reverse'); grid on; xlim([0 4.5]);
    title(sprintf('%.2f Hz', lbs(i)));
end
sgtitle(sprintf('voxel (%d,%d) — real spectrum vs line broadening', vX, vY));

% --- (2) overlay: all Hz on one axes ---
figure('Color','w','Name',sprintf('voxel (%d,%d): filter sweep (overlay)', vX, vY));
hold on; leg = cell(n,1);  cmap = parula(n);
for i = 1:n
    voxf = op_filter(vox0, lbs(i));
    plot(voxf.ppm, real(voxf.specs), 'Color', cmap(i,:), 'LineWidth', 1);
    leg{i} = sprintf('%.2f Hz', lbs(i));
end
set(gca,'XDir','reverse'); grid on; xlim([0 4.5]);
xlabel('ppm'); ylabel('real'); colormap(parula(n));
cb = colorbar; cb.Ticks = []; cb.Label.String = sprintf('%.2f -> %.2f Hz', lbs(1), lbs(end));
title(sprintf('voxel (%d,%d) — narrow (sharp, noisy) -> broad (smooth, lower peak)', vX, vY));

fprintf('Swept %d filters on voxel (%d,%d): %.2f .. %.2f Hz\n', n, vX, vY, lbs(1), lbs(end));
fprintf('Pick the smallest LB that still gives clean, non-noisy peaks (more LB = more SNR\n');
fprintf('but broader lines + more overlap), then set lb to that in run_uline_full.\n');
