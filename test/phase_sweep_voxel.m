%% phase_sweep_voxel.m
% Look at ONE voxel: show the CURRENT spectrum + the spectrum after each
% first-order phase from 1.2e-4 to 16.2e-4 s, so you can pick the upright one.
%
% PREREQ: run run_uline_full FIRST with foPhase = 0 (so the sweep phase is the
% ABSOLUTE phase, not stacked on a pipeline phase). ftSpec_smooth in workspace.
%
% Usage:  phase_sweep_voxel
% ---------------------------------------------------------------------------

vX = 14;  vY = 16;                       % voxel in the bottle: op_CSItoMRS(struct,x,y)
phases = (1.2 : 1.0 : 16.2) * 1e-4;      % 1.2e-4 ... 16.2e-4 s

assert(exist('ftSpec_smooth','var')==1, ...
    'ftSpec_smooth not found. Run run_uline_full first (recommend foPhase=0).');

vox0 = op_CSItoMRS(ftSpec_smooth, vX, vY);   % single-voxel FID-A struct (CURRENT)
n = numel(phases);

% --- (1) tiled: panel 1 = current, then one per phase ---
figure('Color','w','Name',sprintf('voxel (%d,%d): current + phase sweep', vX, vY));
nc = ceil(sqrt(n+1));  nr = ceil((n+1)/nc);
subplot(nr, nc, 1);
plot(vox0.ppm, real(vox0.specs), 'k', 'LineWidth', 1.2);
set(gca,'XDir','reverse'); grid on; xlim([0 4.5]);
title('current (ph1 = 0)');
for i = 1:n
    voxp = op_addphase(vox0, 0, phases(i), 4.65, true);   % suppressPlot = true
    subplot(nr, nc, i+1);
    plot(voxp.ppm, real(voxp.specs), 'LineWidth', 1);
    set(gca,'XDir','reverse'); grid on; xlim([0 4.5]);
    title(sprintf('ph1 = %.1e s', phases(i)));
end
sgtitle(sprintf('voxel (%d,%d) — current + real spectrum vs first-order phase', vX, vY));

% --- (2) overlay: current (thick black) + all phases ---
figure('Color','w','Name',sprintf('voxel (%d,%d): overlay', vX, vY));
hold on; cmap = lines(n);  leg = [{'current'}, cell(1,n)];
plot(vox0.ppm, real(vox0.specs), 'k', 'LineWidth', 2);
for i = 1:n
    voxp = op_addphase(vox0, 0, phases(i), 4.65, true);
    plot(voxp.ppm, real(voxp.specs), 'Color', cmap(i,:), 'LineWidth', 1);
    leg{i+1} = sprintf('%.1e', phases(i));
end
set(gca,'XDir','reverse'); grid on; xlim([0 4.5]);
xlabel('ppm'); ylabel('real'); legend(leg, 'Location','best');
title(sprintf('voxel (%d,%d) — black = current; pick the phase with upright peaks', vX, vY));

fprintf('Voxel (%d,%d): current + %d phases swept %.1e .. %.1e s\n', vX, vY, n, phases(1), phases(end));
