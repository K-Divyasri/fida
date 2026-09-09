%APPLY_WATER_SCALING_ALL  Rescale the workspace `map` (LCModel raw) to mM, for the
%   3-compartment ROSETTE phantom (Small Inner / Large Inner / Outer), metabolites
%   Acetate (Act) / Creatine (Cr) / Choline (Cho) / Lactate (Lac).
%
%   Assumes map, crlb, LW, SNR are in the workspace. Rescales `map` in place;
%   crlb/LW/SNR untouched. Run ONCE (re-running double-scales).
%
%   LCModel already water-scaled each metabolite (Smet/Swater * WCONC). This
%   OPTIONAL step divides out the residual T1/T2 relaxation difference between
%   water and metabolites. Only trustworthy with MEASURED relaxation; OFF by
%   default. Ground-truth concentrations for each compartment are in GT below,
%   for per-ROI validation (the scaling itself is uniform per voxel).

%% ---- GROUND TRUTH (mM) — from final_phantom sheet (actual achieved conc) ----
% Molar concentrations (mol/L * 1000). Use these to validate per-compartment ROIs.
%   Small Inner  70 mL   |  Large Inner 650 mL  |  Outer 3800 mL
GT.SmallInner = struct('Act',40.42, 'Cr',32.60, 'Cho',10.08, 'Lac', 4.09);
GT.LargeInner = struct('Act', 6.09, 'Cr',16.03, 'Cho',14.48, 'Lac', 8.01);
GT.Outer      = struct('Act',20.01, 'Cr',16.00, 'Cho', 4.81, 'Lac', 0.00);
assignin('base','GT',GT);
fprintf('\nGROUND TRUTH (mM):\n');
fprintf('  %-11s %6s %6s %6s %6s\n','compartment','Act','Cr','Cho','Lac');
cn = fieldnames(GT);
for i=1:numel(cn)
    g = GT.(cn{i});
    fprintf('  %-11s %6.2f %6.2f %6.2f %6.2f\n', cn{i}, g.Act, g.Cr, g.Cho, g.Lac);
end

%% ---- CONFIG (edit / MEASURE) -------------------------------------------
% doRelax = apply the relaxation correction below. FALSE = leave map as the raw
% LCModel water-scaled output. This phantom is UNDOPED aqueous (no Gd in recipe)
% -> long water relaxation, unmeasured metab relaxation. Keep FALSE until you
% MEASURE T1w/T2w + per-metabolite T1/T2, else you inject guessed bias.
doRelax = false;   % UNMEASURED T1/T2 -> guessed bias (~0.5-0.6x). Keep raw LCModel mM until measured.
if ~doRelax
    fprintf(['\napply_water_scaling_all: relaxation SKIPPED (doRelax=false). ' ...
             'map kept as LCModel water-scaled (raw mM).\n']);
    return
end

TE = ftSpec_smooth.te;      % ms
TR = ftSpec_smooth.tr;      % ms

% Undoped aqueous solution at 3T. Pure-water T1~3-4 s, T2~1-2 s; salts/O2 shorten
% T2. PLACEHOLDERS -> MEASURE (single-voxel IR for T1, multi-TE for T2).
T1w = 3000;  T2w = 1000;
% Metabolite [T1 T2] ms. UNMEASURED literature-ish guesses for aqueous solution.
relax = struct( ...
    'Act', [1600 1000], ...   % acetate CH3 singlet, long T2 in solution
    'Cr',  [1400  400], ...    % creatine (undoped, longer than in vivo)
    'Cho', [1300  400], ...    % choline
    'Lac', [1500  600]);       % lactate CH3
defaultT1T2 = [1400 500];
%% -----------------------------------------------------------------------

ATTH2O = (1 - exp(-TR/T1w)) * exp(-TE/T2w);

mets = fieldnames(map);
for i = 1:numel(mets)
    m = mets{i};
    if isfield(relax, m), t = relax.(m); else, t = defaultT1T2; end
    ATTMET  = (1 - exp(-TR/t(1))) * exp(-TE/t(2));
    factor  = ATTH2O / ATTMET;
    map.(m) = map.(m) * factor;
    fprintf('%-6s  T1=%4g T2=%4g  ATTMET=%.3f  factor=%.3f\n', m, t(1), t(2), ATTMET, factor);
end
fprintf('ATTH2O=%.3f.  map is now in mM (crlb/LW/SNR unchanged).\n', ATTH2O);
