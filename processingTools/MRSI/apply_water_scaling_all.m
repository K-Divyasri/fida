%APPLY_WATER_SCALING_ALL  Rescale the workspace `map` (LCModel raw) to mM.
%
%   Assumes map, crlb, LW, SNR are already in the workspace. Rescales `map` to
%   mM in place; crlb/LW/SNR are untouched. Run ONCE (re-running double-scales).
%
%   WHAT THIS DOES
%   LCModel already turned each metabolite into a concentration by dividing its
%   signal by the water-reference signal and multiplying by WCONC (=55556, the
%   right value for a pure-water phantom). But it assumed both signals were
%   FULLY RELAXED. They are not: at finite TE and TR the signal is shrunk by
%   T1 (recovery) and T2 (decay). So we divide that shrinkage back out.
%
%       ATTH2O = (1 - exp(-TR/T1w)) * exp(-TE/T2w)     % water shrinkage
%       ATTMET = (1 - exp(-TR/T1m)) * exp(-TE/T2m)     % metabolite shrinkage
%       true mM = raw * ATTH2O / ATTMET
%
%   Works for any basis: it loops over whatever metabolites are in `map` and
%   looks up each one's relaxation in the table below. Unlisted metabolites get
%   the default and print a warning.

%% ---- CONFIG (edit / MEASURE) -------------------------------------------
% doRelax = apply the T1/T2 relaxation correction below.
%   false (default): leave map as LCModel's water-scaled output (raw). Use this
%     until T1w/T2w and per-metabolite T1/T2 are MEASURED for THIS phantom.
%     With unmeasured placeholders the correction multiplies every metabolite by
%     ATTH2O/ATTMET (~0.66 for the pure-water guesses -> values ~halved), even
%     though the raw values already match ground truth.
%   true: apply the correction (only trustworthy with measured relaxation).
doRelax = false;   % raw LCModel water-scaled output was closest to GT; guessed
                   % relaxation over-shrinks (~0.83x). Set true ONLY with
                   % MEASURED phantom T1w/T2w + per-metabolite T1/T2.
if ~doRelax
    fprintf(['apply_water_scaling_all: relaxation correction SKIPPED ' ...
             '(doRelax=false). map kept as LCModel water-scaled (raw).\n' ...
             'Set doRelax=true only after measuring T1w/T2w + per-metabolite ' ...
             'T1/T2 for this phantom.\n']);
    return
end

TE = ftSpec_smooth.te;      % ms  (or type a number, e.g. 35)
TR = ftSpec_smooth.tr;      % ms  (or type a number, e.g. 2000)

% Water reference relaxation (ms). GE Braino is GADOLINIUM-doped so its
% relaxation is set NEAR BRAIN TISSUE (per GE), NOT pure water. Gd content
% varies per batch and DEGRADES with age -> MEASURE for publication. These
% brain-3T-like values give factor ~1 (matches the raw==GT observation).
T1w = 1200;  T2w = 250;

% Per-metabolite [T1 T2] in ms. GE Braino is Gd-doped to ~in-vivo brain 3T, so
% these are literature 3T brain values (Mlynarik/Wyss). MEASURE/cite for
% publication; batch + age shift them.
relax = struct( ...
    'NAA', [1350 295], ...
    'NAAG',[1350 295], ...
    'Cr',  [1240 152], ...
    'PCr', [1240 152], ...
    'Cho', [1080 218], ...
    'Ins', [1010 197], ...
    'Glu', [1270 180], ...
    'Gln', [1270 180], ...
    'Lac', [1500 250], ...
    'Act', [1500 1000], ...
    'Tau', [1200 120], ...
    'GABA',[1310  88], ...
    'Scyllo',[1010 197]);
defaultT1T2 = [1250 200];   % used for any metabolite not in `relax`
%% -----------------------------------------------------------------------

ATTH2O = (1 - exp(-TR/T1w)) * exp(-TE/T2w);   % one water number

mets = fieldnames(map);                        % scale whatever is in the map
for i = 1:numel(mets)
    m = mets{i};
    if isfield(relax, m)
        t = relax.(m);
    else
        t = defaultT1T2;
        fprintf('  (!) %s not in relax table -- using default [%g %g]\n', m, t(1), t(2));
    end
    ATTMET  = (1 - exp(-TR/t(1))) * exp(-TE/t(2));
    factor  = ATTH2O / ATTMET;
    map.(m) = map.(m) * factor;                % raw -> mM
    fprintf('%-6s  T1=%4g T2=%4g  ATTMET=%.3f  factor=%.3f\n', m, t(1), t(2), ATTMET, factor);
end

fprintf('ATTH2O=%.3f.  map is now in mM (crlb/LW/SNR unchanged).\n', ATTH2O);
