function run_lcm_rosette_portable(fileName, ftSpec_smooth, ftSpec_smooth_w, mask, varargin)



    %---------------------------------------------------------------------%
    % 0. Parse name-value inputs
    %---------------------------------------------------------------------%
    cfg = parseInputs(varargin{:});

    %---------------------------------------------------------------------%
    % 1. Dims fix
    %---------------------------------------------------------------------%
    if ftSpec_smooth.dims.t == 0
        ftSpec_smooth.dims.t   = 1;  ftSpec_smooth.dims.f   = 0;
        ftSpec_smooth_w.dims.t = 1;  ftSpec_smooth_w.dims.f = 0;
    end

    %---------------------------------------------------------------------%
    % 2. Resolve working directory from fileName
    %---------------------------------------------------------------------%
    [workDir, fileBase, ~] = fileparts(fileName);
    if isempty(workDir), workDir = pwd; end
    cd(workDir);

    %---------------------------------------------------------------------%
    % 3. Resolve rawFolder — default to <workDir>/<fileBase>_LCModel_RAW
    %---------------------------------------------------------------------%
    if isempty(cfg.rawFolder)
        cfg.rawFolder = fullfile(workDir, [fileBase '_LCModel_RAW']);
        fprintf('rawFolder not specified. Defaulting to:\n        %s\n', cfg.rawFolder);
    end

    rawFolder = cfg.rawFolder;
    if ~exist(rawFolder, 'dir')
        fprintf('Creating RAW folder: %s\n', rawFolder);
        mkdir(rawFolder);
    end

    %---------------------------------------------------------------------%
    % 4. Generate bash runner script into rawFolder
    %---------------------------------------------------------------------%
    bashScript_path = fullfile(rawFolder, 'run_lcmodel_generated.sh');
    writeBashScript(bashScript_path, cfg);

    if cfg.useWSL
        bashScript_exec = win2wsl(bashScript_path);
        fixLineEndings_WSL(bashScript_exec);
    else
        bashScript_exec = bashScript_path;
        system(sprintf('chmod +x %s', shellQuote(bashScript_exec)));
    end

    fprintf('Generated bash script: %s\n', bashScript_path);

    %---------------------------------------------------------------------%
    % 5. Loop over masked voxels
    %---------------------------------------------------------------------%
    [loc_y, loc_x] = find(mask == 1);
    nvox = numel(loc_x);
    fprintf('Processing %d voxels.\n', nvox);

    for k = 1:nvox
        vx = loc_x(k);
        vy = loc_y(k);

        tag      = sprintf('%dx%d_ftSpec_smooth_lcm', vx, vy);
        rawfile1 = fullfile(rawFolder, tag);
        rawfile2 = [rawfile1 '_w'];

        fprintf('\n=== Voxel %d / %d  (%d,%d) ===\n', k, nvox, vx, vy);
        fprintf('  Writing RAW: %s\n', rawfile1);

        io_writelcm(op_CSItoMRS(ftSpec_smooth,   vx, vy), rawfile1, ftSpec_smooth.te);
        io_writelcm(op_CSItoMRS(ftSpec_smooth_w, vx, vy), rawfile2, ftSpec_smooth_w.te);

        %--- Build + run shell command ---%
        if cfg.useWSL
            rawfile1_exec = win2wsl(rawfile1);
            rawfile2_exec = win2wsl(rawfile2);
            inner = sprintf('%s %s %s %d %.8f', ...
                bashDQ(bashScript_exec), ...
                bashDQ(rawfile1_exec), ...
                bashDQ(rawfile2_exec), ...
                cfg.nunfil, cfg.deltat);
            cmd = sprintf('wsl bash -lc %s', bashSQ(inner));
        else
            cmd = sprintf('bash %s %s %s %d %.8f', ...
                shellQuote(bashScript_exec), ...
                shellQuote(rawfile1), ...
                shellQuote(rawfile2), ...
                cfg.nunfil, cfg.deltat);
        end

        fprintf('  CMD: %s\n', cmd);
        [status, out] = system(cmd);
        if status ~= 0
            warning('LCModel failed for voxel (%d,%d).\nOutput:\n%s', vx, vy, out);
        end

        %--- Move outputs into rawFolder ---%
        oldlocation = [rawfile1 '_out'];
        fileExt = {'.table', '.coord', '.control', '.ps', '.print', '.csv', '.log'};
        for e = 1:numel(fileExt)
            src = fullfile(oldlocation, [tag fileExt{e}]);
            if exist(src, 'file')
                movefile(src, fullfile(rawFolder, [tag fileExt{e}]));
            end
        end
        if exist(oldlocation, 'dir')
            rmdir(oldlocation, 's');
        end
    end

    fprintf('\nAll voxels processed. Outputs in: %s\n', rawFolder);
end


function cfg = parseInputs(varargin)

    % ---- Detect OS -------------------------------------------------------
    if ispc()
        os = 'Windows';
    elseif ismac()
        os = 'macOS';
    else
        os = 'Linux';
    end

    % ---- Defaults --------------------------------------------------------
    cfg.rawFolder  = '';
    cfg.basisFile  = '';
    cfg.lcmodelBin = '~/.lcmodel/bin/lcmodel';
    cfg.ownerStr   = '';
    cfg.hzpppm     = 123.25;
    cfg.deltat     = 0.00063;
    cfg.nunfil     = 576;
    cfg.echot      = 15;
    cfg.ppmst      = 4.25;
    cfg.ppmend     = 0.2;
    cfg.wconc      = 55556;   % water concentration for water scaling (mM); pure water
    cfg.useWSL     = ispc();

    % ---- Parse name-value pairs ------------------------------------------
    validKeys = fieldnames(cfg);
    i = 1;
    while i <= numel(varargin)
        key = varargin{i};
        if ~ischar(key) && ~isstring(key)
            error('Expected a parameter name (string) at position %d, got %s.', ...
                  i + 4, class(key));
        end
        key = char(key);
        if ~any(strcmpi(key, validKeys))
            error('Unknown parameter ''%s''.\nValid parameters: %s', ...
                  key, strjoin(validKeys, ', '));
        end
        if i + 1 > numel(varargin)
            error('Parameter ''%s'' has no value.', key);
        end
        idx = find(strcmpi(key, validKeys), 1);
        cfg.(validKeys{idx}) = varargin{i+1};
        i = i + 2;
    end


    if isempty(cfg.basisFile)
        error(['basisFile is required. Pass it as a name-value pair, e.g.:\n' ...
               '  run_lcm_rosette_portable(..., ''basisFile'', ''/path/to/file.basis'')']);
    end
    if isempty(cfg.ownerStr)
        error(['ownerStr is required. Pass it as a name-value pair, e.g.:\n' ...
               '  run_lcm_rosette_portable(..., ''ownerStr'', ''Your Name, Institution'')']);
    end

    cfg.licenseKey = '210387309';

    if cfg.useWSL
        cfg.basisFile_exec  = win2wsl(cfg.basisFile);
        cfg.lcmodelBin_exec = cfg.lcmodelBin;   % always a WSL/POSIX path
    else
        if looksLikeWindowsPath(cfg.basisFile)
            warning(['basisFile looks like a Windows path (%s) ' ...
                     'but OS is %s. Did you mean to set useWSL=true?'], cfg.basisFile, os);
        end
        cfg.basisFile_exec  = cfg.basisFile;
        cfg.lcmodelBin_exec = cfg.lcmodelBin;
    end

    fprintf('OS: %s  |  WSL mode: %d\n', os,                  cfg.useWSL);
    fprintf('basisFile  → %s\n',           cfg.basisFile_exec);
    fprintf('lcmodelBin → %s\n',           cfg.lcmodelBin_exec);
    fprintf('ownerStr   → %s\n',           cfg.ownerStr);
end


%=========================================================================%
%  writeBashScript  
%=========================================================================%
function writeBashScript(outPath, cfg)

    fid = fopen(outPath, 'w');
    if fid == -1
        error('Cannot write bash script to: %s', outPath);
    end

    wl = @(s) fprintf(fid, '%s\n', s);

    wl('#!/usr/bin/env bash');
    wl('set -e');
    wl('if ( set -o 2>/dev/null | grep -q pipefail ); then set -o pipefail; fi');
    wl('');
    wl('# Auto-generated by run_lcm_rosette_portable.m — do not edit');
    wl('# Args: $1=WS stem  $2=H2O stem  $3=NUNFIL(opt)  $4=DELTAT(opt)');
    wl('');

    wl(sprintf('LCMODEL_BIN="%s"',  cfg.lcmodelBin_exec));
    wl(sprintf('BASIS="%s"',        cfg.basisFile_exec));
    wl(sprintf('LICENSE_KEY="%s"',  cfg.licenseKey));
    wl(sprintf('OWNER="%s"',        strrep(cfg.ownerStr, '"', '\"')));
    wl(sprintf('HZPPPM="%g"',       cfg.hzpppm));
    wl(sprintf('ECHOT="%g"',        cfg.echot));
    wl(sprintf('PPMST="%g"',        cfg.ppmst));
    wl(sprintf('PPMEND="%g"',       cfg.ppmend));
    wl(sprintf('WCONC="%g"',        cfg.wconc));
    wl('');

    wl('WS_STEM="${1:?need WS stem}"');
    wl('W_STEM="${2:?need water stem}"');
    wl(sprintf('NUNFIL="${3:-%d}"',   cfg.nunfil));
    wl(sprintf('DELTAT="${4:-%.8f}"', cfg.deltat));
    wl('');

    wl('datadir="$(dirname "$WS_STEM")"');
    wl('series="$(basename "$WS_STEM")"');
    wl('outdir="${datadir}/${series}_out"');
    wl('mkdir -p "$outdir"');
    wl('');

    wl('pick_path() {');
    wl('  local stem="$1"');
    wl('  if   [[ -f "${stem}"     ]]; then printf "%s\n" "${stem}"');
    wl('  elif [[ -f "${stem}.RAW" ]]; then printf "%s\n" "${stem}.RAW"');
    wl('  elif [[ -f "${stem}.raw" ]]; then printf "%s\n" "${stem}.raw"');
    wl('  else return 1; fi');
    wl('}');
    wl('WS_PATH="$(pick_path "$WS_STEM")" || { echo "ERROR: no RAW for ${WS_STEM}"; exit 1; }');
    wl('if W_PATH="$(pick_path "$W_STEM")"; then ECC=1; else ECC=0; W_PATH=""; fi');
    wl('');

    wl('controlfile="${outdir}/${series}.control"');
    wl('logfile="${outdir}/${series}.log"');
    wl('');

    wl('{');
    wl('  echo ''$LCMODL''');
    wl('  echo "OWNER=''${OWNER}''"');
    wl('  echo "KEY=${LICENSE_KEY}"');
    wl('  echo "TITLE=''onepulse ${series}''"');
    wl('  echo "FILRAW=''${WS_PATH}''"');
    wl('  echo "FILPS=''${outdir}/${series}.ps''"');
    wl('  echo "FILTAB=''${outdir}/${series}.table''"');
    wl('  echo "FILCOO=''${outdir}/${series}.coord''"');
    wl('  echo "FILCSV=''${outdir}/${series}.csv''"');
    wl('  echo "FILBAS=''${BASIS}''"');
    wl('  echo "FILPRI=''${outdir}/${series}.print''"');
    wl('  echo "HZPPPM=${HZPPPM}"');
    wl('  echo "DELTAT=${DELTAT}"');
    wl('  echo "NUNFIL=${NUNFIL}"');
    wl('  echo "ECHOT=${ECHOT}"');
    wl('  echo "LPRINT=6"');
    wl('  echo "LTABLE=7"');
    wl('  echo "LCOORD=9"');
    wl('  echo "LCSV=11"');
    wl('  echo "DOREFS(1)=T"');
    wl('  echo "LPS=8"');
    wl('  echo "PPMST=${PPMST}"');
    wl('  echo "PPMEND=${PPMEND}"');
    wl('  echo "RFWHM=0.15"');
    wl('  echo "WDLINE(1)=0.025"');
    wl('  echo "DKNTMN=0.3"');
    wl('  echo "NEACH=999"');
    wl('  echo "NOMIT=14"');
    wl('  echo "CHOMIT(1)=''bHB''"');
    wl('  echo "CHOMIT(2)=''bHG''"');
    wl('  echo "CHOMIT(3)=''Lip13a''"');
    wl('  echo "CHOMIT(4)=''Lip13b''"');
    wl('  echo "CHOMIT(5)=''Lip13c''"');
    wl('  echo "CHOMIT(6)=''Lip13d''"');
    wl('  echo "CHOMIT(7)=''Lip09''"');
    wl('  echo "CHOMIT(8)=''MM09''"');
    wl('  echo "CHOMIT(9)=''Lip20''"');
    wl('  echo "CHOMIT(10)=''MM20''"');
    wl('  echo "CHOMIT(11)=''MM12''"');
    wl('  echo "CHOMIT(12)=''MM14''"');
    wl('  echo "CHOMIT(13)=''MM17''"');
    wl('  echo "CHOMIT(14)=''-CrCH2''"');
    wl('  echo "NNORAT=7"');
    wl('  echo "NORATO(1)=''Lip09/Lip13*''"');
    wl('  echo "NORATO(2)=''Lip20/Lip13*''"');
    wl('  echo "NORATO(3)=''MM20/MM09*''"');
    wl('  echo "NORATO(4)=''MM12/MM09*''"');
    wl('  echo "NORATO(5)=''MM14/MM09*''"');
    wl('  echo "NORATO(6)=''MM17/MM09*''"');
    wl('  echo "NORATO(7)=''-CrCH2/totCr''"');
    wl('  echo "NRATIO=8"');
    wl('  echo "CHRATO(1)=''MM2/MM1 = 0.51 +- 0.17''"');
    wl('  echo "CHRATO(2)=''MM3/MM1 = 1.05 +- 0.63''"');
    wl('  echo "CHRATO(3)=''MM4/MM1 = 1.29 +- 1.03''"');
    wl('  echo "CHRATO(4)=''MM5/MM1 = 3.16 +- 0.79''"');
    wl('  echo "CHRATO(5)=''MM6/MM1 = 0.63 +- 0.16''"');
    wl('  echo "CHRATO(6)=''MM7/MM1 = 0.54 +- 0.27''"');
    wl('  echo "CHRATO(7)=''MM8/MM1 = 0.33 +- 0.17''"');
    wl('  echo "CHRATO(8)=''MM9/MM1 = 1.40 +- 0.7''"');
    wl('  if [[ "$ECC" -eq 1 ]]; then');
    wl('    echo "DOECC=T"');
    wl('    echo "DOWS=T"');
    wl('    echo "FILH2O=''${W_PATH}''"');
    wl('    echo "ATTH2O=1.0"');
    wl('    echo "WCONC=${WCONC}"');
    wl('  fi');
    wl('  echo ''$END''');
    wl('} > "$controlfile"');
    wl('');
    wl('echo "Control file: $controlfile"');
    wl('"${LCMODEL_BIN}" < "$controlfile" > "$logfile" 2>&1 || true');
    wl('echo "LCModel done"');
    wl('');
    wl('missing=0');
    wl('for f in "${outdir}/${series}.table" "${outdir}/${series}.coord" "${outdir}/${series}.csv"; do');
    wl('  if [[ ! -f "$f" ]]; then echo "WARN: missing $(basename "$f")"; missing=1; fi');
    wl('done');
    wl('if [[ "$missing" -eq 1 ]]; then');
    wl('  echo "----- LCModel log tail -----"');
    wl('  tail -n 120 "$logfile" || true');
    wl('  echo "----------------------------"');
    wl('fi');

    fclose(fid);
    fprintf('[INFO] Bash script written: %s\n', outPath);
end


%=========================================================================%
%  Path / quoting helpers
%=========================================================================%
function linuxPath = win2wsl(winPath)
    if isstring(winPath), winPath = char(winPath); end
    if strncmp(winPath, '/mnt/', 5)
        linuxPath = winPath; return;
    end
    wp = strrep(winPath, '\', '/');
    if numel(wp) >= 2 && wp(2) == ':'
        linuxPath = ['/mnt/' lower(wp(1)) wp(3:end)];
        return;
    end
    linuxPath = wp;
end

function tf = looksLikeWindowsPath(p)
    if isstring(p), p = char(p); end
    tf = numel(p) >= 2 && isletter(p(1)) && p(2) == ':';
end

function q = bashDQ(s)
    if isstring(s), s = char(s); end
    q = ['"' strrep(s, '"', '\"') '"'];
end

function q = bashSQ(s)
    if isstring(s), s = char(s); end
    q = ['''' strrep(s, '''', '''"''"''') ''''];
end

function q = shellQuote(s)
    if isstring(s), s = char(s); end
    q = ['''' strrep(s, '''', '''\''''''') ''''];
end

function fixLineEndings_WSL(scriptWSL)
    tryCmd = sprintf('wsl bash -lc %s', bashSQ( ...
        sprintf('command -v dos2unix >/dev/null 2>&1 && dos2unix %s >/dev/null 2>&1; chmod +x %s >/dev/null 2>&1', ...
        bashDQ(scriptWSL), bashDQ(scriptWSL))));
    system(tryCmd);
end