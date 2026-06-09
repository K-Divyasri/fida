function MRSIStruct = setDimension(MRSIStruct, dimLabel, value)
%SETDIMENSION  Assign a new axis-index to one of the .dims fields.
%
%   MRSIStruct = setDimension(MRSIStruct, dimLabel, value)
%
%   • dimLabel  – string / char, case-insensitive (e.g. 'kx', 'kpetal')
%   • value     – positive integer (or 0 to clear)
%
%   All canonical dimension names are accepted; anything else throws
%   an error so typos are caught immediately.

    switch lower(dimLabel)

        % --- temporal / spectral ---------------------------------------
        case {'t','time'}
            MRSIStruct.dims.t  = value;
        case {'f','freq','frequency'}
            MRSIStruct.dims.f  = value;

        % --- spatial image axes ----------------------------------------
        case 'x'
            MRSIStruct.dims.x  = value;
        case 'y'
            MRSIStruct.dims.y  = value;
        case 'z'
            MRSIStruct.dims.z  = value;

        % --- k-space (Cartesian) ---------------------------------------
        case 'kx'
            MRSIStruct.dims.kx = value;
        case 'ky'
            MRSIStruct.dims.ky = value;
        case 'kz'
            MRSIStruct.dims.kz = value;

        % --- non-Cartesian extras --------------------------------------
        % case 'kpetal'
        %     MRSIStruct.dims.kpetal = value;
        case 'kshot'
            MRSIStruct.dims.kshot  = value;
        case 'no_loops'
            MRSIStruct.dims.no_loops = value;

        % --- acquisition loops -----------------------------------------
        case {'coils','coil'}
            MRSIStruct.dims.coils    = value;
        case {'averages','avg'}
            MRSIStruct.dims.averages = value;
        case 'timeinterleave'
            MRSIStruct.dims.timeinterleave = value;

        % --- spectroscopy-specific -------------------------------------
        case {'subspec','subspecs'}
            MRSIStruct.dims.subSpecs  = value;

        % --- catch-all -------------------------------------------------
        case 'extras'
            MRSIStruct.dims.extras = value;
        case 'kpts'
            MRSIStruct.dims.kpts = value;

        otherwise
            error('setDimension: unknown dimension label "%s".', dimLabel);
    end
end
