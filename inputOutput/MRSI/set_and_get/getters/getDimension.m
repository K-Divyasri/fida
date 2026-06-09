function dimNumber = getDimension(MRSIStruct, dimLabel)
%GETDIMENSION   Return the axis index for a given dimension label.
%
%   dim = getDimension(MRSIStruct, dimLabel)
%
%   • dimLabel is case-insensitive; common aliases are accepted.
%   • If the field was never assigned (==0) you’ll get 0 back.
%   • An unrecognised label triggers an error, so typos don’t go unnoticed.

    arguments
        MRSIStruct (1,1) struct
        dimLabel           char
    end

    switch lower(dimLabel)

        % --- temporal / spectral ---------------------------------------
        case {'t','time'}
            dimNumber = MRSIStruct.dims.t;
        case {'f','freq','frequency'}
            dimNumber = MRSIStruct.dims.f;

        % --- image-space axes ------------------------------------------
        case 'x'
            dimNumber = MRSIStruct.dims.x;
        case 'y'
            dimNumber = MRSIStruct.dims.y;
        case 'z'
            dimNumber = MRSIStruct.dims.z;

        % --- k-space (Cartesian) ---------------------------------------
        case 'kx'
            dimNumber = MRSIStruct.dims.kx;
        case 'ky'
            dimNumber = MRSIStruct.dims.ky;
        case 'kz'
            dimNumber = MRSIStruct.dims.kz;

        % --- non-Cartesian extras --------------------------------------
        % case 'kpetal'
        %     dimNumber = MRSIStruct.dims.kpetal;
        case 'kshot'
            dimNumber = MRSIStruct.dims.kshot;
        case 'no_loops'
            dimNumber = MRSIStruct.dims.no_loops;

        % --- acquisition loops -----------------------------------------
        case {'coils','coil'}
            dimNumber = MRSIStruct.dims.coils;
        case {'averages','average','avg'}
            dimNumber = MRSIStruct.dims.averages;
        case 'timeinterleave'
            dimNumber = MRSIStruct.dims.timeinterleave;

        % --- spectroscopy-specific -------------------------------------
        case {'subspec','subspecs'}
            dimNumber = MRSIStruct.dims.subSpecs;

        % --- catch-all -------------------------------------------------
        case {'extras','extra'}
            dimNumber = MRSIStruct.dims.extras;
        case 'kpts'
            dimNumber = MRSIStruct.dims.kpts;

        otherwise
            error('getDimension: unknown dimension label "%s".', dimLabel);
    end
end
