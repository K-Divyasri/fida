function dComp = op_CSIPSFCorrection(MRSIStruct, kFile, method, varargin)
%OP_CSIPSFCORRECTION  Unified density-compensation dispatcher for CSI.
%
%   dComp = op_CSIPSFCorrection(MRSIStruct, kFile, method)
%   dComp = op_CSIPSFCorrection(MRSIStruct, kFile, method, Name, Value, ...)
%
%   method:
%     'nn'         -> op_CSIPSFCorrection_nn  (k-nearest-neighbour)
%     'voronoi'    -> op_CSIPSFCorrection_v   (Voronoi cell area)
%     'pipe_menon' -> op_CSIPSFCorrection_pm  (Pipe-Menon iterative)
%
%   Any extra Name/Value arguments are forwarded to the chosen backend
%   (e.g. 'numNeighbors', 'numIterations', 'modelType', 'isPlotWeights').

    if nargin < 3 || isempty(method)
        method = 'nn';
    end

    methodLower = lower(string(method));
    switch methodLower
        case "nn"
            dComp = op_CSIPSFCorrection_nn(MRSIStruct, kFile, varargin{:});
        case {"voronoi", "v"}
            dComp = op_CSIPSFCorrection_v(MRSIStruct, kFile, varargin{:});
        case {"pipe_menon", "pipemenon", "pm"}
            dComp = op_CSIPSFCorrection_pm(MRSIStruct, kFile, varargin{:});
        otherwise
            error('op_CSIPSFCorrection:badMethod', ...
                'Unknown DCF method "%s". Use ''nn'', ''voronoi'', or ''pipe_menon''.', method);
    end
end
