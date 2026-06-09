function MRSIStruct = op_CSIAverage(MRSIStruct)
    flag = checkArgumentsEnhanced(MRSIStruct);
    %if check arguments returns true, abort
    if(flag); return; end

    fprintf('Starting averaging operation...\n');
    fprintf('Input data size: %s\n', mat2str(MRSIStruct.sz));
    fprintf('Input dims:\n');
    disp(MRSIStruct.dims);

    data = getData(MRSIStruct);

    averageDimension = getDimension(MRSIStruct, 'averages');
    fprintf('Averaging along dimension %d\n', averageDimension);
    
    data = squeeze(mean(data, averageDimension));
    fprintf('Data size after averaging: %s\n', mat2str(size(data)));

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = setFlags(MRSIStruct, 'averaged',  1);
    MRSIStruct = removeDimension(MRSIStruct, 'averages');
    
    fprintf('Final data size: %s\n', mat2str(MRSIStruct.sz));
    fprintf('Final dims:\n');
    disp(MRSIStruct.dims);
    fprintf('Averaging completed.\n\n');
end

% Enhanced check function with better error reporting
function flag = checkArgumentsEnhanced(in)
    flag = false;
    
    % Check if averages dimension exists
    if ~isfield(in.dims, 'averages') || in.dims.averages == 0
        fprintf('No averages dimension found!\n');
        fprintf('Current dims structure:\n');
        disp(in.dims);
        fprintf('Data size: %s\n', mat2str(in.sz));
        disp('Aborting averaging operation');
        flag = true;
        return;
    end
    
    % Check if already averaged
    if isfield(in.flags, 'averaged') && in.flags.averaged
        disp('Already averaged! Aborting')
        flag = true;
        return;
    end
    
    % Check if averages dimension has valid size
    avgDim = in.dims.averages;
    if avgDim > length(in.sz) || in.sz(avgDim) <= 1
        fprintf('Invalid averages dimension: dim=%d, size=%d\n', avgDim, in.sz(avgDim));
        fprintf('Cannot average dimension with size <= 1\n');
        flag = true;
        return;
    end
    
    fprintf('Averages dimension check passed: dim=%d, size=%d\n', avgDim, in.sz(avgDim));
end