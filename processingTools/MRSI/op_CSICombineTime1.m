function MRSIStruct = op_CSICombineTime1(MRSIStruct, dimensionLabel)
    % If Cartesian, no need to combine dimensions — return unchanged
    if isfield(MRSIStruct, 'flags') && isfield(MRSIStruct.flags, 'isCartesian') && MRSIStruct.flags.isCartesian
        return;
    end

    % Compute new time dimension
    newTimeLength = prod(getSizeFromDimensions(MRSIStruct, {'t', dimensionLabel}));

    % Reshape dimensions
    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t', dimensionLabel});
    data = getData(MRSIStruct);

    % Collapse along time and specified dimension
    data = reshape(data, newTimeLength, []);
    MRSIStruct = setData(MRSIStruct, data);

    % Adjust sizes and permute back
    prevSize(2) = [];
    prevSize(1) = newTimeLength;
    prevPermute = removeDimPrevPermute(prevPermute, 2);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end
