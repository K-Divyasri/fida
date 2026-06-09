    function MRSIStruct = reshape_twix_data(MRSIStruct,kFile)
        data = MRSIStruct.data;
        sz = MRSIStruct.sz;
        [kTable, ~] = readKFile(kFile);
        % Spectral updates
        kPtsPerCycle = getKPtsPerCycle(kTable);
        if numel(sz) == 4
            dim1 = kPtsPerCycle;
            dim2 = sz(1) / kPtsPerCycle;
            dim3 = sz(2);
            dim4 = sz(3);
            dim5 = sz(4);
    
            data_resh = reshape(data, [dim1, dim2, dim3, dim4, dim5]);
            data_resh = permute(data_resh, [2, 3, 4, 1, 5]);
            data_resh = reshape(data_resh, [dim2, dim3, dim4, dim1, dim5]);
            MRSIStruct.dims.kpts =4;
            MRSIStruct.dims.kshot =5;
        elseif numel(sz) == 3
            dim1 = kPtsPerCycle;
            dim2 = sz(1) / kPtsPerCycle;
            dim3 = sz(2);
            dim4 = sz(3);
    
            data_resh = reshape(data, [dim1, dim2, dim3, dim4]);
            data_resh = permute(data_resh, [2, 3, 1, 4]);
            data_resh = reshape(data_resh, [dim2, dim3, dim1, dim4]);
            MRSIStruct.dims.kpts =3;
            MRSIStruct.dims.kshot =4;
        else
            data_resh = data;
        end
    
        MRSIStruct.data = data_resh;
        MRSIStruct.sz = size(data_resh);
    end