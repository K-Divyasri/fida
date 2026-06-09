%Make metabolite maps again using svs method of lcmodel processing

function [map,crlb,LW,SNR] = op_CSILCModelMaps(size_x,size_y,lcmodel_output_location,options)
    arguments
        size_x (1,:) double
        size_y (1,:) double
        lcmodel_output_location (1,:) char
        options.figure_folder_name (1,:) char = "false"
    end

   %Create the structure of the maps we want
    cd(lcmodel_output_location)
    %map_list = {'Cr+PCr' 'GPC+PCh' 'NAA+NAAG' 'Glu+Gln' 'Ins'};
    map_list = {'Lac' 'Ace' 'Cr' 'Cho'};
    map_list_readable = erase(map_list,'+');

    % Initialise one map + crlb field per metabolite in map_list, so this
    % works for ANY basis (brain set, phantom Lac/Ace/Cr/Cho, etc.).
    % Pre-sizing to [size_x size_y] is essential: otherwise the first
    % index-assignment below creates an undersized field and the later
    % reshape() in the figure section fails ("Number of elements must
    % not change").
    map  = struct();
    crlb = struct();
    for m = 1:numel(map_list_readable)
        map.(map_list_readable{m})  = zeros(size_x, size_y);
        crlb.(map_list_readable{m}) = zeros(size_x, size_y);
    end

    LW = zeros(size_x,size_y);
    SNR = zeros(size_x,size_y);

    %Import the lcmodel output .table files
    fileList = dir('*.table');
    for s=1:size(fileList,1)%For each voxel inside of the brain
        table_name = fileList(s).name;

        %Clarify which voxel we're at
        location = fileList(s).name(1:5);
        location = erase(location,'_');%Get rid of variables
        location = erase(location,'f');%Get rid of variables
        location = split(location,'x');

        location_x = str2num(location{1});
        location_y = str2num(location{2});
        %% For SNR and LW

        raw_text = fileread(table_name);

        % Skip voxels where LCModel failed: a FATAL run writes a .table with
        % only $$DIAG/$$INPU sections — no concentration ($$CONC) or misc
        % ($$MISC) tables. There is nothing to read, so leave the voxel at 0
        % instead of crashing on an empty parse result.
        if ~contains(raw_text,'$$CONC')
            warning('op_CSILCModelMaps:failedVoxel', ...
                'No concentration table in %s (LCModel failed) — voxel %dx%d left as 0.', ...
                table_name, location_x, location_y);
            continue
        end

        start_integer = strfind(raw_text,'FWHM =');   % $$MISC line, not 'RFWHM='
        end_integer   = strfind(raw_text,'shift');
        if ~isempty(start_integer) && ~isempty(end_integer) ...
                && start_integer(1) < end_integer(1)
            misc_text = raw_text(start_integer(1):end_integer(1));
            misc_text = split(misc_text,["=","ppm","Data"]);
            if numel(misc_text) >= 4
                LW(location_y,location_x)  = str2num(misc_text{2});
                SNR(location_y,location_x) = str2num(misc_text{4});
            end
        end
        %% For metabolite maps and CRLB
        try
            T_vox = importdata(table_name);
        catch ME
            warning('op_CSILCModelMaps:importdataFailed', ...
                'importdata failed on %s (%s) — voxel %dx%d left as 0.', ...
                table_name, ME.message, location_x, location_y);
            continue
        end
        if ~isstruct(T_vox) || ~isfield(T_vox,'textdata')
            continue   % no parseable text block — nothing more to extract
        end
        for l=1:size(T_vox.textdata,1)
            lines_in_table = split(T_vox.textdata{l},["  "," "]);
            lines_in_table(ismember(lines_in_table,'')) = [];
            for tf=1:size(lines_in_table,1)
                for map_item=1:size(map_list,2)
                    yn = strcmp(lines_in_table{tf},map_list{map_item});
                    if yn==1 && numel(lines_in_table) >= 2
                        map.(map_list_readable{map_item})(location_y,location_x)=str2num(lines_in_table{1});
                        crlb.(map_list_readable{map_item})(location_y,location_x)=str2num(erase(lines_in_table{2},'%'));
                    end
                end
            end
        end
    end



    %Save these figures
    if(options.figure_folder_name)
        cd ../..
        for f=1:size(map_list_readable,2)
            f1=figure('visible','off');
            imagesc(map.(map_list_readable{f}))
            colorbar('Color','w','FontSize',20)
            colormap hot
            title( ["Metabolite Map for" map_list{f}],'Color','w')
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf,'Color',[0 0 0]);
            axis = gca;
            axis.YColor = 'w';
            axis.XColor = 'w';
            %saveas(f1,fullfile(pwd,options.figure_folder_name,map_list_readable{f}),'jpg');
            f2=figure('visible','off');
            crlb_rs = reshape(crlb.(map_list_readable{f}),[1,size_x*size_y]);
            crlb_pass = ((size(crlb_rs(crlb_rs>0 & crlb_rs < 20))) / size(crlb_rs(crlb_rs>0)) )*100;
            crlb_pass = ['Percent CRLB that passes is ' num2str(crlb_pass)];
            imagesc(crlb.(map_list_readable{f}))
            colorbar('Color','w','FontSize',20)
            colormap hot
            title( ["CRLB Map for" map_list{f}],'Color','w')
            subtitle(crlb_pass,'Color','w')
            clim([0 50]);
            crlb_name = strcat('crlb_', map_list_readable{f});
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf,'Color',[0 0 0]);
            axis = gca;
            axis.YColor = 'w';
            axis.XColor = 'w';
            %saveas(f2,fullfile(pwd,options.figure_folder_name,crlb_name),'jpg');
        end
        f3=figure('visible','off');
        LW_rs = reshape(LW,[1,size_x*size_y]);
        LW_mean=((size(LW_rs(LW_rs < 0.1 & LW_rs>0))) / size(LW_rs(LW_rs>0)) )*100;
        LW_mean = ['Percent LW that passes is ' num2str(LW_mean)];
        imagesc(LW)
        colorbar('Color','w','FontSize',20)
        colormap hot
        clim([0 0.5]);
        title("LW Map in ppm, from LCModel",'Color','w')
        subtitle(LW_mean,'Color','w')
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf,'Color',[0 0 0]);
        axis = gca;
        axis.YColor = 'w';
        axis.XColor = 'w';
        %saveas(f3,fullfile(pwd,options.figure_folder_name,'LW_map_lcm'),'jpg');

        f4=figure('visible','off');
        SNR_rs = reshape(SNR,[1,size_x*size_y]);
        SNR_pass=((size(SNR_rs(SNR_rs > 3))) / size(SNR_rs(SNR_rs>0)) )*100;
        SNR_pass = ['Percent SNR that passes is ' num2str(SNR_pass)];
        imagesc(SNR)
        colorbar('Color','w','FontSize',20)
        colormap hot
        clim([0 10]);
        title("SNR Map, from LCModel",'Color','w')
        subtitle(SNR_pass,'Color','w')
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf,'Color',[0 0 0]);
        axis = gca;
        axis.YColor = 'w';
        axis.XColor = 'w';
        %saveas(f4,fullfile(pwd,options.figure_folder_name,'SNR_map_lcm'),'jpg');
        cd(lcmodel_output_location);
    end
end
