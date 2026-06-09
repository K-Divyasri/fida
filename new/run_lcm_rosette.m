if round(ftSpec_smooth.te)==9
    run_script_location = '~/Documents/Lubna/LCModel_run_scripts/3T_SpinEchoRosetteMRSI_svs_TE8p5.sh';
elseif round(ftSpec_smooth.te)==20
    run_script_location = '~/Documents/Lubna/LCModel_run_scripts/sneha3T_svs_TE20.sh';
end

lcmfoldername = ['FIDA_WRITELCM_',fileName];
    lcmfoldername = erase(lcmfoldername,'.dat');
    if ~exist(lcmfoldername,'dir')
        mkdir(lcmfoldername)
    end
    addpath(genpath(fullfile(pwd,'lcmfoldername')));
    cd(lcmfoldername)
    lcmoutput = erase(fileName,'.dat');
    lcmoutput = [lcmoutput '_out'];
    if ~exist(lcmoutput,'dir')
        mkdir(lcmoutput)
    end
    [loc_y,loc_x] = find(brain_area_raw==1);
    for k=1:size(loc_x)
        vox = op_CSItoMRS(ftSpec_smooth,loc_x(k),loc_y(k));
        vox_w = op_CSItoMRS(ftSpec_smooth_w,loc_x(k),loc_y(k));
        lcm_svs_file_name = [num2str(loc_x(k)) 'x' num2str(loc_y(k)) '_ftSpec_smooth_lcm'];
        lcm_svs_file_name = erase(lcm_svs_file_name,'.dat');
        lcm_svs_file_name_w = [lcm_svs_file_name '_w'];
        io_writelcm(vox,lcm_svs_file_name,vox.te);
        io_writelcm(vox_w,lcm_svs_file_name_w,vox_w.te);
        
        vox_lcm = [pwd '/' lcm_svs_file_name];
        vox_lcm_w = [pwd '/' lcm_svs_file_name_w];
        command = [run_script_location ' ' vox_lcm ' ' vox_lcm_w];
        system(command)
        %Plop it back into one folder and then you can make the crlb, met maps, and
        %lw and snr from .coord files!
        oldlocation = [lcm_svs_file_name '_out'];
        newlocation = [pwd '/' lcmoutput];
        variable_type = {'.table' '.coord' '.control' '.ps'};
        for i=1:size(variable_type,2)
            oldfile = [oldlocation '/' lcm_svs_file_name variable_type{i}];
            movefile(oldfile, newlocation)
        end
        rmdir(oldlocation)
    end

    final_location = [lcmfoldername '/' lcmoutput];

    cd ..
    [map,crlb,LW,SNR] = op_CSILCModelMaps(48,48,final_location,figure_folder_name=figurefoldername);
