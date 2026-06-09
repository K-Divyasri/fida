function io_CSIwritelcmcontrolfile_separate(refcontrolfile,LCDir,lcmodel_name,MRSIStruct,basisfile,LCM_indices,cfile_nb)
    fid = fopen(refcontrolfile) ; 
    S = textscan(fid,'%s','delimiter','\n') ; 
    fclose(fid) ; 
    S = S{1} ; 

    % Find the string you want to change 
    idx = contains(S,'TITLE') ; 
    % Replace 
    S(idx) =  {strcat('TITLE= ',char(39),'TE/TR/NS=',num2str(MRSIStruct.te),'/',num2str(MRSIStruct.tr),'/',num2str(MRSIStruct.sz(1)),' Scan=',MRSIStruct.seq,' Scan Date=',char(MRSIStruct.scanDate),char(39))} ; 
    % Find the string you want to change 
    idx = contains(S,'NUNFIL') ; 
    % Replace 
    S(idx) =  {strcat('NUNFIL= ',num2str(MRSIStruct.sz(1)))} ; 
    % Find the string you want to change 
    idx = contains(S,'DELTAT') ; 
    % Replace 
    S(idx) =  {strcat('DELTAT= ',num2str(MRSIStruct.spectralDwellTime))} ; 
    % Find the string you want to change 
    idx = contains(S,'ECHOT') ; 
    % Replace 
    S(idx) =  {strcat('ECHOT= ',num2str(MRSIStruct.te))} ; 
    % Find the string you want to change 
    idx = contains(S,'NDCOLS') ; 
    % Replace 
    S(idx) =  {strcat('NDCOLS= ',num2str(MRSIStruct.sz(2)))} ; 
    % Find the string you want to change 
    idx = contains(S,'NDROWS') ; 
    % Replace 
    S(idx) =  {strcat('NDROWS= ',num2str(MRSIStruct.sz(3)))} ;     
    % Find the string you want to change 
    idx = contains(S,'ICOLST') ; 
    % Replace 
    %S(idx) =  {strcat('ICOLST= ',num2str(1))} ; 
    S(idx) =  {strcat('ICOLST= ',num2str(LCM_indices(1,1)))} ; 
    % Find the string you want to change 
    idx = contains(S,'IROWST') ; 
    % Replace 
    %S(idx) =  {strcat('IROWST= ',num2str(1))} ; 
    S(idx) =  {strcat('IROWST= ',num2str(LCM_indices(2,1)))} ; 
    % Find the string you want to change 
    idx = contains(S,'ICOLEN') ; 
    % Replace 
    %S(idx) =  {strcat('ICOLEN= ',num2str(MRSIStruct.sz(2)))} ; 
    S(idx) =  {strcat('ICOLEN= ',num2str(LCM_indices(1,2)))} ; 
    % Find the string you want to change 
    idx = contains(S,'IROWEN') ; 
    % Replace 
    %S(idx) =  {strcat('IROWEN= ',num2str(MRSIStruct.sz(3)))} ; 
    S(idx) =  {strcat('IROWEN= ',num2str(LCM_indices(2,2)))} ; 
    % Find the string you want to change 
    idx = contains(S,'FILBAS') ; 
    % Replace 
    S(idx) =  {strcat('FILBAS= ',char(39),basisfile,char(39))} ;     
    % Find the string you want to change 
    idx = contains(S,'FILRAW') ; 
    % Replace 
    S(idx) =  {strcat('FILRAW= ',char(39),strcat(LCDir,'/',lcmodel_name),char(39))} ; 
    % Find the string you want to change 
    idx = contains(S,'FILH2O') ; 
    % Replace 
    S(idx) =  {strcat('FILH2O= ',char(39),strcat(LCDir,'/',lcmodel_name,'_w'),char(39))} ; 
    % Find the string you want to change 
    idx = contains(S,'FILCSV') ; 
    % Replace 
    S(idx) =  {strcat('FILCSV= ',char(39),strcat(LCDir,'/',lcmodel_name,'_out','/',lcmodel_name,num2str(cfile_nb),'.csv'),char(39))} ; 
    % Find the string you want to change 
    idx = contains(S,'FILTAB') ; 
    % Replace 
    S(idx) =  {strcat('FILTAB= ',char(39),strcat(LCDir,'/',lcmodel_name,'_out','/',lcmodel_name,'.table'),char(39))} ; 
    % Find the string you want to change 
    idx = contains(S,'FILPS') ; 
    % Replace 
    S(idx) =  {strcat('FILPS= ',char(39),strcat(LCDir,'/',lcmodel_name,'_out','/',lcmodel_name,'.ps'),char(39))} ; 
    % Find the string you want to change 
    idx = contains(S,'FILCOO') ; 
    % Replace 
    S(idx) =  {strcat('FILCOO= ',char(39),strcat(LCDir,'/',lcmodel_name,'_out','/',lcmodel_name,'.coord'),char(39))} ; 

    % Write to text file 
    fid = fopen(strcat(LCDir,'/',lcmodel_name,num2str(cfile_nb),'.control'),'w');
    fprintf(fid,'%s\n',S{:});
    fclose(fid);
end