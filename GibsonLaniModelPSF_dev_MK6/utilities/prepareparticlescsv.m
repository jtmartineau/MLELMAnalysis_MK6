% this function prepares the particles.xls. It needs to be converted to a
% csv to be read by srx.

function ret=prepareparticlescsv(analysisfile,localizationstructfile)
load(fullfile(analysisfile,localizationstructfile));
vutaraHeaderInfo = {'image-ID','time-point','cycle','z-step','frame','accum','probe','photon-count','photon-count11',...
        'photon-count12','photon-count21','photon-count22','psfx','psfy','psfz','psf-photon-count','x','y','z','stdev',...
        'amp','background11','background12','background21','background22','maxResidualSlope','chisq','log-likelihood',...
        'llr','accuracy','fiducial','valid','precisionx','precisiony','precisionz'};
vutaraCSVData=[];
for ii=1:numel(localization_struct)
    if isempty(localization_struct(ii).stats_array)
        continue
    else
        number_ROIs=size(localization_struct(ii).stats_array,1);
        image_ID=(ii-1);
        time_points=1;
        cycle=1;
        z_step=1;
        frame=image_ID;
        accum=1;
        for jj=1:number_ROIs
            psfx=localization_struct(ii).stats_array(jj,7);
            psfy=localization_struct(ii).stats_array(jj,8);
            psfz=localization_struct(ii).stats_array(jj,9);
            x0=localization_struct(ii).stats_array(jj,11);
            y0=localization_struct(ii).stats_array(jj,12);
            z0=localization_struct(ii).stats_array(jj,13);
            CRLBx=localization_struct(ii).stats_array(jj,27);
            CRLBy=localization_struct(ii).stats_array(jj,28);
            CRLBz=localization_struct(ii).stats_array(jj,29);
            %custom_confidence1=zeros(number_ROIs,1);
            %custom_confidence2=zeros(number_ROIs,1);
            %custom_confidence3=zeros(number_ROIs,1)
            if (abs(psfx)<1000) && (abs(psfy)<1000) && (abs(psfz)<1000) && (abs(x0)<20000) && (abs(y0)<20000) && (abs(z0)<3000) && (abs(CRLBx)<500) && (abs(CRLBy)<500) && (abs(CRLBz)<1000)
                array_temp=cat(2,image_ID,time_points,cycle,z_step,frame,accum,localization_struct(ii).stats_array(jj,:));
                vutaraCSVData=cat(1,vutaraCSVData,array_temp);
            end
        end
    end
end
vutaraCSVData=cast(vutaraCSVData,'single');
[~,nColVutaraCSV] = size(vutaraHeaderInfo);

% Create String for use with "fprintf" with the number of appropriate columns in "vutaraHeaderInfo"
% Use "for" loop to iterate over size of "vutaraHeaderInfo" to generate string of proper size
nColVutaraCell = cell(1,nColVutaraCSV-1);
for i = 1:nColVutaraCSV-1
    nColVutaraCell{1,i} = '%s,';
end
% Append string to proper termination line call for "fprintf"
nColVutaraCell = [nColVutaraCell{1,:},'%s\n'];

% Choose desired file name and destination folder location for file saving
%[saveFileName,saveDirectory,saveFilterIndex] = uiputfile('.csv','Save CSV File As...','particles.csv');
%fullSavePathName = fullfile(saveDirectory,saveFileName);
%if isequal(saveFileName,0) || isequal(saveDirectory,0) || isequal(saveFilterIndex,0)
%    disp('User Cancelled File Save...')
%else
%    disp(['User saved file as:',fullSavePathName])
%end

% Save current working/home directory to use at end of function
%homeDir = pwd;
% Change to save directory
%cd(saveDirectory);

saveFileName='particles.csv';

% Use "fprintf" to write the appropriate column headers to a CSV file called "particles.csv"
fID = fopen(fullfile(analysisfile,saveFileName),'w');
fprintf(fID,nColVutaraCell,vutaraHeaderInfo{1,:});
fclose(fID);

% Use "dlmwrite" to append ThunderStorm data to the CSV file with appropriate column headers
% Variable "data" is the input data here, passed into the function call
% Use 'roffset' value of '0' to avoid systemtatic output of ',' strings in first row from 'append' flag
dlmwrite(fullfile(analysisfile,saveFileName),vutaraCSVData,'delimiter',',','-append','roffset',0,'coffset',0);

ret=saveFileName;
end