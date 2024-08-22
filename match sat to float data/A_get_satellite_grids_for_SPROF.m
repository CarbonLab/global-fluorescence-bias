% Program to add satellite data to float data files
%
% INPUT:
% Float data
%  Processed .mat float files created in
%  "A_format_SPROF_data.m"
% Satellite data
%  Extracts satellite chl created in B_NASAOC_nc2mat.m over a grid area 
%  (multiple pixels) for each float profile
%
%  Satellite chl is MODIS 8 day, 4km resolution downloaded from NASA OC
%  Data file located under:
%  MBGC_SatData/Data/Satellite/Downloaded/MODIS/MODIS_8Day....
%
% Functions called:
% extract_sat_data - match up satellite data to
% each float profile using res (resolution in degrees) defined below
%
% Created: 02/06/2023 By Jacki Long (MBARI)

clear;close all;clc

%% Put the date that you processed last here, any new profiles since will be added
% last_download_date = datenum('January 01 2002');

%% Opts for get_satchl_to_floatODfluor_GlobalFL2Chl_V3.m function below
% What satellite data product to use?
vnm = 'chl_4km';
tres = '8day'; % daily
% Resolution buffer for match-ups, in degrees
res = 1;
% Name for output data '_[sat product]_[matchup res]_[sat res]'
if res == 1
    fvnm = ['_ChlOC_1deggrid_' vnm(end-2:end) tres];
elseif res == 0.5
    fvnm = ['_ChlOC_halfdeggrid_' vnm(end-2:end) tres];
else
    disp('Make new fvnm option in if statement')
end


%% Create and define your paths
if filesep == '/'
    dpath = ['/Volumes/MBGC_SatData/Data/SPROF/mat_files/'];
    savepath = ['/Volumes/MBGC_SatData/Data/SPROF/mat_files_' fvnm '/'];
else
    dpath = 'X:\Data\SPROF\mat_files\';
    savepath = ['X:\Data\SPROF\mat_files_' fvnm '\'];
end
% Get list of float .mat files that have been processed
cd(dpath)
float = dir('*_Sprof.mat');
float = {float.name};
float(contains(float,'._')) = [];

% Grab satellite data for variable "vnm"
if filesep == '/'
    satpath = ['/Volumes/MBGC_SatData/Data/Satellite/Downloaded/MODIS/MODIS_' tres];
else
    satpath = ['X:/Data/Satellite/Downloaded/MODIS/MODIS_' tres];
end
clear sat; sat = matfile([satpath '/' vnm '.mat' ],'Writable',false);
% sat.date is the mean of the start time and end time from the satellite
% .nc file
sat_time = sat.date;
sat_lat = sat.lat;
sat_lon = sat.lon;
% Switch to positive lon (0 to 360)
sat_lonE = sat.lon;
sat_lonE(sat_lonE<0) = 360 + sat_lonE(sat_lonE<0);
  
chlopts = {'chla_raw','chla_drk_off','chla_npq'};

% Log problems
ee = 1;% problem_file = [];
for i = 1:length(float)
    fnum = char(float(i));
    load([dpath fnum])
    % Only move forward if the float data have chl values not all NaN
    if sum(~isnan(f.chla_npq)) == 0
        % Save a log
        problem_file(1,ee) = {fnum};
        problem_file(2,ee) =  {'All f.chla_npq data is NaN'};
        problem_file(3,ee) =  {[ 'i = ', num2str(i)]};
        problem_file(4,ee) =  {[ 'nprofs', num2str(length(unique(f.profile)))]};
        disp(['Moving to next file without success: All float ' fnum ' NPQ chl data is NaN']);
        ee = ee + 1;
    else
        try
            % Add satellite data and OD data
            disp(['Starting chl satellite extraction: ' datestr(datetime) ' Float ' char(fnum)])
            extract_sat_data
            disp(['Finished chl satellite extraction: ' datestr(datetime)])
            %Save the float data
            if ~isfolder(savepath)
                mkdir(savepath)
            else
            end
            save([savepath fnum(1:end-4) '.mat'],'f'); clear f
            disp(['Float ' char(fnum) ' saved as .mat'])

        catch ME
            disp(['Moving to next file without success: '  ME.message]);
            % Save a log
            problem_file(1,ee) = {fnum};
            problem_file(2,ee) =  {ME.message};
            problem_file(3,ee) =  {[ 'i = ', num2str(i)]};
            problem_file(4,ee) =  {[ 'nprofs', num2str(length(unique(f.profile)))]};
            ee = ee + 1;
        end
    end
    clearvars -except ee problem_file savepath satpath dpath key1 key2 ...
        float flist float_OD outfls last_download_date subdirs sat res ...
        vnm sat_time sat_lat sat_lon chlopts

end

%Output issues with processing
save([savepath 'problem_floats.mat'],'problem_file')
