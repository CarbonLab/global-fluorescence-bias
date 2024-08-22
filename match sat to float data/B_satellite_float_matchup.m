% Program to add satellite data to float data files
% Takes a median accross already extracted grids of satellite data 
% (completed in A_get_satellite_grids_for_SPROF.m) per float profile, 
% OD and Zeu averages are also created
%
% INPUT:
% Float data
%  Processed .mat float files created in
%  "A_format_SPROF_data.m"
%  that have been matched to a 1°x1° GRID of satellite data using:
%  "A_get_satellite_grids_for_SPROF.m" and function "exctract_sat_data.m"
%
% Satellite data
%  Adds satellite chl, optical depth and zeu estimates, and
%  OD-averaged data.
%  Satellite chl is MODIS 8 day, 4km resolution downloaded from NASA OC and
%  processed using "Process_NASA_nc_hdf_files.m", data file located under:
%  MBGC_SatData/Data/Satellite/Downloaded/MODIS/MODIS_8Day....
%
% Functions called:
% get_satchl_to_floatODfluor_GlobalFL2Chl_V4 - match up satellite data to
% each float profile using res (resolution) defined below, using NPQ
% correction from BGC Argo GitHub
%
% OUTPUT:
% float data in 'f' structure, .mat file saved per float with extracted
% area of satellite data included per profile
%
% Created: 05/30/2023 By Jacki Long (MBARI)

clear;close all;clc

%% Opts for get_satchl_to_floatODfluor_GlobalFL2Chl_V3.m function below
% What satellite data product to use?
vnm = 'chl_4km';
tres = '8day'; % daily
% Input grid size, in degrees
grid = 1;
% Name for output data '_[sat product]_[matchup res]_[sat res]'
if grid == 1
    fvnm = ['_ChlOC_1deggrid_' vnm(end-2:end) tres];
    svnm = ['_ChlOC_1deg_' vnm(end-2:end) tres];
elseif grid == 0.5
    fvnm = ['_ChlOC_halfdeggrid_' vnm(end-2:end) tres];
    svnm = ['_ChlOC_halfdeg_' vnm(end-2:end) tres];
else
    disp('Make new fvnm option in if statement')
end
% Resolution for buffer area - grids that will be averaged, in degrees
% These are defined and looped through in the function below
% 1/2, 1/4, 1/6, 1/8 degrees
% 0.5, 0.25, 0.17, 0.125, 0.036
% 55.5, 27.75, 18.5, 13.875, 4km (minimum resolution, use exact match function instead)

%% Create and define your paths
if filesep == '/'
    dpath = ['/Volumes/MBGC_SatData/Data/SPROF/mat_files_' fvnm '/'];
    savepath = ['/Volumes/MBGC_SatData/Data/SPROF/mat_files_' svnm '/'];
else
    dpath = ['X:\Data\SPROF\mat_files_' fvnm '\'];
    savepath = ['X:\Data\SPROF\mat_files_' svnm '\'];
end
% Get list of float .mat files that have been processed
cd(dpath)
float = dir('*_Sprof.mat');
float = {float.name};
float(contains(float,'._')) = [];
% find(contains(float,'2902239')==1) ; %6901511
% 2903451 i = 42, profile 1 at +179.9386°, profile 2 = <0, 
% also last three profiles are after satellite record, and is near 0deg lat

odvars = {'chla_raw','chla_drk_off','chla_npq','no3','doxy','pH','bbp',...
    'chla_adj','temp','sal'};

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
        disp(['Moving to next file without success: All float NPQ chl data is NaN']);
        ee = ee + 1;
    else
        try
            % Add median of satellite grid data and OD data
            disp(['Starting chl satellite extraction: ' datestr(datetime) ' Float ' char(fnum)])
            get_satchl_to_floatODfluor_GlobalFL2Chl_V5
            disp(['Finished chl satellite extraction: ' datestr(datetime)])
            if isnan(mean(f.chla_sat,'omitnan'))
                % Save a log
                problem_file(1,ee) = {fnum};
                problem_file(2,ee) =  {'All sat_chla data is NaN'};
                problem_file(3,ee) =  {[ 'i = ', num2str(i)]};
                problem_file(4,ee) =  {[ 'nprofs', num2str(length(unique(f.profile)))]};
                disp(['Moving to next file without success: All sat_chla data is NaN']);
                ee = ee + 1;
            else
                %Save the float data
                if ~isfolder(savepath)
                    mkdir(savepath)
                else
                end

                plot(f.udate,f.mld(1,:),'k.'); hold on
                plot(f.udate,f.od_sat(1,:),'g.'); hold on
                plot(f.udate,f.zeu_sat(1,:),'b.'); hold on
                datetick('x','mm yyyy')

                save([savepath fnum(1:end-4) '.mat'],'f'); clear f
                disp(['Float ' char(fnum) ' saved as .mat'])
            end

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
        vnm sat_time odvars

end

%Output issues with processing
if exist('problem_file')
    save([savepath 'problem_floats.mat'],'problem_file')
else
end
