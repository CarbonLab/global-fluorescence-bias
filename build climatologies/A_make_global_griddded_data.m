% Make globally gridded climatolgies and other statistics
%
% Function called: make_global_gridded_monthly_clim and filter_float_data
%
% INPUT:
% 'm' struct created in C_make_m_struct.m, see Processing Steps for
% Global Analysis 
%
% OUPUT:
% 'G' structs, representing globally gridded (1deg or X deg res) data
% averaged over DOY, month, or season, and climatological amplitude,
% e.g., Gdoy Gclim Gseas GseasAMP, respectively, represent
% the 1° resolution data and GR* are the regridded data to the optional
% second reolution (likely 5° or 10°, can tell by size of array)
%
% Seasons are ordered as Winter, Spring, Summer, Fall in the GClim struct
%
% For each calculation above, an individual float is averaged for that lat,
% lon, and day, then any following averages are calculated. This is to
% avoid temporal profiling bias. 
%
% DATA QC:
% Minor QC occurred in the previous step with generating 'm'. 
% See previous code for details
%
% Fluor:Chl is calculated on original float profile data first, prior to averaging

% Written By: Jacki Long (MBARI) 2021

clear;close all;clc

%% Choose which mat_files to use (varies in match-up resolution and OC product)
% Will save with this name appended as well
nm = 'ChlOC_1deg_4km8day';

% Choose buffer resolution, options = [32	24	16	12	8	4.63]
buff_idx = 8; % Choosing 8 km based on Nina's threshold analysis

% Set paths and load float data:
if filesep == '/'
    dpath = ['/Volumes/MBGC_SatData/Data/Fluor2Chl/mstructs/'];
    savepath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Gstructs/';
    figpath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Figures/';
else
    dpath = ['X:/Data/Fluor2Chl/mstructs/'];
    savepath = ['X:/Data/Fluor2Chl/Gstructs/'];
    figpath = ['X:/Data/Fluor2Chl/Figures/'];
end
% savepath = '/Users/jlong/Documents/GitHub/Lets-look/global fluor2chl/';
% dpath = '/Volumes/CarbonLab/Data/Float data/Global_Fluor2Chl_related/';
load([dpath 'compiled_SPROF_data_' nm '.mat'])


%% Begin Program

% STEP 1: QC DATA
% RATIO, POC, and Cpytho are also calculated in the function below

% Enter Date of latest satellite data - this may not be necessary anymore
sat_dt_end = datenum('29-March-2023'); 
mqc = filter_float_data(m,sat_dt_end);
% Save the filtered data
save([dpath '/compiled_SPROF_data_' nm 'QCd.mat'],'mqc')
clear m err

%% STEP 2: Make gridded climatologies
% Load mqc data
load([dpath '/compiled_SPROF_data_' nm 'QCd.mat'])

%% Make input data for climatologies
% Select column (buffer region) of data you want to build climatology from
idx = find(mqc.sat_buffer_km == buff_idx);

fn = fieldnames(mqc);
for i = 1:length(fn)
    cvar = char(fn(i));
    if size(mqc.(cvar),2) == 1
        m_in.(cvar) = mqc.(cvar);
    else
        m_in.(cvar) = mqc.(cvar)(:,idx);
    end
end


%% STEP 2a: Make calculations at small resolution
% Calculate climatologies at smaller resolution
clearvars -except nm savepath dpath figpath buff_idx m_in mqc; close all

res = 5;
[Gclim, Gstd, Gstat, Ginter_var, GclimAMP, GclimSTD, Gseas] = make_global_gridded_monthly_clim(m_in,res);
save([savepath '/G_clim_' nm num2str(res,2) 'deg.mat'],'Gclim','Gstd','Gstat','GclimAMP','GclimSTD', 'Gseas','Ginter_var','mqc','m_in')
% load([savepath '/G_clim_' nm num2str(res,2) 'deg.mat'])
% % Save as xls data in vector format - option for PlotlyDASH
% [m, n] = size(Gclim.lat);
% lat = reshape(Gclim.lat,[m*n 1]);
% lon = reshape(Gclim.lon,[m*n 1]);
% for i = 1:12
%     D(:,i) = reshape(Gclim.fluor2chl(:,:,i),[m*n 1]);
% end
% % csvwrite([bpath '/Gstructs/G_clim_' num2str(res,2) 'deg.csv'], [lat,lon,D])
% % csvwrite(['/Users/jlong/Desktop/G_clim_' num2str(res,2) 'deg.csv'], [Gclim.fluor2chl(:,:,6)])


%% STEP 2b: Make calculations at larger resolution if you want...
% Calcualte climatologies at 1deg resolution
clearvars -except nm savepath dpath figpath buff_idx m_in mqc; close all
% load([dpath '/compiled_SPROF_data_' nm 'QCd.mat'])
res = 10;
[Gclim, Gstd, Gstat, Ginter_var, GclimAMP, GclimSTD, Gseas] = make_global_gridded_monthly_clim(m_in,res);
save([savepath '/G_clim_' nm num2str(res,2) 'deg.mat'],'Gclim','Gstd','Gstat','GclimAMP','GclimSTD', 'Gseas','mqc')
% Make Figures
% clc; clearvars -except savepath res
% load([savepath '/G_clim_' num2str(res,2) 'deg.mat'])
% figpath = ['/Users/jlong/Documents/MATLAB/PUBS/2022_fluor2chl/res_' num2str(res) '_FL2SatCHL/'];




