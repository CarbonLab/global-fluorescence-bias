% Make globally gridded climatolgies and other statistics
% This _Kd version was made to run using Kd float data which estimated
% Chl from Kd using irradiance
%
% Function called: make_global_gridded_monthly_clim_Kd
%
% INPUT:
% 'm' struct created in D_match_Kd_data.m, see Processing Steps for
% Global Analysis 
%
% OUPUT:
% 'G' structs, representing globally gridded (1deg or X deg res) data
% averaged over DOY, month, or season, and climatological amplitude,
% e.g., Gdoy Gclim Gseas GseasAMP, respectively, represent
% the 1° resolution data and GR* are the regridded data to the optional
% second reolution (likely 5° or 10°, can tell by size of array)
%
% For each calculation above, an individual float is averaged for that lat,
% lon, and day, then any following averages are calculated. This is to
% avoid temporal profiling bias.
%
% Fluor:Chl is calculated on original float data first, before any averaging

% Written By: Jacki Long (MBARI) 2021

clear;close all;clc

nm = 'ChlOC_1deg_4km8day';

% Set paths and load float data:
if filesep == '/'
    dpath = ['/Volumes/MBGC_SatData/Data/Fluor2Chl/mstructs/'];
    savepath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Gstructs/';
    figpath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Figures/';
else
    dpath = ['X:/Data/Fluor2Chl/mstructs/'];
    savepath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Gstructs/';
    figpath = ['X:/Data/Fluor2Chl/Figures/'];
end

load([dpath '/compiled_SPROF_data_' nm 'QCd_Kd.mat'])

%% Begin Program %%

%% STEP 1: FIX UP DATA, define new variables and some minor QC

% Set bad data (QF = 4) to NaN. Data were flagged in "match_Kd_data.m".
% Also cut data less than 30° sun angle
% 6,148 additional profiles remove with low sun angle
fn = fieldnames(mqc_kd);
mkd_in = mqc_kd;
badidx = find(mqc_kd.QF == 4 | mqc_kd.SunAngle < 30); 
for v = 1:length(fn)
    if ~contains(char(fn(v)),'QF') || ~contains(char(fn(v)),'float_ID') || ...
            ~contains(char(fn(v)),'dt') || ~contains(char(fn(v)),'lon') || ~contains(char(fn(v)),'lat')
        mkd_in.(char(fn(v)))(badidx) = NaN;
    else
    end
end


%% STEP 2: Make gridded climatologies

%% STEP 2a: Make calculations at small resolution
% Calculate climatologies at smaller resolution
clearvars -except nm savepath dpath figpath buff_idx mqc_kd mqc mkd_in; close all

res = 5;
[Gclim, Gstd, Gstat, Ginter_var, GclimAMP, GclimSTD, Gseas] = make_global_gridded_monthly_clim_Kd(mkd_in,res);
save([savepath '/G_clim_' nm num2str(res,2) 'deg_Kd.mat'],'Gclim','Gstd','Gstat','GclimAMP','GclimSTD', 'Gseas','mqc','mqc_kd','mkd_in')




