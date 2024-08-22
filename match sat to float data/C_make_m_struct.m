% Loop through all available float data and concatnate them all together
% Logs an estimate of the gain offset
%
% This version loops through the multiple buffer regions and outputs them
% into one m structure

clear;close all;clc



%% Create and define your paths
vnm = 'chl_4km';
tres = '8day'; % daily
nm = ['_ChlOC_1deg_' vnm(end-2:end) tres];

if filesep == '/'
    dpath = ['/Volumes/MBGC_SatData/Data/SPROF/mat_files_' nm '/'];
    savepath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/mstructs/';
else
    dpath = ['X:\Data\SPROF\mat_files_' nm '\'];
    savepath = 'X:/Data/Fluor2Chl/mstructs/';
end

% Variables from float data to output (must be vector data, no depth-data)
% Horizon depths averaged over 
depths = [{'zeumn'},{'mldmn'},{'odmn'}];
avg_vars = {'chla_raw','chla_drk_off','chla_npq','chla_adj','temp','sal',...
    'no3','doxy','pH','bbp'};
gv = [];
for i = 1:length(avg_vars)
    for j = 1:length(depths)
        gv = [gv, {[char(avg_vars(i)), '_', char(depths(j))]}];
    end
end
% Add other variables
gv = [gv, {'lat'},{'lon'},{'lon_W'},{'profile'},{'dt'},{'fl_sat_time_diff'},...
    {'od_sat'},{'zeu_sat'},{'mld'},{'chla_raw_surf'},{'chla_drk_off_surf'},...
    {'chla_npq_surf'},{'surf_depth'}...
    {'chla_sat'},{'chla_std'},{'nsat_pixels'},{'nsat_pixels_valid'}];

% Error count
err = [];
% Set up m vectors (length of all appended float data)
for v = 1:length(gv)
    m.(char(gv{v}))  = [];
end
m.float_ID = [];

cd(dpath)
flist = dir('*Sprof.mat');
float = {flist.name};
float(contains(float,'._')) = [];
% float = {'5904685_Sprof.mat'};

mm = 1; % Trouble shooting loop
no_npq_fl = {};
for i = 1:length(float)
    load([dpath '/' char(float(i))])

    % Reshape some fields to consistent dimensions
    [~,pidx] = unique(f.profile);
    f.dt = f.date(pidx);
    f.lat = f.lat(pidx);
    f.lon = f.lon(pidx);
    f.lon_W = f.lon_W(pidx);
    f.profile = f.profile(pidx);
    fn = fieldnames(f);
    for j = 1:length(fn)
        cvar = char(fn(j));
        cdata = f.(cvar);
        % If a vector and row instead of column, switch orientation
        if sum(size(cdata) > 1 ) == 1 && size(cdata,1) < size(cdata,2)
            f.(cvar) = f.(cvar)';
            % If array and profiles are on row dimension, switch to column
        elseif sum(size(cdata) > 1 ) == 2 && size(cdata,1) == length(f.sat_buffer_km)
            f.(cvar) = f.(cvar)';
        end
    end

    if i == 1
        % Save satellite buffer distance info
        m.sat_buffer_km = f.sat_buffer_km;
    else
    end

    % Some floats have no valid chl_npq_odmn data, but all floats have
    % valid chl sat data
    if ~isnan(mean(f.chla_npq_odmn,'all','omitnan'))
    else
        no_npq_fl(mm,1)  = float(i);
        no_npq_fl(mm,2)  = {'no ~isnan chl data'};
        mm = mm + 1;
    end

    clear tmp; tmp = char(float(i));
    tmp = str2num(tmp(1:7));

    m.float_ID = [m.float_ID;repmat(tmp,length(f.lon),1)];

    for v = 1:length(gv)
        m.(char(gv{v}))  = [m.(char(gv{v}));f.(char(gv{v}))];
    end

    clear f upro pidx yy doy
    if i == length(float)/4
        disp('1/4 done')
        disp(datetime)
    elseif i == length(float)/2
        disp('1/2 done')
        disp(datetime)
    end
end

% Show floats not included
disp(['These floats had no valid NPQ_chla data: ' ])
disp([no_npq_fl])

% Get local time values and sun angle
for i = 1:length(m.lat)
    gps = [m.dt(i), m.lon_W(i), m.lat(i)];
    [Az,El] = SolarAzElq(gps(:,1),gps(:,3),gps(:,2),0); % sdn lat lon alt
    m.Az(i,1) = Az;
    m.El(i,1) = El;
    if El < 0 %Then it's at night
        m.night(i,1) = 1;
    elseif El >= 0
        m.night(i,1) = 0;
    end
    clear El Az
end
disp(['Finished making m vector! '])

save([savepath 'compiled_SPROF_data' nm '.mat'],'m','err')




