% Program to match up Kd data with my already processed float data
% outputs structure 'mqc_kd' with original 'mqc' too

clear;close all;clc

nm = 'ChlOC_1deg_4km8day';

% Set paths and load float data:
if filesep == '/'
    dpath = ['/Volumes/MBGC_SatData/Data/Fluor2Chl/mstructs/'];
    figpath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Figures/';
else
    dpath = ['X:/Data/Fluor2Chl/mstructs/'];
    figpath = ['X:/Data/Fluor2Chl/Figures/'];
end
% Load my float data
load([dpath '/compiled_SPROF_data_' nm 'QCd.mat'])
% New path to updated data:
% jpath = '/Volumes/Chem/ARGO/Kd490_Fluor_Chl/data/';
jpath = '/Volumes/MBGC_SatData/Data/Fluor2Chl/Josh data Oct 10/josh_Kd490_chl_101023/';
load([jpath 'kd_vs_fluor_LRdP2m_4JL.mat'])

filter = 1; % Option to filter the Kd data
% Options to filter out data
% Filter out specific floats idenfitied in ODV
if filter == 1
    OGIRvsFL = IRvsFL;
    % Make a QF field
    IRvsFL.hdr(:,end+1) = {'QF'};
    IRvsFL.data(:,end+1) = 2.*ones(length(IRvsFL.data),1);
    clear QFidx; QFidx = find(matches(IRvsFL.hdr,'QF')==1);
    T = readtable('/Volumes/MBGC_SatData/Data/Fluor2Chl/floats_to_ignore.xlsx');
    for i = 1:length(T.FloatWMO)
        clear idx; idx = find(IRvsFL.data(:,1) == T.FloatWMO(i));
        IRvsFL.data(idx,QFidx) = 4;
    end
    % Ignore a regional boundary above 50deg latitude cuz CDOM/Case 2
    %     clear idx; idx = find(IRvsFL.data(:,5) > 50);
    %     IRvsFL.data(idx,:) = NaN;
    % Apply Josh's suggestions to filter data further
    clear varidx; varidx = find(matches(IRvsFL.hdr,'RR490')==1);
    % RR fit less than 80% -> might be low light or high cloud
    clear badidx; badidx = find(IRvsFL.data(:,varidx) <= 0.8);
    clear varidx; varidx = find(matches(IRvsFL.hdr,'Kd490 rel std %')==1);
    % Only use data with st dev < 10
    badidx = [badidx; find(IRvsFL.data(:,varidx) >= 10)];
    badidx = unique(badidx);
    IRvsFL.data(badidx,QFidx) = 4;
else
end
% Make mat structure instead of using table
fn = IRvsFL.hdr;
for i = 1:length(fn)
    clear cvar; cvar = char(fn(i));
    cvar = cvar(isstrprop(cvar,'alphanum')); % Only extract letters
    F.(cvar) = IRvsFL.data(:,i);
end
% Convert lon
F.Lon(find(F.Lon>180)) = F.Lon(F.Lon>180) - 360;

save([jpath '/Josh_Kd_only_data_QF.mat'],'F')


clear IRvsFL log_file filter fn badidx cvar idx list night_ct profile_ct varidx i T

%% Match Kd data to existing data
% Variables you want from new 'F' structure
fn = fieldnames(F);
uf = unique(F.WMO); % 290 floats total, 49 floats are not found in my data, maybe some floats had IRR and no FL?

% Variables to include from my dataset in mqc
OGvars = {'chla_raw_odmn','chla_drk_off_odmn','chla_npq_odmn','chla_sat',...
    'od_sat','zeu_sat','chla_std','nsat_pixels','nsat_pixels_valid',...
    'fluor2chl'};

% Loop through fieldnames to pre-size arrays - this avoids adding zeros if
% just looping through
for j = 1:length(fn)
    cvar = char(fn(j));
    mqc_kd.(cvar) = nan(size(mqc.lat));
end
for j = 1:length(OGvars)
    cvar = char(OGvars(j));
    mqc_kd.(cvar) = nan(size(mqc.lat));
end

% Loop through unique floats and put Josh's data
for i = 1:length(uf)
    clear idx_kd; idx_kd = find(F.WMO == uf(i)); % Where you'll pull data from

    % See if this float is in mqc structure
    clear tmpidx;tmpidx = find(mqc.float_ID == uf(i));
    if ~isempty(tmpidx)
        for k = 1:length(idx_kd)
            % Find index in mqc where you'll put the data
            clear idx_mqc;idx_mqc = find(mqc.float_ID == uf(i) & mqc.dt == F.SDN(idx_kd(k)));
            % Loop through fieldnames and add data
            for j = 1:length(fn)
                cvar = char(fn(j));
                mqc_kd.(cvar)(idx_mqc) = F.(cvar)(idx_kd(k));
            end
            % Loop through my variables to save matching data set in time
            for j = 1:length(OGvars)
                cvar = char(OGvars(j));
                mqc_kd.(cvar)(idx_mqc,1) = mqc.(cvar)(idx_mqc,5);
            end        
        end
    else
        % disp(uf(i)) % No. floats
    end

end
% Add standard vars
mqc_kd.dt = mqc.dt;
mqc_kd.lat = mqc.lat;
mqc_kd.lon = mqc.lon;
mqc_kd.lon_W = mqc.lon_W;
mqc_kd.float_ID = mqc.float_ID;
mqc_kd.profile = mqc.profile;

% Define ratios
mqc_kd.fluor2Kd = mqc_kd.chla_npq_odmn./mqc_kd.AvgKd490JODCHLM07;
mqc_kd.Kd2sat = mqc_kd.AvgKd490JODCHLM07./mqc_kd.chla_sat;

save([dpath '/compiled_SPROF_data_' nm 'QCd_Kd.mat'],'mqc','mqc_kd')


