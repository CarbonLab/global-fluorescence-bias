%% Convert .nc BGC Argo files to .mat with selected variables
%
%
% Version 3 was copied from V2 to seperate the float
% processing from the optical depth averaging and the satellite extraction
% in effort to save some processing times. Previously, if you wanted to
% change the resolution for match-up satellite data you would need to
% re-process all the float data as well, which was very inefficient.
%
% Version (V2) uses a satellite estimate of OD rather than the
% float-estimated OD. I decided this may be better since the OD estimated
% by the float is using the fluoresence estimate, which is biased. 
% And if we  can't estimate OD using the satellite due to no
% satellite data, it doesn't matter since the chla_sat will be NaN as well
% To see all variable options do ncinfo('floatname.nc')
%
% Input:
% SPROF Argo float data downloaded as monthly snapshot from SEANOE:
% https://www.seanoe.org/data/00311/42182/
%  Skips floats that don't have CHLA (raw) data since this was made for the
%  global fluor:chl estimates
%
% Output:
% 'f' structure containing float data. Each variable is a 1xN vector saved
% as a unique field
% All data are QC'd as following (Josh Plant suggestions):
% *_ADJUSTED DATA:
% Good data = QF 1, 2, 5, and 8
% TEMP & PSAL - Using raw if there are not yet ADJ data processed
% RAW DATA:
% CHL
% Good data = QF 1, 5, and 8
% Additional adjustments:
% Chl_raw values > 50 are  set to NaN
% Dark offset applied - finds lowest value of first five casts
% that go deeper than 900 m and applies the median of these as an offset to
% entire profile
% Step 2: NPQ correction - if daytime cast, will find the chl max value
% above the MLD (Xing 2012) and apply through to surface. Chl data are
% filtered to remove spikes prior to searching for max, however spikes are
% added back in if they exist.
% MLD is estimated using a density threshold
% Zeu is estimated using Kd from Kim et al 2015
%
% Functions called:
% get_NPQcorr_currentBGCversion - NPQ correction
% SolarAzElq - Estimate sun angle
% dark_offset - Fluorescence dark offset
%
% Created: 01/9/2020 By Jacki Long (MBARI)

clear;close all;clc

%% Create and define your paths
if filesep == '/'
    dpath = '/Volumes/CarbonLab/Data/Argo_Sprof_Snapshots/May_2023/202305-BgcArgoSprof/dac/';
    % Use SPROF snapshot with doi
    % Local path
%     dpath = '/Users/jlong/Documents/Data/float_data/Argo_Sprof_Snapshots/202305-BgcArgoSprof/dac/';
    savepath = '/Volumes/MBGC_SatData/Data/SPROF/mat_files/';
else
%     dpath = 'Y:\SPROF\';
    % Use SPROF snapshot with doi
    dpath = 'S:\Data\Argo_Sprof_Snapshots\May_2023\202305-BgcArgoSprof\dac\';
    savepath = 'X:\Data\SPROF\mat_files\';
end

% subdirectories under SPROF
    subdirs = {'aoml','bodc','coriolis','csio','csiro','incois','jma','kma',...
        'kordi','meds'};

% Get list of current float .mat files that have been processed
cd(savepath)
outfls = dir('*.mat');
outfls = {outfls.name};


% Variable names of float data under sprof
key1 = {'TEMP','TEMP_QC','TEMP_ADJUSTED','TEMP_ADJUSTED_QC','JULD',...
    'CYCLE_NUMBER','LONGITUDE','LATITUDE','PRES','PRES_ADJUSTED',...
    'PRES_QC','PRES_ADJUSTED_QC','PSAL',...
    'PSAL_QC','PSAL_ADJUSTED','PSAL_ADJUSTED_QC',...
    'CHLA','CHLA_QC','CHLA_ADJUSTED','CHLA_ADJUSTED_QC','POSITION_QC',...
    'NITRATE_ADJUSTED','NITRATE_ADJUSTED_QC','DOXY_ADJUSTED',...
    'DOXY_ADJUSTED_QC','BBP700','BBP700_QC',...
    'PH_IN_INSITU_TOTAL_ADJUSTED','PH_IN_INSITU_TOTAL_ADJUSTED_QC',...
    'DIRECTION'};
% Variable names we will use for the above names
key2 = {'temp','temp_QC','temp_adj','temp_adj_QC','date',...
    'profile','lon','lat','press','press_adj',...
    'press_QC','press_adj_QC','sal',...
    'sal_QC','sal_adj','sal_adj_QC',...
    'chla_raw','chla_raw_QC','chla_adj','chla_adj_QC','position_QC',...
    'no3','no3_QC','doxy','doxy_QC','bbp','bbp_QC','pHin','pHin_QC',...
    'direction'};


%% Part 1: Get data
notunique_log = [];

ee = 1; % Iteration for problem files
cd(dpath)
% Loop through subdirectories in SPROF
count = [];
latvec = [];
lonvec = [];
for sub = 1:length(subdirs)
    cd([dpath char(subdirs(sub))])
    flist = dir('*.nc');
    float = {flist.name};
    float(contains(float,'._')) = [];
    DIR = subdirs(sub);
    disp(['Starting directory ' char(DIR)])
    for i = 1:length(float) 
        fnum = char(float(i));
        tt = ncinfo(fnum);
        vars = {tt.Variables.Name};
        m = 1;
        zoo = 1;
        % Convert Julian date to matlab date time
        clear x;x = find(strcmp(char(vars),{'REFERENCE_DATE_TIME'})==1);
        tmp = (ncread(fnum,char(vars(x))));
        tmp = tmp';
        startdt = tmp(1:8);
        dn = datenum(startdt,'yyyymmdd'); %Get the start datenum
        tmpdate = ncread(fnum,'JULD');
        tmpdate = dn + tmpdate; %Add the "days from start date" to datenum startdate
        % Must have a CHLA field, otherwise don't bother
        if contains('CHLA',vars)
                % Continue only if the chla data is not all NaN
                if ~isnan(mean(mean(ncread(fnum,'CHLA'),'omitnan'),'omitnan'))
                    %Get array size first by looking at salinity data array
                    clear x;x = find(strcmp(char(vars),{'PSAL'})==1);
                    clear tmp; tmp = ncread(fnum,char(vars(x)));
                    array_size = size(tmp); % array_size(1) = samples, array_size(2) = profiles
                    % Profile counter
                    cnt1 = array_size(2); % length of profiles

                    clear x;x = find(strcmp(char(vars),{'PROJECT_NAME'})==1);
                    Project = (ncread(fnum,char(vars(x))));
                    Project = {Project(:,1)'};
                    clear x;x = find(strcmp(char(vars),{'DATA_CENTRE'})==1);
                    DAC = (ncread(fnum,char(vars(x)))); 
                    DAC = {DAC(:,1)'};
                    WMO = {fnum};
                    start_dt = tmpdate(1);
                    end_dt = tmpdate(end);
                    nprofs = array_size(2);
                    nrow = [nprofs,WMO,Project,DAC,DIR,start_dt,end_dt,1,{'wCHLA'}];
                    if isempty(count)
                        count = table(nprofs,WMO,Project,DAC,DIR,start_dt,end_dt,1,{'wCHLA'});
                    else
                        count = [count;nrow];
                    end

                    if array_size(2) > 4 %Else only 4 profiles collected from this float - skip it
                        clear ff; disp(['Starting Float ' char(fnum) ' in ' char(DIR) ', ' char(datetime)])
                        
                        for j = 1:length(vars)
                            clear x;x = find(strcmp(char(vars(j)),key1)==1); %Find variables that the current float has
                            %If it found a matching variable, then save
                            if ~isempty(x)
                                clear current_var; current_var = ncread(fnum,char(vars(j)));
                                if size(current_var,2) > 1 %Then it is an array of data - reformat it to be one column
                                    %Flip data so it's going from bottom to surface
                                    current_var = flipud(current_var);
                                    current_var = reshape(current_var,[(array_size(1)*array_size(2)),1]);
                                else % It's a single value per profile, need to make it repeat through depth
                                    current_var = current_var'; %First transform into a single row
                                    current_var = repmat(current_var, array_size(1),1);%Then repete the data through depth
                                    current_var = reshape(current_var,[(array_size(1)*array_size(2)),1]);
                                end
                                vars_out(m) = key2(x);
                                vname(m) = key1(x);
                                eval(['ff.' char(vars_out(m)) '= current_var;']);
                                m = m+1;
                            else % If did not find, output as all NANs
                            end
                            clear current_var
                        end
                        % Make QF fields numerical
                        % Setting empty rows to QF = 9
                        vv = {'temp_QC','temp_adj_QC','sal_QC','sal_adj_QC',...
                            'chla_raw_QC','chla_adj_QC','press_adj_QC',...
                            'press_QC','position_QC','no3_QC','doxy_QC',...
                            'bbp_QC','pHin_QC'};
                        for v = 1:length(vv)
                            cvar = char(vv(v));
                            if isfield(ff,cvar)
                                clear tmp; tmp = find(ff.(cvar)==' ');
                                ff.(cvar)(tmp) = '9';
                                ff.(cvar) = str2num(ff.(cvar));
                            else
                            end
                        end

                        %% Housekeeping
                        % Remove data where lat, lon, press, or date = NaN, and
                        % where position quality flag is good (=1) and
                        % where it is only ascending profiles (this is to
                        % try to avoid the issue with same cycle number
                        % given to ascending and descending data)
                        keepdata = find(~isnan(ff.lat) & ...
                            ~isnan(ff.lon) & ~isnan(ff.date) & ...
                            ~isnan(ff.press) & ff.position_QC == 1 & ...
                            (ff.direction == 'A'));
                       
                        % Log data removed/set to NaN
                        cnt = 1;
                        data_bad_count(cnt) = length(ff.profile) - length(keepdata);
                        
                        % Loop through variables and only keep good data
                        fvars = fieldnames(ff);
                        for k = 1:length(fvars)
                            var = string(fvars(k));
                            ff.(var) = ff.(var)(keepdata);
                        end

                        %% Pre-define some f parameters
                        f.lon = ff.lon;
                        f.lat = ff.lat;
                        f.profile = ff.profile;
                        f.press = ff.press_adj;
                        f.chla_raw = ff.chla_raw;
                        f.temp = ff.temp_adj;
                        f.sal = ff.sal_adj;
                        f.chla_adj = ff.chla_adj;
                        f.chla_adj_QC = ff.chla_adj_QC;

                        % Convert Julian date to matlab date time
                        clear x;x = find(strcmp(char(vars),{'REFERENCE_DATE_TIME'})==1);
                        tmp = (ncread(fnum,char(vars(x))));
                        tmp = tmp';
                        startdt = tmp(1:8);
                        dn = datenum(startdt,'yyyymmdd'); %Get the start datenum
                        f.date = dn + ff.date; %Add the "days from start date" to datenum startdate

                        % Convert longitude to degrees East if any are negative
                        clear tmp; tmp = find(ff.lon<0);
                        f.lon(tmp) = ff.lon(tmp) + 360;
                        % Keep a degrees West version (-180 to 180)
                        f.lon_W = f.lon;
                        % clear tmp; tmp = find(ff.lon > 180);
                        clear tmp; tmp = find(f.lon > 180);
                        f.lon_W(tmp) = f.lon(tmp) - 360;

                        % Good profile log
                        % Number of profiles that aren't all NaNs
                        clear nprofs; nprofs = sum(~isnan(mean(ncread(fnum,'CHLA'),'omitnan')));
                        clear nrow; nrow = [nprofs,WMO,Project,DAC,DIR,start_dt,end_dt,2,{'wCHLAvalid'}];
                        count = [count; nrow];
                      
                        % Good data log
                        good_data(1) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));

                        %% Use Quality Flags to set data to NaN
                        % TEMPERATURE
                        % Set bad flags to NaN, Good data = QF 1, 2, 5, and 8
                        clear nanidx; nanidx = find(ff.temp_adj_QC ~= 1 & ...
                            ff.temp_adj_QC ~= 2 & ff.temp_adj_QC ~= 5 & ...
                            ff.temp_adj_QC ~= 8);
                        % Log data removed/set to NaN
                        cc = sum(~isnan(f.temp));
                        % Set bad QF for adjusted temp data to NaN
                        f.temp(nanidx) = NaN;
                        % Good data log
                        good_data(2) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));
                        % If there's unprocessed (NaN) adj data, set it to the raw data
                        clear idx; idx = find(isnan(ff.temp_adj));
                        if ~isempty(idx)
                            f.temp(idx) = ff.temp(idx); %set to equal raw temp
                            % Set bad flags to NaN for this raw data, Good data = QF <= 3 or =8
                            clear nanidx; nanidx = find(ff.temp_QC(idx) > 3 & ...
                                ff.temp_QC(idx) ~= 8 & ff.temp_QC(idx) ~= 5);
                            % Log data removed/set to NaN
                            cc = sum(~isnan(f.temp));
                            % Set bad QF for raw temp data to NaN
                            f.temp(idx(nanidx)) = NaN;
                        else
                        end
                        % Good data log
                        good_data(3) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));

                        % SALINITY
                        % Set bad flags to NaN, Good data = QF 1, 2, 5, and 8
                        clear nanidx; nanidx = find(ff.sal_adj_QC ~= 1 & ...
                            ff.sal_adj_QC ~= 2 & ff.sal_adj_QC ~= 5 & ...
                            ff.sal_adj_QC ~= 8);
                        f.sal(nanidx) = NaN;
                        % Good data log
                        good_data(5) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));
                        % If there's unprocessed adj data, set it to the raw data
                        clear idx; idx = find(isnan(ff.sal_adj));
                        if ~isempty(idx)
                            f.sal(idx) = ff.sal(idx);
                            % Set bad flags to NaN for this *_ADJ data, Good data = QF <= 3 or =8
                            clear nanidx; nanidx = find(ff.sal_QC(idx) > 3 & ...
                                ff.sal_QC(idx) ~= 8 & ff.sal_QC(idx) ~= 5);
                            f.sal(idx(nanidx)) = NaN;
                        else
                        end

                        % PRESSURE
                        % Set bad flags to NaN, Good data = QF 1, 2, 5, and 8
                        clear nanidx; nanidx = find(ff.press_adj_QC ~= 1 & ...
                            ff.press_adj_QC ~= 2 & ff.press_adj_QC ~= 5 & ...
                            ff.press_adj_QC ~= 8);
                        f.press(nanidx) = NaN;
                        % Good data log
                        good_data(5) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));
                        % If there's unprocessed adj data, set it to the raw data
                        clear idx; idx = find(isnan(ff.press_adj));
                        if ~isempty(idx)
                            f.press(idx) = ff.press(idx);
                            % Set bad flags to NaN for this *_ADJ data, Good data = QF <= 3 or =8
                            clear nanidx; nanidx = find(ff.press_QC(idx) > 3 & ...
                                ff.press_QC(idx) ~= 8 & ff.press_QC(idx) ~= 5);
                            f.press(idx(nanidx)) = NaN;
                        else
                        end

                        % Good data log
                        good_data(4) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));

                         % BBP
                        % Set bad flags to NaN, Good data = QF 1, 2, 3, 5, and 8
                        if isfield(ff,'bbp')
                            f.bbp = ff.bbp;
                            f.bbp_QC = ff.bbp_QC;
                            clear nanidx; nanidx = find(ff.press_QC(idx) > 3 & ...
                                ff.press_QC(idx) ~= 8 & ff.press_QC(idx) ~= 5);
                            f.bbp(nanidx) = NaN;
                        else
                            f.bbp = NaN(size(ff.sal));
                            f.bbp_QC = NaN(size(ff.sal));
                        end

                        % pH
                        % Set bad flags to NaN, Good data = QF 1, 2, 5, and 8
                        if isfield(ff,'pHin')
                            f.pHin = ff.pHin;
                            f.pHin_QC = ff.pHin_QC;
                            clear nanidx; nanidx = find(ff.pHin_QC ~= 1 & ...
                                ff.pHin_QC ~= 2 & ff.pHin_QC ~= 5 & ...
                                ff.pHin_QC ~= 8);
                            f.pHin(nanidx) = NaN;
                        else
                            f.pHin = NaN(size(ff.sal));
                            f.pHin_QC = NaN(size(ff.sal));
                        end

                        % NO3
                        % Set bad flags to NaN, Good data = QF 1, 2, 5, and 8
                        if isfield(ff,'no3')
                            f.no3 = ff.no3;
                            f.no3_QC = ff.no3_QC;
                            clear nanidx; nanidx = find(ff.no3_QC ~= 1 & ...
                                ff.no3_QC ~= 2 & ff.no3_QC ~= 5 & ...
                                ff.no3_QC ~= 8);
                            f.no3(nanidx) = NaN;
                        else
                            f.no3 = NaN(size(ff.sal));
                            f.no3_QC = NaN(size(ff.sal));
                        end

                        % OXY
                        % Set bad flags to NaN, Good data = QF 1, 2, 5, and 8
                        if isfield(ff,'doxy')
                            f.doxy = ff.doxy;
                            f.doxy_QC = ff.doxy_QC;
                            clear nanidx; nanidx = find(ff.doxy_QC ~= 1 & ...
                                ff.doxy_QC ~= 2 & ff.doxy_QC ~= 5 & ...
                                ff.doxy_QC ~= 8);
                            f.doxy(nanidx) = NaN;
                        else
                            f.doxy = NaN(size(ff.sal));
                            f.doxy_QC = NaN(size(ff.sal));
                        end

                        % CHLA RAW
                        clear nanidx; nanidx = find(ff.chla_adj_QC ~= 1 & ...
                            ff.chla_adj_QC ~= 8 & ff.chla_adj_QC ~= 5);
                        f.chla_raw(nanidx) = NaN;
                        
                        % There's some very large values that need to be
                        % set to NaN
                        f.chla_raw(f.chla_raw > 50) = NaN;
                      
                        % Good data log
                        good_data(6) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw));

                        clear ff
                        %% Start Adjusting Chla_raw data
                        % Dark offset
                        dark_offset
                        % If "No deep profiles available for in situ CHL dark estimate"
                        % it still moves to chl extraction but CHL_Dark = NaN,
                        % so all chla_drk_off data = NaN, then all chl data
                        % = NaN
                        f.chla_drk_off = f.chla_raw - CHL_Dark;

                        % NPQ correction
                        f.chla_npq = f.chla_drk_off;
                        up = unique(f.profile);
                        % MLD and Zeu estimated in NPQ function below
                        f.Dmld = NaN(1,length(up)); % MLD estimated in NPQ function below
                        f.Zeu = NaN(1,length(up)); % Kim et al 2015
                        % Good data log
                        good_data(7) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw) & ~isnan(f.chla_npq));
                        dirs = [];
                        for p = 1:length(up)
                            pidx = find(f.profile == up(p));
                            gps = [f.date(pidx(1)), f.lon_W(pidx(1)), f.lat(pidx(1))];
                            [Az,El] = SolarAzElq(gps(:,1),gps(:,3),gps(:,2),0); % sdn lat lon (in E or W) alt
                            % Night time cast, option to run this section
                            % by uncommenting out below
                            if El < 0 
                                % Not doing MLD estimate here anymore,
                                % switched to estimated in next step so
                                % that all profiles are done exactly the
                                % same and all depth horizons are
                                % calculated in the same step
%                                 % MLD estimate - copied from section in
%                                 % get_NPQcorr_currentBGCversion.m -
%                                 % FIND DENSITY MLD
%                                 mld_den_threshold = 0.03; % Dong et al 2008
%                                 den = density(f.sal(pidx),f.temp(pidx));
% 
%                                 % There will be at least one data point <=25m
%                                 % in original code the data is going from
%                                 % shallow to deep, so I switched 'first' to
%                                 % 'last' here to get the shallow reference - Jacki Long
%                                 goodidx = ~isnan(den); % Some profiles might have NaN salinity vals (e.g., WMO 1902303)
%                                 ref_ind  = find(f.press(pidx(goodidx))<=25, 2, 'last'); % find ref points for mld calc
%                            
%                                 tDMLD = den - mean(den(ref_ind),'omitnan') <= mld_den_threshold; % mixed layer
% 
%                                 if all(tDMLD ==0) || all(tDMLD ==1) || sum(~isnan(den)) == 0
%                                     % f.Dmld(p,1) = NaN; - already
%                                     % pre-defined
%                                 else
%                                     f.Dmld(p,1) = max(f.press(pidx(tDMLD))); % depth of mixed layer
%                                 end
%                                 % f.Zeu(p,1) = NaN;

                                clear tDMLD mld_den_threshold den

                            else % day-time cast, apply npq correction
                                data = [f.press(pidx),f.temp(pidx),f.sal(pidx),f.chla_drk_off(pidx)];
                                if sum(isnan(mean(data,1,'omitnan'))) == 0
                                    clear NPQ; NPQ = get_NPQcorr_currentBGCversion(gps, data, dirs);
                                    if ~isempty(NPQ.data)
                                        t_nan = isnan(f.chla_npq(pidx)); % Keep log of NaN indeces
                                        iXing   = find(strcmp('Xing_MLD',NPQ.hdr) == 1); % chl NPQ based on MLD
                                        iSPIKE  = find(strcmp('CHLspike',NPQ.hdr) == 1); % The filtered chla data - orginal
                                        tNPQ = f.press(pidx) <= NPQ.XMLDZ;% find data above the cmax found in ML
                                        % Set this depth range of data to the corrected NPQ data but add back
                                        % in the spikes
                                        f.chla_npq(pidx(tNPQ & ~t_nan)) = ...
                                            NPQ.data(tNPQ & ~t_nan,iXing) + ...
                                            NPQ.data(tNPQ & ~t_nan,iSPIKE);
                                        f.Dmld(p,1) = NPQ.Dmld; % Density threshold
                                    else % NPQ correction failed, save profile as NaNs
                                        f.chla_npq(pidx) = NaN;
                                    end
                                else % One input variable is all NaN
                                    f.chla_npq(pidx) = NaN;
                                end
                            end
                            clear data gps pidx
                        end
                        % Good data log
                        good_data(8) = sum(~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw) & ~isnan(f.chla_npq));

                        % Profiles with a valid chl data point above 10 m
                        surfidx = find(f.press < 10 & ~isnan(f.press) & ~isnan(f.temp) & ...
                            ~isnan(f.sal) & ~isnan(f.chla_raw) & ~isnan(f.chla_npq));
                        % What profiles are these?
                        clear profs; [profs,pidx] = unique(f.profile(surfidx));
                        clear nprofs; nprofs = length(profs);
                        clear nrow; nrow = [nprofs,WMO,Project,DAC,DIR,start_dt,end_dt,3,{'wCHLAsurf'}];
                        count = [count; nrow];

                        % Map the valid data
                        latvec = [latvec;f.lat(surfidx(pidx))];
                        lonvec = [lonvec;f.lon_W(surfidx(pidx))];
                       
                        % Save data as .mat file with structure 'f'
                        save([savepath fnum(1:end-3) '.mat'], 'f')

                        if length(unique(f.date)) ~= length(unique(f.profile))
                            disp('LENGTH OF UNIQUE DATES NOT EQUAL TO UNIQUE PROFILES')
                            notunique_log = [notunique_log, str2num(fnum(1:end-9))];
                        else
                        end
                    else% Less than 4 profiles collected
                    end
                else% All chla_raw data is NaN
                end
                clearvars -except notunique_log latvec lonvec count DIR ee problem_file savepath dpath key1 key2 float flist outfls subdirs
        else
        end
    end % Loop for floats
end % Loop for sub dirs

% Save count
save([savepath 'valid_profile_count.mat'], 'count','notunique_log')
writetable(count,[savepath 'valid_profile_count.xlsx'])

den = sum(count.nprofs(count.Var8 == 2)); % Total valid Chl profiles
num = sum(count.nprofs(count.Var8 == 3)); % Total valid Chl profiles with > 10m data
disp(['Percent of profile data retained after QC:' num2str((num/den)*100,3) '%'])

% Map all profile locations
figure(2);hold on
m_proj('Robinson','long',[-180 180],'lat',[-90 90]);hold on
h=m_plot(lonvec,latvec,'r.');
m_gshhs_l('patch',rgb('gray'),'edgecolor','k');hold on
% set(h,'edgecolor','none')
m_grid('box','fancy','tickdir','in','linestyle','none','backcolor',rgb('white'));
title('Valid profiles (wCHLAsurf)')
% set(findall(gcf,'-property','Fontname'),'Fontname',figset.font_type)

exportgraphics(figure(2),[savepath 'valid_profile_count_wCHLAsurf.pdf'],'Resolution',200)

% Map day vs night locations
figure(1)
ylabel('Chl')
title('Blue = NPQ (day time), Black = dark off only (night time)')
exportgraphics(figure(1),[savepath 'above10m_data.pdf'],'Resolution',200)

if exist('problem_file','var')
    save([savepath '/problem_floats.mat'],'problem_file');
else
end


%% If you want to load one of the f structs you just made...
% fnum = '3900521';
% load([savepath '/' fnum '_Sprof.mat'])
%
%% And check it out:
%
% figure(1); hold on
% plot(f.chla,-f.press,'o')
% title('chla')
%
% figure(2); hold on
% plot(f.bp700,-f.press,'o')
% title('bbp')
% xlim([0 1])
% ylim([-2000 0])


% f.udate = unique(f.date);
%
% figure(2)
% plot(f.udate, f.chla_sat,'r.');hold on
% plot(f.udate, f.chla_raw_odmn,'b.');hold on
% plot(f.udate, f.chla_drk_off_odmn,'mx');hold on
% plot(f.udate, f.chla_npq_odmn,'g.');hold on
% datetick('x')
% legend('Sat','chla_raw','chla_drk_off','chla_npq')
% title(['Float ' float])



