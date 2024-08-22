% Climatologies are calculated as follows:
% Data from 'm' are extracted first, within a lat, lon,
% and month. Then a single year within this index is extracted
% individual floats are averaged first, and then all float averages
% for this month, lat and lon are averaged.
% MEDIAN is used rather than MEAN
%
% OUTPUTS:
% MAIN STRUCTS:
% Gclim = monthly climatologies (using median)
% Gstd = the stdev at each grid, representing the intra-annual variability
% Gstat = multiple structures within for each statistic
% Number of unique years of data that went into the median/std calc
% What years are included in the grid
% Number of unique floats
% What floats went into a grid
% Number of profiles per year
% Total N profiles in that grid
% GclimAMP = monthly maxima - monthly minima in each gridded climatology,
% if there are less than 7 months of data then it is NaN
% GclimSTD = standard deviation of the monthly climatology at each grid,
% if there are less than 7 months of data then it is NaN
% Gseas = seasonal climatolgy 
% Ginter_var = the inter-annual variability of each grid, calcaulted by
% looping through full years (at least 350 days) for each float, and
% calculating the stdev of that year's data after applying a moving mean of
% 5 to the data (smoothing filter). Then the average for all floats and
% years of stdev is taken
% STRUCTURE FIELDS:
% [{'chla_raw_odmn'},{'chla_npq_odmn'},{'chla_drk_off'},{'chla_sat'},...
% {'fluor2chl'},{'rel_diff'},{'diff'}];
%
% INPUTS:
% m = structure with vectors of float data
% res = resolution to grid data

function [Gclim, Gstd, Gstat, Ginter_var, GclimAMP, GclimSTD, Gseas] = make_global_gridded_monthly_clim_Kd(m_buff,res)

% Prep input data from m structure
% Get DOY vector
[yy, mm, dd, doy] = get_doy(m_buff.dt);
% convert lon
lontmp = m_buff.lon;
tt = find((lontmp)>180==1);
lontmp(tt) = lontmp(tt) - 360; clear tt

% Pre-make G arrays
gv = [{'chla_raw_odmn'},{'chla_npq_odmn'},{'chla_drk_off_odmn'},{'chla_sat'},...
    {'fluor2chl'},{'AvgFCHLNPQ'},{'Kd490'},...
    {'AvgKd490JODCHLM07'},{'AvgKd490JZeuCHLM07'},{'AvgKd490JMLDCHLM07'},...
    {'Kd2sat'},{'fluor2Kd'}];
gvv = [{'Gclim'}, {'Gstd'}, {'Gstat'},{'GclimAMP'}, ...
    {'GclimSTD'},{'Gseas'},{'Ginter_var'}];
statvars = [{'Nufl'},{'Nuyr'},{'uyr'},{'ufl'},{'Nobs'}];
for v = 1:length(gv)
    for vv = 1:length(gvv)
        if contains('Gclim',char(gvv(vv))) || contains('Gstd',char(gvv(vv)))
            eval([(char(gvv(vv))) '.' (char(gv(v))) '=  nan([180/res,360/res,12])']);
        elseif contains('Gseas',char(gvv(vv)))
            eval([(char(gvv(vv))) '.' (char(gv(v))) '=  nan([180/res,360/res,4])']);
        elseif contains('Gstat',char(gvv(vv)))
            for j =1:length(statvars)
                if contains('uyr',char(statvars(j))) || contains('ufl',char(statvars(j)))
                    % eval([(char(gvv(vv))) '.' char(statvars(j)) '.' (char(gv(v))) '=  nan([180/res,360/res,12])']);
                    eval([(char(gvv(vv))) '.' char(statvars(j)) '.' (char(gv(v))) '=  cell([180/res,360/res,12])']);
                else
                    eval([(char(gvv(vv))) '.' char(statvars(j)) '.' (char(gv(v))) '=  nan([180/res,360/res,12])']);
                end
            end
        else
            eval([(char(gvv(vv))) '.' (char(gv{v})) '=  nan([180/res,360/res,1])']);
        end
    end
end

for vv = 1:length(gvv)
    eval( ['[' (char(gvv(vv))) '.lon,'  (char(gvv(vv))) '.lat]' '=meshgrid(-180+res/2:res:180-res/2,-90+res/2:res:90-res/2)']);
end

% Loop through months
for h = 1:12
    for i = 1:size(Gclim.lat,1)
        for j = 1:size(Gclim.lon,2)
            clear idx
            idx = find(mm == h & m_buff.lat <= Gclim.lat(i,1) + res/2 & ...
                m_buff.lat > Gclim.lat(i,1) - res/2 & ...
                lontmp <= Gclim.lon(1,j) + res/2 & lontmp > Gclim.lon(1,j) - res/2);
            if ~isempty(idx)
                uyr = unique(yy(idx));
                for v = 1:length(gv)
                    if isnan(mean(m_buff.(char(gv{v}))(idx),'omitnan'))
                        % Skip, all NaN values
                    else
                        clear tmp2
                        for k = 1:length(uyr)
                            uyridx = find(yy(idx) == uyr(k));
                            ufl = unique(m_buff.float_ID(idx(uyridx)));
                            clear tmp
                            for fl = 1:length(ufl)
                                flidx = find(m_buff.float_ID(idx(uyridx)) == ufl(fl));
                                idx2 = idx(uyridx(flidx));
                                % Average same-day data so that high
                                % resolution sampling isn't biasing the data
                                clear data idata; data = m_buff.(char(gv{v}))(idx2)';
                                uday = unique(dd(idx2));
                                if length(uday) ~= length(dd(idx2)) % Average same-day data, if it exists
                                    for jj = 1:length(uday)
                                        idx3 = find(dd(idx2) == uday(jj));
                                        idata(jj) = median(data(idx3),2,'omitnan'); clear idx3
                                    end
                                else
                                    idata = data;
                                end
                                tmp(fl) = median(idata,'omitnan'); clear idata idx2 flidx uday% Average this float's days in this month/year
                            end % Done with unique float loop
                            tmp2(k) = median(tmp,'omitnan'); clear tmp ufl uyridx% Average accross floats for this year
                        end % Done with unique year loop
                        Gstat.Nuyr.(char(gv{v}))(i,j,h) = sum(~isnan(tmp2)); % number of unique years averaged
                        Gstat.uyr.(char(gv{v}))(i,j,h) = {uyr(~isnan(tmp2))}; % unique years averaged
                        Gstat.Nobs.(char(gv{v}))(i,j,h) = sum(~isnan(m_buff.(char(gv{v}))(idx))); % Total data points that went into this median grid

                        % Index of valid data
                        valididx = ~isnan(m_buff.(char(gv{v}))(idx));
                        Gstat.Nufl.(char(gv{v}))(i,j,h) = length(unique(m_buff.float_ID(idx(valididx)))); % number of unique floats in data
                        Gstat.ufl.(char(gv{v}))(i,j,h) = {(unique(m_buff.float_ID(idx(valididx))))}; % WMOs of unique floats in data
                        Gclim.(char(gv{v}))(i,j,h) = median(tmp2,'omitnan'); % Average accross all years for this month
                        Gstd.(char(gv{v}))(i,j,h) = std(tmp2,'omitnan');% Stdev accross all years for this month
                        clear valididx tmp2
                    end
                end % Done with variable loop
            else
            end; clear uyr
        end % Done with longitude loop
    end % Done with latitude loop
end % Done with month loop


%% Get an estimate of how large the annual variability is in fluor:chl at each grid
for i = 1:size(Gclim.lat,1)
    for j = 1:size(Gclim.lon,2)
        clear idx
        idx = find(m_buff.lat <= Gclim.lat(i,1) + res/2 & ...
            m_buff.lat > Gclim.lat(i,1) - res/2 & ...
            lontmp <= Gclim.lon(1,j) + res/2 & lontmp > Gclim.lon(1,j) - res/2);
        if ~isempty(idx)
            for v = 1:length(gv)
                m = 1;
                tmp = NaN;
                % Loop through floats in this grid
                ufl = unique(m_buff.float_ID(idx));
                for fl = 1:length(ufl)
                    % Loop through unique years
                    flidx = find(m_buff.float_ID(idx) == ufl(fl));
                    clear uyr; uyr = unique(yy(idx(flidx)));
                    % Loop through years
                    for k = 1:length(uyr)
                        uyridx = find(yy(idx(flidx)) == uyr(k));
                        if max(doy(flidx(uyridx))) -  min(doy(flidx(uyridx))) >= 340
                            % Full year represented -> estimate variabiltity after applying a moving
                            % mean to data
                            idx2 = idx(flidx(uyridx));
                            clear data; data = m_buff.(char(gv{v}))(idx2)'; clear idx2
                            clear mdata; mdata = movmean(data,5);
                            % Stdev of this year
                            tmp(m) = std(mdata,'omitnan');
                            m = m + 1;
                        else % not a full year
                        end
                    end
                end
                Ginter_var.(char(gv{v}))(i,j) = mean(tmp,'omitnan'); clear tmp;
            end
        else
        end
    end
end

%% Estimate amplitude of monthly climatology & stdev of monthly climatology
for v = 1:length(gv)
    GclimAMP.(char(gv{v})) = max(Gclim.(char(gv{v})),[],3) - min(Gclim.(char(gv{v})),[],3);
    GclimSTD.(char(gv{v}))  = std(Gclim.(char(gv{v})),0,3,'omitnan');
    % If less than 7 months of observations occurred, set the Amp to NaN
    clear tmp; tmp = sum(~isnan(Gclim.(char(gv{v}))),3);
    GclimAMP.(char(gv{v}))(tmp < 7)  = NaN;
    GclimSTD.(char(gv{v}))(tmp < 7)  = NaN;
end


%% Seasonal climatology
[~,mm,~,~] =  get_doy([1:365]);
ss = [12,1,2;3,4,5;6,7,8;9,10,11]; % Months per season in N. Hemi
ss_south =  [6,7,8;9,10,11;12,1,2;3,4,5];% Months per season in S. Hemi

for i = 1:4
    for v = 1:length(gv)
        nidx = find(Gclim.lat(:,1) > 0);
        sidx = find(Gclim.lat(:,1) <= 0);
        Gseas.(char(gv{v}))(nidx,:,i)  = median(Gclim.(char(gv{v}))(nidx,:,ss(i,:)),3,'omitnan');
        Gseas.(char(gv{v}))(sidx,:,i)  = median(Gclim.(char(gv{v}))(sidx,:,ss_south(i,:)),3,'omitnan');
        Gseas.Nobs.(char(gv{v}))(nidx,:,i)  = sum(Gstat.Nobs.(char(gv{v}))(nidx,:,ss(i,:)),3,'omitnan');
        Gseas.Nobs.(char(gv{v}))(sidx,:,i)  = sum(Gstat.Nobs.(char(gv{v}))(sidx,:,ss_south(i,:)),3,'omitnan');
        for k = 1:length(statvars)
            % if contains('Nuyr',char(statvars(k))) || contains('Nufl',char(statvars(k)))
            %     Gseas.Gstat.(char(statvars{k})).(char(gv{v}))(nidx,:,i)  = median(Gstat.(char(statvars{k})).(char(gv{v}))(nidx,:,ss(i,:)),3,'omitnan');
            %     Gseas.Gstat.(char(statvars{k})).(char(gv{v}))(sidx,:,i)  = median(Gstat.(char(statvars{k})).(char(gv{v}))(sidx,:,ss_south(i,:)),3,'omitnan');
            % else % cell array
            %     Gseas.Gstat.(char(statvars{k})).(char(gv{v}))(nidx,:,i)  = median(Gstat.(char(gv{v}))(nidx,:,ss(i,:)),3,'omitnan');
            %     Gseas.Gstat.(char(statvars{k})).(char(gv{v}))(sidx,:,i)  = median(Gstat.(char(gv{v}))(sidx,:,ss_south(i,:)),3,'omitnan');
            % end
        end
    end
end



return