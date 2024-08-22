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
% Gstd = the st dev across all years of data within a month grid
% (one option for intra-annual variability)
% Gstat = multiple structures within for each statistic
% Number of unique years of data that went into the median/std calc
% What years are included in the grid
% Number of unique floats
% What floats went into a grid
% Number of profiles per year
% Total N profiles in that grid
% GclimAMP = monthly maxima - monthly minima in each gridded climatology,
% if there are less than 7 months of data then it is NaN. From Feb 06 2024
% onward, this has been updated so that bias corrections < 1 are converted
% to a multiplicative integer, otherwise the amplitude was only higher for
% regions with a bias correction > 1
% GclimSTD = standard deviation of the monthly climatology at each grid,
% if there are less than 7 months of data then it is NaN
% Gseas = seasonal climatolgy - at the moment just calcaulted as the
% median of the monthly climatology for each season
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

function [Gclim, Gstd, Gstat, Ginter_var, GclimAMP, GclimSTD, Gseas] = make_global_gridded_monthly_clim(m_buff,res)
% Save a QF field to write to later to flag data outside 2 sigma of the
% climatology
% m.chla_raw_odmn_QF = 2.*ones(size(m.lat));
% m.chla_adj_odmn_QF = 2.*ones(size(m.lat));
% m.chla_sat_QF = 2.*ones(size(m.lat));
% m.fluor2chl_QF = 2.*ones(size(m.lat));

% Prep input data from m structure
% Get DOY vector
[yy, mm, dd, doy] = get_doy(m_buff.dt);
% convert lon
lontmp = m_buff.lon;
tt = find((lontmp)>180==1);
lontmp(tt) = lontmp(tt) - 360; clear tt

% Pre-make G arrays
gv = [{'chla_raw_odmn'},{'chla_npq_odmn'},{'chla_npq_zeumn'},...
    {'chla_npq_mldmn'},{'chla_drk_off_odmn'},{'chla_drk_off_zeumn'},...
    {'chla_drk_off_mldmn'},{'chla_sat'},{'sal_odmn'},...
    {'no3_odmn'},{'no3_zeumn'},{'no3_mldmn'},{'doxy_odmn'},{'doxy_zeumn'},...
    {'doxy_mldmn'},{'pH_odmn'},{'pH_zeumn'},{'pH_mldmn'},{'bbp_odmn'},...
    {'bbp_zeumn'},{'bbp_mldmn'},{'temp_odmn'},{'temp_zeumn'},{'temp_mldmn'},...
    {'poc_odmn'},{'poc_zeumn'},{'poc_mldmn'},{'fluor2chl'},{'rel_diff'},...
    {'diff'}];
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
% Note - the lat matrix will index from south to north, if you change this
% you'll need to update any functions that extract the ratio from the grids
% so that they are gathering data >= the upper bounds and < the lower
% bounds appropriately
for vv = 1:length(gvv)
    eval( ['[' (char(gvv(vv))) '.lon,'  (char(gvv(vv))) '.lat]' '=meshgrid(-180+res/2:res:180-res/2,-90+res/2:res:90-res/2)']);
end

% Loop through months
for h = 1:12
    for i = 1:size(Gclim.lat,1)
        for j = 1:size(Gclim.lon,2)
            clear idx
            % Was rounding lat & lon m_buff values, just removed July 11 2023
            idx = find(mm == h & m_buff.lat <= Gclim.lat(i,1) + res/2 & ...
                m_buff.lat > Gclim.lat(i,1) - res/2 & ...
                lontmp <= Gclim.lon(1,j) + res/2 & lontmp > Gclim.lon(1,j) - res/2);
            %             if length(idx)> 2
            %                 [i j] % 160 185
            %             else
            %             end
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
                        % Flag data in m that is > 2 sigma
                        %                     clear badidx; badidx = find(m.(char(gv{v}))(idx)) > 2*Gstd.(char(gv{v}))(i,j,h) + median(tmp2,2,'omitnan') || ...
                        %                         m.(char(gv{v}))(idx) < median(tmp2,2,'omitnan') - 2*Gstd.(char(gv{v}))(i,j,h);
                        %                     m.([char(gv{v}) '_QF'])(idx(badidx)) = 4;
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
    clear cmat; cmat = Gclim.(char(gv{v}));
    if strcmp(char(gv{v}),'fluor2chl')
        % Find clim values < 1 and convert before estimating amplitude
        clear idxl; idx = find(cmat < 1);
        cmat(idx) = -1./cmat(idx);
    else
    end
    GclimAMP.(char(gv{v})) = abs(max(cmat,[],3) - min(cmat,[],3));
    GclimSTD.(char(gv{v}))  = std(cmat,0,3,'omitnan');
    % If less than 6 months of observations occurred, set the Amp to NaN
    clear tmp; tmp = sum(~isnan(cmat),3);
    GclimAMP.(char(gv{v}))(tmp < 7)  = NaN;
    GclimSTD.(char(gv{v}))(tmp < 7)  = NaN;
end


%% Seasonal climatology - at this point just taking median of monthly clim
[~,mm,~,~] =  get_doy([1:365]);
ss = [12,1,2;3,4,5;6,7,8;9,10,11]; % Months per season in N. Hemi
ss_south =  [6,7,8;9,10,11;12,1,2;3,4,5];% Months per season in S. Hemi

for i = 1:4
    % I think this part was from when I was doing it on DOY basis
    %     moidx = find(mm == ss(i,1) | (mm == ss(i,2)) | (mm == ss(i,3)));
    %     modix_south =  find(mm == ss_south(i,1) | (mm == ss_south(i,2)) | (mm == ss_south(i,3)));
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