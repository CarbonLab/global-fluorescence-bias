% Function to filter float data for make_global_gridded_data
%
% Option to loop through each float and sets values outside X sigma to NaN
% Option to set nighttime casts of chla_raw to NAN
%
% Sets float data outside the satellite record to NaN
% Sets negative OD-averaged data to NaN
%
% Set small chl values to NaN
% Small for satellite = < 0.01 based on the minimum discrete data shown in
% Hu et al 2012 for the CI product which is used for the empirical fit
% Small value for float is < 0.015 which is Seabird's listed
% sensitivity for it's mid observation range (0 - 30 ug/L)
% file:///Users/jlong/Desktop/ECOPuckRevC=en.pdf
%
% Save the ratio of fluor2chl as the float chla_npq (based on the chla_raw
% data with get_NPQcorrand dark off applied) and MODIS 8 day 4km satellite data
% Saves a relative difference ratio as well
% Optional: night-only casts from chla_raw to chla_sat also saved
%
function [mqc] = filter_float_data(m,sat_dt_end)

ufl = unique(m.float_ID);

% Horizon depths averaged over 
depths = [{'zeumn'},{'mldmn'},{'odmn'}];
avg_vars = {'chla_raw','chla_drk_off','chla_npq','chla_adj'};
gv = [];
for i = 1:length(avg_vars)
    for j = 1:length(depths)
        gv = [gv, {[char(avg_vars(i)), '_', char(depths(j))]}];
    end
end
gv = [gv, {'chla_raw_surf'},{'chla_drk_off_surf'},...
    {'chla_npq_surf'}];

std_filt = 7;

mqc = m;

%% Chlorophyll
% Set negative profiles to NaN
% We do not do this before the averaging because the negative data are
% "real", but if they result in a negative value for the ratio, then we
% remove these so they do not go into the ratio estimation where satellite
% is never negative and a negative ratio would be given
for v = 1:length(gv)
    mqc.(char(gv(v)))(m.(char(gv(v))) < 0) = NaN;
    clear idx; idx = find(m.(char(gv(v)))  < 0.014); % 2,950 profiles with valid sat lost in chla_npq_odmn
    mqc.(char(gv(v)))(idx) = NaN;
end
% Set low satellite chl to NaN
clear idx; idx = find(m.chla_sat < 0.01); % 356 profiles 
mqc.chla_sat(idx) = NaN;
% clear idx; idx = find(mqc.chla_npq_odmn < 0.007); % 1,324 profiles with valid sat lost
% clear idx; idx = find(mqc.chla_npq_odmn < 0.015); % 4,672 profiles with valid sat lost
% clear idx; idx = find(mqc.chla_npq_odmn < 0.03); % 19,655 profiles with valid sat lost

%% BBP
other_vars = {'bbp'};
gv2 = [];
for i = 1:length(other_vars)
    for j = 1:length(depths)
        gv2 = [gv2, {[char(other_vars(i)), '_', char(depths(j))]}];
    end
end
% Set negative data to NaN
for v = 1:length(gv2)
    mqc.(char(gv2(v)))(m.(char(gv2(v))) < 0) = NaN;
    % clear idx; idx = find(m.(char(gv2(v)))  < 0.014); % What's the lower
    % limit for bbp?
    % mqc.(char(gv2(v)))(idx) = NaN;
end

% Set daytime casts in chla_raw to NaN
% night_only = 0;
% if night_only == 1 
%     m.chla_npq_odmn(m.night == 0) = NaN;
% else
% end

% Set any data after satellite record to NaN - should be obsolete
% ignore satellite var
for v = 1:length(gv)-1
    mqc.(char(gv(v)))(m.dt > sat_dt_end) = NaN;
end

% 
% % Set data outside 'std_filt' stdev to NaN
% for i = 1:length(ufl)
%     idx = find(mqc.float_ID == ufl(i));
%     for v = 1:length(gv)
%         cvar = char(gv{v});
%         mn = mean(mqc.(cvar)(idx),'omitnan');
%         st = std(mqc.(cvar)(idx),'omitnan');
%         clear badidx; badidx = find(mqc.(cvar)(idx) > std_filt.*st+mn | ...
%             mqc.(cvar)(idx) < mn-std_filt.*st);
%         mqc.(cvar)(idx(badidx)) = NaN;
%     end
% end

%% Save other variables
% Choose estimate for ratio and filter outliers as in function above
% ('std_filt' sigma is the limit)
cvar = 'fluor2chl';
mqc.(cvar) = mqc.chla_npq_odmn./mqc.chla_sat;
mqc.rel_diff = ((mqc.chla_npq_odmn - mqc.chla_sat)./mqc.chla_sat).*100;
mqc.diff = mqc.chla_npq_odmn - mqc.chla_sat;

% POC (mg m−3) = 31,200 * bbp_700nm (1/m) + 3.04 (Johnson et al., 2017).
mqc.poc_odmn = 31200.*mqc.bbp_odmn + 3.04;
mqc.poc_zeumn = 31200.*mqc.bbp_zeumn + 3.04;
mqc.poc_mldmn = 31200.*mqc.bbp_mldmn + 3.04;

% mqc.cphyto = 

% for i = 1:length(ufl)
%     idx = find(mqc.float_ID == ufl(i));
%     mn = mean(mqc.(cvar)(idx),'omitnan');
%     st = std(mqc.(cvar)(idx),'omitnan');
%     clear badidx; badidx = find(mqc.(cvar)(idx) > std_filt.*st+mn | ...
%         mqc.(cvar)(idx) < mn-std_filt.*st);
%     mqc.(cvar)(idx(badidx)) = NaN;
% end


% Save the same ratio with drk offset only & night casts
% cvar = 'fluor2chl_drknight';
% mqc.(cvar) = mqc.chla_drk_off./mqc.chla_sat;
% mqc.(cvar)(mqc.night == 0) = NaN;
% for i = 1:length(ufl)
%     idx = find(mqc.float_ID == ufl(i));
%     mn = mean(mqc.(cvar)(idx),'omitnan');
%     st = std(mqc.(cvar)(idx),'omitnan');
%     clear badidx; badidx = find(mqc.(cvar)(idx) > 3.*st+mn | ...
%         mqc.(cvar)(idx) < mn-3.*st);
%     mqc.(cvar)(idx(badidx)) = NaN;
% end



