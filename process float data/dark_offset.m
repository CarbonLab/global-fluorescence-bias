% CORRECT RAW CHL FOR IN SITU DARK IN CHL SPACE 
% Finds at least 5 profiles that go deeper than 900m
% Find min chl value for entire profile
% Takes the median of these 5 mins from each profile
%
% Code adapted from a method sent by Josh Plant
%
% If a negative min value is found in the profile, it is used

% Find data deeper than 900m
tP = find(f.press > 900 & ~isnan(f.chla_raw));
% Find number of profile with deep data
uC = unique(f.profile(tP));

% tP  = d.raw(:,iP) > 900 & ~isnan(d.raw(:,iCHL)); % data greater than 900m
% tmp = d.raw(tP,:); % C
% uC = unique(tmp(:,iSTA)); % Cycles
if isempty(uC) || size(uC, 1) < 5
    disp(['No deep profiles available for in situ CHL ',...
        'dark estimate']);
    CHL_Dark = NaN;
else
    chl_min = ones(5,1) * NaN;
    for i = 1:5 % step through cycles
        tC = find(f.profile == uC(i));
        chl_min(i) = min(f.chla_raw(tC),[],'omitnan'); % min over whole profile, min is usually not at depth due to non chl fluor
    end
    CHL_Dark = median(chl_min,'omitnan');
    clear chla_dark_corr; chla_dark_corr = f.chla_raw - CHL_Dark; % correct all chla
end
clear  chl_min tC tP uC