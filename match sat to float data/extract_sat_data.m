% Program to match float data in latitude, longitude, and time to sat data
% A grid of satellite data surrounding the float position is extracted to
% be used later on so that we can change the buffer extraction region
% without re-processing. OD-estimates are not done in this function 
%
% Called in: "A_get_satellite_grids_for_SPROF.m"
%
% INPUT:
% Satellite data, temporal and spatial resolution are chosen in main
% program
% Float data
%
% OUTPUT:
% f.sat.chl: satellite data matrixof extracted pixels near float position
% and time
% f.sat.lat: latitude value (at center) per pixel 4km= [-89.979 to 89.979]
% f.sat.lon: logitude value (at center) per pixel 4km= [-179.979 to 179.979]
% f.sat.dt: datenum of satellite data (mean of start and end time listed
% in .nc file)
%
% Written by Jacki Long (MBARI)

% %% User choices - better to pull this chunk out of this function and into
% the main code so it only runs once
if ~isvarname('sat')
    %     % What satellite data product to use?
    %     vnm = 'chl_4km';
    %     tres = 'daily';
    %     % Resolution buffer for match-ups, in degrees
    %     res = 0.125;
    %% Program start
    % Grab satellite data for variable "vnm"
    if filesep == '/'
        satpath = ['/Volumes/MBGC_SatData/Data/Satellite/Downloaded/MODIS/MODIS_' tres];
    else
        satpath = ['X:/Data/Satellite/Downloaded/MODIS/MODIS_' tres];
    end
    sat = matfile([satpath '/' vnm '.mat' ],'Writable',false);
    % sat.date is now the mean of the start and end dates in the satellite
    % nc file
    sat_time = sat.date;
    sat_lat = sat.lat;
    sat_lon = sat.lon;
    %     chlopts = {'chla_raw','chla_drk_off','chla_npq'};
else
end

%% Float data
%- had to switch to using profile index below because there
% are sometimes more profiles than unique dates (e.g., two profiles ran
% under same date)
% [f.udate,ind] = unique(f.date);
% lat = f.lat(ind); 
% lon = f.lon(ind);
% % Satellite longitude grid is from -180 to 180, float data are from 0 to
% % 360 so convert first
% lon(lon>180) = lon(lon>180)-360;

%% Pull Satellite Data Along Float Path
clear qlon qlat qdate
[prof,ind] = unique(f.profile);
lat = f.lat(ind);
lon = f.lon(ind);
f.udate = f.date(ind);
% Satellite longitude grid is from -180 to 180, float data are from 0 to
% 360 so convert first
lon(lon>180) = lon(lon>180)-360;
% Flag for near 180deg profiles
f.near_180 = zeros(length(prof),1);


% Loop through profiles
for p = 1:length(prof)
    % If there's no processed satellite data, set NaN
    clear current_prof; current_prof = find(f.profile == prof(p));
    % Find closest lat/lon and time to float within defined res
    clear Del; Del = abs(lon(p) - sat_lon);
    [qlon] = find(Del < res);
    clear Del; Del = abs(lat(p) - sat_lat);
    [qlat] = find(Del < res);

    [f.fl_sat_time_diff(p),qdate] = min(abs(f.udate(p) - sat_time));

    % Because we are in degrees space, the length of qlat should = qlon
    if length(qlon) ~= length(qlat) % Then float is near 180deg boundary
        % Number of pixels missing
        d = length(qlat) - length(qlon);
        % Flag this profile
        f.near_180(p) = 1;
        % If on negative side, you want the last indeces added (positive side)
        if lon(p) < 0
            qlon2 = [length(sat_lon)-(d-1):length(sat_lon)]';
            % If on positive side, you want the first indeces added (negative side)
        elseif lon(p) >= 0
            qlon2 = [1:d]';
        end
    else
        qlon2 = [];
    end

    if f.udate(p) < max(sat_time)
        % If float is near 180deg boundary append satellite pixels
        if exist('d','var')
            f.sat.chl(p,1) = {[sat.var(qlat,qlon,qdate),sat.var(qlat,qlon2,qdate)]};
        else
            f.sat.chl(p,1) = {sat.var(qlat,qlon,qdate)};
        end
    else % Satellite record needs to be updated, set to NaN
        % Assumes there has been a previous profile to set size, but if
        % first profile is before sat record, then entire float is, so it's
        % okay if this breaks and we move to next float
        f.sat.chl(p,1) = {NaN(size(f.sat.chl{p-1,1}))};
    end

    f.sat.lat(p,:) = sat_lat(qlat);
    f.sat.lon(p,:) = sat_lon([qlon;qlon2]);
    f.sat.dt(p,1) = sat_time(qdate);

    clear d qdate qlon qlat qlon2
end




