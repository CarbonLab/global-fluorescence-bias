% V4 uses the extracted satellite GRID data previosuly matched in
% "get_satellite_grids_for_SPROF.m" and "extract_sat_data.m" to each
% profile and takes the median in the user defined buffer region.
% The buffer region is defined in km and converted to degrees using cosine
%
% Uses the satellite chl to estimate OD rather than the float chl
% Sets data to NaN if no satellite data available over same time as float
%
% OUTPUT:
% f.sat_chl: 8-day averaged satellite chl extracted at float position and near time
% f.chl_od: OD-averaged float fluoresence
% f.od: optical depth estimate
% f.zeu: euphotic depth estimate
%
% Written by Jacki Long (MBARI)

% Degrees to loop through as buffer distance for satellite data averaging,
% NOTE: The extracted satellite grid is 0.5 degrees. This means that at our
% largest search buffer (1/4 deg), any profiles above 60 degrees latitude
% will exceed our grid size (27.27km/(111.32*cos(60)) = 0.49 deg
% 1/4, 1/6, 1/8 degrees
% 0.25, 0.17, 0.125, 0.036 degrees
% 27.75, 18.5, 13.875, 4 km (converted at lat = 0, 111km/deg)
f.sat_buffer_km = [32, 24, 16, 12, 8, 4.63, NaN]; % Last option will find overlapping grid


%% Float data
[f.udate,ind] = unique(f.date);
f.udate = f.udate';
lat = f.lat(ind);
lon = f.lon(ind);
% Satellite longitude grid is from -180 to 180, float data are from 0 to
% 360 so convert first
lon(lon>180) = lon(lon>180)-360;

%% Pull Satellite Data Along Float Path
clear qlon qlat qdate
[prof] = unique(f.profile);
% Loop through profiles
for p = 1:length(prof)
    clear current_prof; current_prof = find(f.profile == prof(p));

    % Get delta in degrees
    clear dlat; dlat = abs(lat(p) - f.sat.lat(p,:));
    % Using absolute values so that near +/-180° profiles index correctly
%     clear dlon; dlon = abs(abs(lon(p)) - abs(f.sat.lon(p,:)));

    % If the profile is within 3deg of +/-180 border, switch to 0-360deg
    if 180 - abs(lon(p)) <  3
        % Switch satellite lon vec to 0-360
        clear lontmp; lontmp = f.sat.lon(p,:); 
        lontmp(lontmp<0) = 360 + lontmp(lontmp<0);
        % Switch float lon too if it is <0
        if lon(p)<0
            flon = 360 + lon(p);
        else
            flon = lon(p);
        end
        clear dlon; dlon = abs(flon - lontmp);
    else
        clear dlon; dlon = abs(lon(p) - f.sat.lon(p,:));
    end

    for j = 1:length(f.sat_buffer_km)

        if j < length(f.sat_buffer_km)
            % Find buffer distance in degrees for lon and lat
            % 1lat deg = 110.574 km
            % 1lon deg = 111.32*cos(lat)
            lon_lim = f.sat_buffer_km(j)./(111.32.*cosd(lat(p)));
            lat_lim = f.sat_buffer_km(j)./110.574;
            % Find lat and lon indeces that are less than the buffer in degs
            clear qlat; qlat = find(dlat <= lat_lim);
            clear qlon; qlon = find(dlon <= lon_lim);
            
            % Find overlapping pixel for nearest option
        else
            % Find exact crossover
           clear qlat; [~, qlat] = min(dlat);
           clear qlon; [~, qlon] = min(dlon);
        end

        % Save sat chla, takes the median of all satellite data found in
        % resolution chosen above, ex: f.sat.chl{1}(1:4,2:3)
        f.chla_sat(j,p) = median(f.sat.chl{p}(qlat,qlon),'all','omitnan');
        f.chla_std(j,p) = std(f.sat.chl{p}(qlat,qlon),[],'all','omitnan');
        f.nsat_pixels(j,p) = length(qlat) * length(qlon);
        f.nsat_pixels_valid(j,p) = sum(~isnan(f.sat.chl{p}(qlat,qlon)),'all','omitnan');

        % Where were data extracted
%         figure(1);clf
%         set(gcf,'Units','Inches');
%         m_proj('Robinson','long',[-179 179],'lat',[-90 90]);hold on
%         m_grid(['box'],'off','tickdir','in','linestyle','none','yticklabels',[],'xticklabels',[]);
%         % m_gshhs_l('patch',rgb('grey'),'edgecolor','k');hold on% - this function no
%         m_coast('patch','w');
%         m_plot(f.sat.lon(qlon),f.sat.lat(qlat),'r.','MarkerSize',12);


        % Get OD-averaged data
        % Kd is estimated using the satellite value for Chl
        clear calChlMLD; calChlMLD = f.chla_sat(j,p);
        clear Kd490_M07; Kd490_M07 = 0.0166 + 0.0773.*(calChlMLD).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
        Kd490_M07(Kd490_M07<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold (1986?), L&W Table 3.16
        clear KdPAR_M07; KdPAR_M07 = 0.0665 + 0.874.*(Kd490_M07) - 0.00121./(Kd490_M07);  %Morel et al. 2007 Eq. 9' (KdPAR over 2 optical depths)
        
        % Save the OD that was estimated based on each input chla type
        f.od_sat(j,p) = 1/Kd490_M07;
        f.zeu_sat(j,p) = 4.6/KdPAR_M07;

        for c = 1:length(chlopts)
            clear cvar; cvar = char(chlopts(c));
            clear gidx; gidx = find(~isnan(f.(cvar)(current_prof)) ...
                & ~isnan(f.press(current_prof)) & f.press(current_prof) >= 0);
            clear good_chl; good_chl = f.(cvar)(current_prof(gidx));
            clear good_press; good_press = f.press(current_prof(gidx));
            clear odidx; odidx = find(good_press <= f.od_sat(j,p));
            % Interpolate before averaging
            if length(odidx) > 1
                minpress = min(good_press(odidx));
                if minpress > 0
                    clear surf_chl; surf_chl = good_chl(odidx((good_press(odidx) == minpress)));
                    f.([(cvar) '_odmn'])(j,p) = median(interp1([good_press(odidx);0],[good_chl(odidx);surf_chl],[0:1:f.od_sat(j,p)]'),'omitnan');
                else
                    f.([(cvar) '_odmn'])(j,p) = median(interp1(good_press(odidx),good_chl(odidx),[0:1:f.od_sat(j,p)]'),'omitnan');
                end
            elseif length(odidx) == 1
                f.([(cvar) '_odmn'])(j,p) = good_chl(odidx);
            else
                % If the OD is shallower than the surface most observation
                %  by the float, we set this to NaN
                f.([(cvar) '_odmn'])(j,p) = NaN;
            end
        end
    end %end buffer option loop
end




