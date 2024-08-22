% V5 copied from V4, satellite array data are now converted to vectors for
% indexing - this removes complications in the previous extract_sat_code of
% needing the longitude to be in a certain order for near 180deg profiles
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
% f.*_(od/mld/zeu)mn: median of variable * through 1st optical depth (OD),
% euthpotic zone (zeu) or mixed layer depth (MLD)
% f.od_sat: optical depth estimate from satellite CHL
% f.zeu_sat: euphotic depth estimate from satellite CHL
%
% Written by Jacki Long (MBARI)

% Degrees to loop through as buffer distance for satellite data averaging
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
    %     clear dlat; dlat = abs(lat(p) - f.sat.lat(p,:));
    % Using absolute values so that near +/-180° profiles index correctly
    %     clear dlon; dlon = abs(abs(lon(p)) - abs(f.sat.lon(p,:)));

    % Switch to vector indexing instead
    clear m n;  [m,n] = size(f.sat.chl{p}); % m = lat, n = lon
    sat_chl_vec = reshape(f.sat.chl{p},1,m*n)';

    clear lat_vec; lat_vec = repmat(f.sat.lat(p,:)',1,length(f.sat.lon(p,:)));
    lat_vec = reshape(lat_vec,m*n,1);
    clear lon_vec; lon_vec = repmat(f.sat.lon(p,:),length(f.sat.lat(p,:)),1);
    lon_vec = reshape(lon_vec,m*n,1);

    % If the profile is within 3deg of +/-180 border, switch to 0-360deg
    if 180 - abs(lon(p)) <  3
        % Switch satellite lon vec to 0-360
        clear lontmp; lontmp = lon_vec;
        lontmp(lontmp<0) = 360 + lontmp(lontmp<0);
        % Switch float lon too if it is <0
        if lon(p)<0
            flon = 360 + lon(p);
        else
            flon = lon(p);
        end
        clear dlon; dlon = abs(flon - lontmp);
    else
        clear dlon; dlon = abs(lon(p) - lon_vec);
    end

    clear dlat; dlat = abs(lat(p) - lat_vec);

    % Save surface values of data, ~44.6 was max OD found in float data
    surfidx = find(f.press(current_prof)< 44 & ~isnan(f.chla_npq(current_prof)),1,'last');
    if ~isempty(surfidx)
        f.chla_raw_surf(1,p) = f.chla_raw(current_prof(surfidx));
        f.chla_drk_off_surf(1,p) = f.chla_drk_off(current_prof(surfidx));
        f.chla_npq_surf(1,p) = f.chla_npq(current_prof(surfidx));
        f.surf_depth(1,p) = f.press(current_prof(surfidx));
    else
        f.chla_raw_surf(1,p) = NaN;
        f.chla_drk_off_surf(1,p) = NaN;
        f.chla_npq_surf(1,p) = NaN;
        f.surf_depth(1,p) = NaN;
    end


    for j = 1:length(f.sat_buffer_km)

        if j < length(f.sat_buffer_km)
            % Find buffer distance in degrees for lon and lat
            % 1lat deg = 110.574 km
            % 1lon deg = 111.32*cos(lat)
            lon_lim = f.sat_buffer_km(j)./(111.32.*cosd(lat(p)));
            lat_lim = f.sat_buffer_km(j)./110.574;
            % Find lat and lon indeces that are less than the buffer in degs
            %             clear qlat; qlat = find(dlat <= lat_lim);
            %             clear qlon; qlon = find(dlon <= lon_lim);
            idx = find(dlat <= lat_lim & dlon <= lon_lim);

            % Find overlapping pixel for nearest option
        else
            % Find exact crossover
            %            clear qlat; [~, qlat] = min(dlat);
            %            clear qlon; [~, qlon] = min(dlon);
            % Find minimum in lat and lon
            tmp = abs(dlat + dlon);
            % Collect minimum index
            [~, idx] = min(tmp);

        end

        % Save sat chla, takes the median of all satellite data found in
        % resolution chosen above, ex: f.sat.chl{1}(1:4,2:3)
        f.chla_sat(j,p) = median(sat_chl_vec(idx),'all','omitnan');
        f.chla_std(j,p) = std(sat_chl_vec(idx),[],'all','omitnan');
        f.nsat_pixels(j,p) = length(sat_chl_vec(idx));
        f.nsat_pixels_valid(j,p) = sum(~isnan(sat_chl_vec(idx)),'all','omitnan');
        %
        %         % Where were data extracted
        %         figure(1);clf
        %         set(gcf,'Units','Inches');
        %         m_proj('Robinson','long',[-179 179],'lat',[-90 90]);hold on
        %         m_grid(['box'],'off','tickdir','in','linestyle','none','yticklabels',[],'xticklabels',[]);
        %         % m_gshhs_l('patch',rgb('grey'),'edgecolor','k');hold on% - this function no
        %         m_coast('patch','w');
        %         m_plot(lon_vec(idx),lat_vec(idx),'r.','MarkerSize',12);


        % Get depth horizon values
        % Kd is estimated using the satellite value for Chl
        clear calChlMLD; calChlMLD = f.chla_sat(j,p);
        clear Kd490_M07; Kd490_M07 = 0.0166 + 0.0773.*(calChlMLD).^0.6715;  %Morel et al. 2007 Eq. 8 empirical fit to NOMAD + LOV data
        Kd490_M07(Kd490_M07<0.0224) = 0.0224;  % Kd490 can't be less than pure-water minimum; Austin & Petzold (1986?), L&W Table 3.16
        clear KdPAR_M07; KdPAR_M07 = 0.0665 + 0.874.*(Kd490_M07) - 0.00121./(Kd490_M07);  %Morel et al. 2007 Eq. 9' (KdPAR over 2 optical depths)

        % Save the OD that was estimated based on each input chla type
        f.od_sat(j,p) = 1/Kd490_M07;
        f.zeu_sat(j,p) = 4.6/KdPAR_M07;

        % MLD estimate - copied from section in
        % get_NPQcorr_currentBGCversion.m -
        % FIND DENSITY MLD
        mld_den_threshold = 0.03; % Dong et al 2008
        den = density(f.sal(current_prof),f.temp(current_prof));
        goodidx = ~isnan(den); % Some profiles might have NaN salinity vals (e.g., WMO 1902303)
        ref_ind  = find(f.press(current_prof(goodidx))<=25, 2, 'last'); % find ref points for mld calc

        tmld = den - mean(den(ref_ind),'omitnan') <= mld_den_threshold; % mixed layer

        if all(tmld ==0) || all(tmld ==1) || sum(~isnan(den)) == 0
            f.mld(j,p) = NaN;
        else
            f.mld(j,p) = max(f.press(current_prof(tmld))); % depth of mixed layer
        end

        depths = [{'zeu_sat'},{'mld'},{'od_sat'}];
        dappend = [{'_zeumn'},{'_mldmn'},{'_odmn'}];
        for c = 1:length(odvars)
            clear cvar; cvar = char(odvars(c));
            if isfield(f,cvar)
                clear gidx; gidx = find(~isnan(f.(cvar)(current_prof)) ...
                    & ~isnan(f.press(current_prof)) & f.press(current_prof) >= 0);
                clear good_vals; good_vals = f.(cvar)(current_prof(gidx));
                clear good_press; good_press = f.press(current_prof(gidx));
                for dd = 1:length(depths)
                    outvar = [(cvar) char((dappend(dd)))];
                    cdepth = char(depths(dd));
                    didx = find(good_press <= f.(cdepth)(j,p));
                    % Interpolate before averaging
                    if length(didx) > 1
                        minpress = min(good_press(didx));
                        if minpress > 0 % Then add 0m to data for interpolation
                            clear surf_val; surf_val = good_vals(didx((good_press(didx) == minpress)));
                            f.(outvar)(j,p) = median(interp1([good_press(didx);0],[good_vals(didx);surf_val],[0:1:f.(cdepth)(j,p)]'),'omitnan');
                        else
                            f.(outvar)(j,p) = median(interp1(good_press(didx),good_vals(didx),[0:1:f.(cdepth)(j,p)]'),'omitnan');
                        end
                    elseif length(didx) == 1 % Then just equals this value
                        f.(outvar)(j,p) = good_vals(didx);
                    else
                        % If the OD is shallower than the surface most observation
                        %  by the float, we set this to NaN
                        f.(outvar)(j,p) = NaN;
                    end

                end % end dd loop for depths
            else % This variable isn't on the float, set to NaN
                 for dd = 1:length(depths)
                    outvar = [(cvar) char((dappend(dd)))];
                    cdepth = char(depths(dd));
                    f.(outvar)(j,p) = NaN;
                 end
            end
        end % end c loop for odvars
    end %end buffer option loop
end


