function [sat] = sat_nc2mat(fname,web_access,fdate,mat_path,mat_name)
%
% [sat] = sat_nc2mat(fname,web_access,fdate,mat_path,mat_name)
% [sat] = sat_nc2mat(fname,web_access_nc(ii),fdate(xx(xind)),mat_path,mn);
% Generic code to convert global satellite data file(s) to a formatted
% matrix. Bad values are set to NaN. Matricies are oriented correctly.
%
% Input Variables:
%    fname = array of strings including full path and file names. only
%           include new files that will be appended to the existing matrix, 
%           though a check is included below.
%    mat_path = path for data matrix that conatins or will contain variable
%    web_access = link to data source
%    mat_name = name of matrix file
%    fdate = date file was downloaded from source
%
% Output Variables:
%   V: data structure that includes
%       .hdr = variable name
%       .units = units
%       .var = variable data (rotated 90 deg)
%       .lon = longitude (centered on pixel)
%       .lat = latitude (centered on pixel)
%       .date_end = datemumber (end date of 8 day composite)
%       .date_start = datemumber (start date of 8 day composite)
%       .date = datemumber (mean of date_start and date_end)
%       .source = source file name.
%       .web_access = website the data were downloaded from
%       .date_downloaded = date file was downloaded from web
%       .date_processed = date flie was processed into .mat file
%
% Functions Called:
%   latlon_find
%
% Authors: Andrea Fassbender & Jacki Long
% Created: 4/18/2020
% Last updated:
% Further updates made by Jacki Long (MBARI) from December 10 2022 on

%% Process NASA .nc Files to one .mat File for Each Variable

tstart = now;
for j = 1:length(fname)
    clear I V z var x tt 
    
    % Extract file source info. variable name, and fill value
    I = ncinfo(char(fname(j)));
    x = strcmp({I.Attributes.Name},'product_name')==1;
    V.source = I.Attributes(x).Value;
    V.hdr = I.Variables(1).Name;
    
    x = strcmp({I.Variables(1).Attributes.Name},'units')==1;
    V.units = I.Variables(1).Attributes(x).Value;
    x = strcmp({I.Variables(1).Attributes.Name},'_FillValue')==1;
    fill = I.Variables(1).Attributes(x).Value;
    
    % Load variable and change fill values to NaNs
    var = ncread(char(fname(j)), V.hdr);
    z = reshape(var,[],1) == fill;
    var(z) = NaN;
    z = reshape(var,[],1) < 0;
    var(z) = NaN;
    
    % Save beginning time
    x = strcmp({I.Attributes.Name},'time_coverage_start')==1;
    tt = char(I.Attributes(x).Value);
    V.date_start = datenum(tt(1:10),'yyyy-mm-dd');

    % Save ending time
    x = strcmp({I.Attributes.Name},'time_coverage_end')==1;
    tt = char(I.Attributes(x).Value);
    V.date_end = datenum(tt(1:10),'yyyy-mm-dd');
    
    % Set date to middle of composite
    V.date = mean([V.date_start,V.date_end]);
   
    % Determine lat and lon of variable
    [lat(:,1),lon(:,1),~,~] = latlon_find(var);
    
    V.lat = lat;
    V.lon = lon;
    V.var = rot90(var);
    
%     figure(1);clf
%     h = pcolor(rot90(var));
%     set(h,'edgecolor','none')
    
    % Save data to temporary matrix until processing the last batch file
    if exist('sat','var') == 0
        sat = V;
    else
        sat.var = cat(3, sat.var, V.var);
        sat.date = [sat.date; V.date];
        sat.source = vertcat(sat.source, V.source);
        sat.date_end =  [sat.date_end; V.date_end];
        sat.date_start =  [sat.date_start; V.date_start];
    end
    
    if j == length(fname)
        % List of .mat files in mat_path
        fl = dir(fullfile(mat_path, '*.mat'));
        if isempty(fl) == 1 || sum(strcmp({fl.name},[mat_name '.mat'])) == 0
            % This matrix doesn't exist yet. Create it!
            sat.date_processed = now;
            sat.date_downloaded = fdate(j);
            sat.web_access = web_access;
            save([mat_path '/' mat_name '.mat'], '-struct', 'sat', '-v7.3')
        elseif sum(strcmp({fl.name},[mat_name '.mat'])) == 1
            tmat = sat; clear sat
            % Matrix already exists. Find new data files
            T = load([mat_path '/' mat_name '.mat'],'date');
            
            [~,new_vals] = setdiff(tmat.date,T.date);
            if isempty(new_vals) == 0
                % Append to matrix without loading it to memory
                sat = matfile([mat_path '/' mat_name '.mat'],'Writable',true);
                
                [nrows, ~] = size(sat,'date');
                clear iy;iy = nrows+1:nrows+length(new_vals); 
                sat.date(iy,1) = tmat.date(new_vals);
                sat.source(iy,:) = tmat.source(new_vals,:);
                sat.date_processed(iy,:) = now; 
                sat.date_downloaded(iy,:) = fdate(j);
                
                [nrows, ncols, ndeep] = size(sat,'var');
                clear id;id = ndeep+1:ndeep+length(new_vals);
                sat.var(1:nrows,1:ncols,id) = tmat.var(:,:,new_vals); 
            else
                clc
                disp('-- no new data to add to the matrix --')
                sat = [];
            end
        end
        tot = now - tstart;
        disp(['Total run time: ' num2str(24*60*tot) ' min'])
    end
    
    if j == round(length(fname)/2)
        tot = now - tstart;
        disp(['Half way. Run time: ' num2str(24*60*tot) ' min'])
    end
end

end