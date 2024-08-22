%% NASA Satellite Data Processing
%
% This version has been updated for the recent filenaming convention change
% from NASA OC
%
% Input:
% .nc files saved from A_get_sat_data.m
%
% Output:
% This code processes NASA .nc files created in "A_get_sat_nc_files.m' and puts
% them into one large matrix for each variable. If the matrix doesn't exist,
% it is created. Once it exists, new satellite data files are appended to the
% existing matrix. Metadata are also saved. The .nc files are processed in
% batches of 10 at a time.
%
% Functions Called:
%   sat_nc2mat
%
% Author: Andrea Fassbender & Jacki Long
% Created: Dec 10 2022


clear;close all;clc

tres = 'monthly';

%% Data Access Path
if filesep == '/'
    dpath =  ['/Volumes/MBGC_SatData/Data/Satellite/Downloaded/MODIS/MODIS_' tres '/'];
else
    dpath = ['X:\Data\Satellite\Downloaded\MODIS\MODIS_' tres '\'];
end

%% .nc Input Fields
upath_nc      = {'bb4704km_nc'};%,'chl4km_nc','gamma_nc','nFLH_nc','pic_nc','poc_nc','zeu_lee_nc'}';
web_access_nc = {'https://oceancolor.gsfc.nasa.gov/l3/order/'};
%                  'https://oceancolor.gsfc.nasa.gov/atbd/pic/',...
%                  'https://oceancolor.gsfc.nasa.gov/atbd/poc/',...
%                  'https://oceancolor.gsfc.nasa.gov/atbd/'}';
mat_name_nc   = {'bb470_4km'};%,'chl4km','pic','poc','zeu_lee','gamma','nFLH'}';

%% Process New Satellite .nc Files & Append to Matricies

for ii = 1:length(upath_nc)
    start_process = now;

    clear f_name fname flist fdate sat xx fl u vv t
    flist = dir(fullfile([dpath filesep char(upath_nc(ii))], '*.nc'));
    fdate = datenum({flist.date});

    f_name = strings(length(flist),1);
    t = NaN(length(flist),1);
    for i = 1:length(flist)
        clear fn vinds yr d
        f_name(i) = [dpath filesep char(upath_nc(ii)) filesep flist(i).name];

        % Extract date from file names - not used anymore
%         fn = char(flist(i).name);
%         if contains(fn(1:2),'AQ') % New format uses full sat name, not abbv
%             yr = fn(12:15);
%             mo = fn(16:17);
%             d  = fn(18:19);
%             t(i,1) = datenum(yr,'yyyy') + str2double(d)-1;
%         else % Old format (e.g., 'A' for Aqua, 'T' for Terra)
%             yr = fn(2:5);
%             d  = fn(6:8);
%             t(i,1) = datenum(yr,'yyyy') + str2double(d)-1;
%         end
    end

    %  Find new data files
    u = char(upath_nc(ii));
    vv = u(1,1:end-3);
    mat_path = dpath;
    fl = dir(fullfile(mat_path, '*.mat'));
    if isempty(fl) == 1 || sum(strcmp({fl.name},[vv '.mat'])) == 0
        clear xx;xx = 1:length(flist);
    elseif find(strcmp({fl.name},[vv '.mat']) == 1)
        T = load([mat_path filesep vv '.mat'],'date');
%         [~,d] = setdiff(t,T.date);

        % I belive this was a minor date issue with the poc/pic files
        % So work around (5/18/2020)
        %         if ii < 3 && d == 5
        %             clear xx;xx = [];
        %         else
        %             clear xx;xx = d;
    end
end
disp([num2str(length(xx)) ' new files for ' char(mat_name_nc(ii))])

if isempty(xx) == 1
    disp(['No new files for ' char(mat_name_nc(ii))])
else
    % Code runs way faster in blocks when making initial matricies
    % (i.e., when processing the first large set of sat files)
    if length(xx) > 10
        clear xs;xs = 1:10:length(xx);
        for jj = 1:length(xs)
            if jj == length(xs)
                xind = xs(jj):length(xx);
            else
                xind = xs(jj):xs(jj+1)-1;
            end
            disp(['Process ' num2str(length(xind)) ' new files for ' char(mat_name_nc(ii))])
            mn = char(mat_name_nc(ii));
            fname = f_name(xx(xind));
            [sat] = sat_nc2mat(fname,web_access_nc,fdate(xx(xind)),mat_path,mn);
            disp(['Block ' num2str(jj) ' of ' num2str(length(xs)) ' complete.'])
        end
    else
        disp([num2str(length(xx)) ' new files for ' char(mat_name_nc(ii))])
        mn = char(mat_name_nc(ii));
        fname = f_name(xx);
        [sat] = sat_nc2mat(fname,web_access_nc,fdate(xx),mat_path,mn);
    end
end
process_time = now - start_process;
disp(['----- Variable ' char(mat_name_nc(ii)) ' done! Total time: ' num2str(24*60*process_time) ' -----'])

if exist('sat','var') == 1
    sat
end
disp(['----- Done processing NASA .nc Files! Total time: ' num2str(24*60*process_time) ' -----'])

