%% Batch save .nc NASA Ocean Color satellite data
%First, place data request at https://oceancolor.gsfc.nasa.gov/l3/order/
%select "standard" (or provisional, whatever), satellite, time and spatial resolution, and "mapped"
%data. There are also options to either extract data (regional extraction),
%or if you just download it, you will get global data
%After ordering data, a list of https' will appear, copy these to the
%"NASA_OC_files_to_update.txt" file and remove the *NRT.nc files
%
%The program below will batch download all the nc files
%
% After running this program: convert to mat files using "Process_NASA_nc_hdf_files.m"
% If you haven't downloaded data in more than a month, you'll need to make
% a new app key:
% Go here: https://oceandata.sci.gsfc.nasa.gov/appkey/
% Log in and request a new key and update the one below
% It may also be necessary to be logged into your NASA ocean data account 
% on your primary web browser for this to work
%
% Written by: Jacki Long (MBARI) December 10 2022

%% User inputs
appkey = 'a925ccc2a09cb32eedc16049ffa4203c2d89fac7';
% Choose variable output directory
var = 'bb4704km_nc'; %chl4km_nc
tres = 'monthly';
% Update your path if need be
if filesep == '/'
    bpath = '/Volumes/MBGC_SatData/Data/Satellite/Downloaded/';
else
    bpath = 'X:/Data/Satellite/Downloaded/';
end

%% Program begins
txt = importdata([bpath 'NASA_OC_files_to_update.txt']);
strtxt = char(txt);
filenames = strtxt(:,49:end);
% Options for websave
opts = weboptions;
opts.Timeout = 10; % Increase wait time
% Save an error log
err = struct('file',[],'url',[]);

%% Start looping through text lines
disp(['Starting file process ' char(datetime('now'))])
for i = 1:length(txt)
    try
        % Add appkey to URL name
        url = [strtxt(i,:) '?appkey=' appkey];
        filename = strcat([bpath 'MODIS' , filesep, 'MODIS_', tres, filesep, var,filesep, filenames(i,:)]);
        % Save .nc file from web URL to filename
        websave(filename, url, opts);
    catch
        % If above fails, log the error message
        err.file = [err.file; filenames(i,:)];
        err.url = [err.url ; strtxt(i,:)];
    end
    if i == length(txt)/4
        disp(['Fourth way done ' char(datetime('now')) ' '])
    elseif i == length(txt)/2
        disp(['Half way done ' char(datetime('now'))])
    end
end
disp(['Finished with satellite files ' char(datetime('now'))])

%% Save error log
save([bpath 'MODIS' filesep 'MODIS_', tres, filesep, var,filesep, 'errlog.mat'],'err')

% DING

