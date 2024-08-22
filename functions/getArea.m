function area = getArea(lat,lon)
%%
% area = getArea(lat,lon)
%
% Input:
%   lat: centered grid points (i.e., -89.5 : 1 : 89.5);
%   lon: centered grid points (i.e.,   0.5 : 1 :359.5);
%   
% Output: 
%   Global area in m^2

% Get grid information
nlon = length(lon); % number of longitude points
dlat = abs(lat(2) - lat(1)); % latitude spacing
Re   = 6.371E6; % radius of Earth

% Pre-allocate array
area = ones(length(lat), length(lon)) * NaN; 

% Loop over each latitude band to get area
for i = 1 : length(lat)
    latArea = 2*pi*Re^2 * abs(sind(lat(i) - dlat/2.) - sind(lat(i) + dlat/2.));
    area(i, :) = latArea / nlon;% m^2
end

end