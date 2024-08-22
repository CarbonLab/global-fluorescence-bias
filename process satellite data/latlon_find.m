function [lat,lon,resx,resy] = latlon_find(var)
%
% [lat,lon,resx,resy] = latlon_find(var)
%
% Determine latitude, longitude, grid resolution of global field
%
% Input:
%   var = matrix of interest
%
% Output
%    lon = longitude (-180 to 180)
%    lat = latitude (-90 to 90)
%    resx = x grid size/resolution
%    resy = y grid size/resolution
%

[n,m]=size(var);
if n > m
    resx = 360./n;
    resy = 180./m;
elseif m > n
    resx = 360./m;
    resy = 180./n;
end
lon = -180+(resx/2) : resx : 180-(resx/2);
lat =  -90+(resy/2) : resy : 90-(resy/2);

end

