function weighted_mean = area_weighted_mean(var,ilat,ilon)
%%
% weighted_mean = area_weighted_mean(var,ilat,ilon)
%
% Weighted global mean using area as the weighting function
%
% Input: 
%   var: variable of interest (2D or 3D w/ lon as x and lat as y dimension)
%   lat: centered grid points (i.e., -89.5:1:89.5);
%   lon: centered grid points (i.e.,  0.5:1:359.5);
%
% Output: 
%   weighted average over lat, lon bounds
%
% Functions Called: 
%   getArea
%

AREA = getArea(ilat,ilon);

%% Make sure Area and Var are the same size/orientation 

[n,m,o] = size(var);
[na,ma] = size(AREA);

if n ~= na & n == ma
    aarea = transpose(AREA);
elseif n == na & m == ma
    aarea = AREA;
else
    disp('Area sizing problem!')
end

w_mean = NaN.*ones(o,1);
for i = 1:o
    clear xvar area xnan ynan ind
    xvar = reshape(var(:,:,i),[],1);
    area = reshape(aarea,[],1);
    
    % Find good values
    ind = find(isnan(xvar)==0 & isnan(area)==0);
    
    w_mean(i,1) = nansum(xvar(ind) .* area(ind)) ./ nansum(area(ind));
end

weighted_mean = nanmean(w_mean);

end