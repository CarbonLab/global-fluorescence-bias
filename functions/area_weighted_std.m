function weighted_std = area_weighted_std(var,ilat,ilon)
%%
% weighted_std = area_weighted_std(var,ilat,ilon)
%
% Weighted global standard deviation using area as the weighting function
%
% Input: 
%   var: variable of interest (2D or 3D w/ lon as x and lat as y dimension)
%   lat: centered grid points (i.e., -89.5:1:89.5);
%   lon: centered grid points (i.e.,  0.5:1:359.5);
%
% Output: 
%    weighted standard deviation over lat, lon bounds
%
% Functions Called: 
%   getArea
%   area_weighted_mean
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

w_std = NaN.*ones(o,1);
for i = 1:o
    clear xvar area xnan ynan ind
    xvar = reshape(var(:,:,i),[],1);
    area = reshape(aarea,[],1);

    % Find good values
    ind = find(isnan(xvar)==0 & isnan(area)==0);
    N   = length(ind);
    
    weighted_mean = area_weighted_mean(var(:,:,i),ilat,ilon);

    numerator   =  N .* nansum(area(ind) .* (xvar(ind) - weighted_mean).^2);
    denominator = (N - 1) .* nansum(area(ind));
    w_std(i,1)  = (numerator ./ denominator).^.5;
end

weighted_std = nanmean(w_std);

end
