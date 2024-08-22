%convert datenum to doy
%[yr mo day doy] = get_doy(date)
%Example:
%date = f.date

function [yr mo day doy] = get_doy(date)

DV = datevec(date);
mo = DV(:,2);
day = DV(:,3);
DV = DV(:,1:3);
DV2 = DV;
DV2(:, 2:3) = 0;
tmp = cat(2,DV(:,1), datenum(DV) - datenum(DV2));
yr = tmp(:,1);
doy = tmp(:,2);
return



