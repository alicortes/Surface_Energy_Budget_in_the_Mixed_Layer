function dnum = doy2dnum(doy,year)
% Converts doy (day of year) to datenumber
% USAGE:
% dnum = doy2dnum(doy,year)

y = datenum(num2str(year),'yyyy');
dnum = y + doy-1;



