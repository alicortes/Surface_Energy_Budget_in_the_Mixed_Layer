function doy = dnum2doy_ow(dnum)
% Converts datenumber to doy (day of year)
%Edited by Adam Crowe, 3/25/2014 for over winter chains. The doy here
%exceeds 365 and continues from the year of the first data point.
[y,m,d,h,mi,s] = datevec(dnum);
ynum = datenum(y,0,0,0,0,0);
doy = dnum-ynum(1);
end