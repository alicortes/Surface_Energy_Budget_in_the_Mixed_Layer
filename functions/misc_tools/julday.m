function jd=julday(stime)
% JULDAY.M
% jd=julday(stime)
% Converts Matlab serial time to julday
% assuming jan 1 is day zero ...

% 12nov03 Brian Emery

jd=stime(:)-datenum(str2num(datestr(stime(:),10)),1,1);