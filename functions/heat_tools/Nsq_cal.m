function [n2] = Nsq(rho,varargin);
%- function [n2] = Nsq(rho,z(opt));
%- Calculates buoyancy frequency from density
%- INPUTS:
%-  1) rho - (N x M) density in kg/m3
%-  2) z - optional (if z is a scalar, it indicates spacing between points.
%-                  if z is a vector it should be the coordinates of rho
%-                  in space (meaning depth).
%-
%- A better version of this code is Phil Morgan's sw_bfrq.m
%- \\scully\user5\matlibrary\users\wshaw\science\seawater\sw_bfrq.m

 mean_rho=nanmean(nanmean(rho));
 
 if nargin>1,
     z = varargin{1};
 else
     z = 1;
 end

if prod(size(rho))==length(rho), %- rho is a vector
     drhodz = gradient(rho,z);
else %- rho is a matrix

    [fx drhodz] = gradient(rho,1,z);
end

 n2 = 9.81*drhodz/mean_rho;
 %n=imag(sqrt(nsq_f))*3600/(2*pi);
return