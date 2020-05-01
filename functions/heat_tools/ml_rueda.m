function mldepth = ml_rueda(doy,RHO,depth,drhodz)
% ML_DEPTH:  This program calculates the mixed layer depth 
%           when low accuracy temperature data is available (HOBO) using a
%           density gradient approach (drhodz = 0.003)
%
% USAGE:
%       mldepth = ml_rueda(doy,RHO,depth,drhodz)
% 
% INPUT:
%       RHO          [doy x depth]
%       depth        Horizontal
%       doy          vertical
%       drhodz         Density gradient. Default = 0.003kg/m3/m (Rueda et al. 2007)
%
% OUTPUT:
%       Mixed Layer Depth -- used in many SEML calculations
%
% NOTES:
%       This approach tends to select the top of the metalimnion with the 
%       exception of the strongest diurnal stratification
%
% VERSION
%      1.0 -- ACC Jan 2018 - Cleaned up and renamed variable in MLD_FR.m function
%

d_10c = [depth(1):.1:max(depth)];
mldepth = nan*ones(length(doy),1);   % doy for seml

for i = 1:length(doy) 
    
    B = unique(RHO(i,:));
    %if length(B) == length(RHO(i,:))
    RHO_10c = interp1(depth, RHO(i,:)', d_10c);
    rho = RHO_10c';
    delta_rho = gradient(rho);
    delta_rho = abs(delta_rho);
      if any(delta_rho > drhodz)
      ixdr = find(delta_rho> drhodz,1);
      mldepth(i) = d_10c(ixdr);
      else
       mldepth(i) = nan;
      end
    %else 
    %mldepth(i) = nan;
    %end
    
end
mldepth = interpNaN_works(mldepth);
dumi = find(mldepth < 0.15);
mldepth(dumi) = 0.5; %force the MLD to be 50 cm

clear i ixdr B T_10c rho delta_rho d_10c
