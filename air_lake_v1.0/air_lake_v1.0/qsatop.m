  function q=qsatop(Ta, op, Pa)

% QSATOP: computes specific humidity at saturation. 
% q=QSAT(Ta) computes the specific humidity (kg/kg) at satuation at
% air temperature Ta (deg C) using Tetens' formula for saturation vapor
% pressure from Buck (1981), J. App. Meteor., 1527-1532. An option has been 
% added to calculate the specific humidity (kg/kg) from an updated 
% reference by Buck (1996). The dependence % on pressure is small (<0.5%) 
% and has been removed using a mean pressure of 1020 mb.  
%
%    INPUT:   Ta - air temperature  (deg C)
%             op - Options: 1 - Buck (1981); 2 - Buck (1996)
%             Pa - (optional) Pressure (mbars)
%
%    OUTPUT:  q  - saturation specific humidity  (kg/kg)
%
% Original code by Rich Pawlowicz (3/8/1997)
% Version 1.0 by Alex Forrest (3/23/2007)


if nargin==2,
  as_consts;
  Pa=P_default; % pressure in mb
end;

if op == 1
    a=(1.004.*6.1121*0.6220)./Pa;           % Changed 6.112 to 6.1121
    q=a.*exp((17.502.*Ta)./(240.97+Ta));
    else if op == 2
        a=(1.004.*6.1121*0.6220)./Pa; 
        q=a.*exp(((18.678 - Ta/234.5).*Ta)./(257.14+Ta));
        else
            disp('This is not a valid option, SH not calculated') 
        end
    end
end
  
        


