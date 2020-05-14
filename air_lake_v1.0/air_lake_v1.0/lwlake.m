function qlw = lwlake(Ts, Tz, cloud)

% Set default for cloud cover. 

if nargin == 2
    cloud = 0;
end 

if (cloud > 1 || cloud < 0)
    error('cloud parameter is out of range, choose a value between 0-1')
end

alb = 0.03;
Ce = 0.937e-5;
boltz = 5.67e-8;        % Stefan-Boltzman constant
em_water = 0.97;        % Constant used for freshwater as compared to 0.985
                        % for seawater.
                        
em_air = Ce .* (Tz + 273.15)^2;
qlwin = - em_air * boltz * ( 1+0.17*cloud^2)*(Tz + 273.15)^4;
qlwout = em_water * boltz * (Ts + 273.15)^4;
qlwref = - alb * qlwin;
qlwbal = qlwin + qlwout + qlwref;

qlw = [qlwin, qlwout, qlwref, qlwbal];

