function A =lakeflux(Tz,Ts,Uz,z,rh,Pa);
% LAKEFLUX: computes sensible and latent heat fluxes and other variables.
% [LE,H,Cd,Ct,Cq]=LAKEFLUX(Tz,Ts,Uz,z,rh) computes the following variables:
%
%           LE = latent heat flux into ocean (W/m^2)
%           H = sensible heat flux into ocean (W/m^2)
%           Cd = drag coefficient
%           Ct = temperature transfer coefficient (Stanton number)
%           Cq = moisture transfer coefficient (Dalton number)
%           
% Based on the following metereological station data:
%
%           Tz = air temperature (deg C) measured at height zt (m)
%           Ts = air temperature (deg C) of the water surface
%           Uz = mean wind speed (m/s) measured at height zr (m) 
%           z = height (m) of the weather station
%           rh = relative humidity (%) measured at height zq (m)
%           Pa = air pressure (mb) (optional)
%                (note: If not input, Pa will default to 1020 mbar)
%           
% This code was inspired by the air-sea flux calculations (hfbulktc.m) 
% file found in air_sea_v1.1 written by Rich Pawlowicz (ref:
% http://www.eos.ubc.ca/~rich/). The variation used here is that the
% iterations are based on optimizing the integrated universal function by
% altering Cd, Ct and Cq after the method detailed in Launiainen, 1995 and
% Heikenheimo et al., 1999. The intent of this code is to apply these
% calculations to lakes in the absence usable data for the eddy correlation
% method. It should be noted that this will only work for in boundary layer
% of near-neutral stability (e.g. zeta < 0.5). The eventual intention is to
% upgrade this for other regions. 
% 
% References:
%
% Launianen, J., (1995). "Derivation of the relationship between the 
% Obukhov Stability Parameter and the Bulk Richardson Number for Flux-
% Profile Studies." Boundary-Layer Meteorology, vol. 76, 165-179.
%
% Heikinheimo, M., Kangas, M., Tourula, T., Venalainen, A. & Tattari, S.,
% (1999). "Momentum and heat fluxes over lakes Tamnaren and Raksjo
% determined by the bulk-aerodynamic and eddy-correlation methods. 

%% Calculation of sensible and latent heat flux

% Last updated on March 23, 2007

if nargin == 5,
  as_consts;
  Pa=P_default;     % pressure in mb
end;


%% Define constants

% Bring in pre-defined values

as_consts;

% Define constants

tol = 1e-15;        % Tolerance for convergence
g = 9.81;           % m/s^2
k = 0.42;           % von Karmen constants (Hogstrom, 1996; Heikinheimo et 
                    % al., 1999) 

K = 273.15;         % celcius 
eps_air = 0.62197;  % Used in calculation of specific humidity
cp = 1000;          % J/kg K 
b = 5.2;            % Webb, 1970			
zq	= 3.00e-6;      % m	Heikinheimo et al., 1999
zt	= 3.00e-6;      % m	Heikinheimo et al., 1999
alpha = 0.035;      % Charnock constant (Garratt, 1992; Heikinheimo et al.,
                    % 1999). This differs from open ocean with Charnock 
                    % values from 0.011 to 0.018 used for open oceans.
beta = 0.11;        % limiting Roughness Reynolds # for aerodynamically 
                    % smooth flow (Smith et al., 1996; Heikinheimo et al., 
                    % 1999) 
qs = qsatop(Tz, 2); % Calculates saturation specific humidity using qsatop
qz = (0.01.*rh).*qs;
nu = viscair(Tz);
                                   


%% Define reference values

Cinit = [1.40e-3, 1.80e-3, 1.80e-3];

                    % Where the first term is the drag coefficient, 
                    % Cd (Heikinheimo et al., 1999). The second is the 
                    % Dalton number (Imberger, 1983) and the 3rd is the
                    % Stanton number (Imberger, 1983)

%% Start calculating values

C = Cinit;

TzK = Tz + K;
ToK = TzK*(1+((1/eps_air)-1)*qz);
To = ToK - K;

rho_air=(100*Pa)./(gas_const_R*ToK);  

L = 2.5007996e6 - (0.0023536e6 * Ts);
LE = rho_air * C(2) * (qs - qz) * Uz * L;
H = rho_air * cp * C(3) * (Ts - Tz) * Uz;
ustar = (C(1) * Uz^2)^(1/2);
zeta = (-z * g * k * (1 + 0.61 * ToK * cp * (LE/L)/H))/((ustar^3)...
    * ToK * rho_air * cp);

% Assume a near-neutral stability (zeta < 0.5)						

if zeta >= 0.5
    error('In a region of not near-neutral stability (zeta > 0.5)')
end

univ = -b * zeta;   

% Intergrated universal function (Launiainen, 1995)			

zm = (alpha * ustar^2)/g + beta * nu / ustar;

% Using formula from Heikinheimo et al., 1999	

Cnew = [];

Cnew(1) = (k^2)*(log(z/zm)-univ)^(-2);
Cnew(2) = (k^2)*((log(z/zm)-univ)^(-1))*((log(z/zq)-univ)^(-1));
Cnew(3) = (k^2)*((log(z/zm)-univ)^(-1))*((log(z/zt)-univ)^(-1));

%% Use the objective function to minimize the error

ob = @(C)((C(1)-Cnew(1))^2 + (C(2)-Cnew(2))^2 + (C(3)-Cnew(3))^2);

C = fminsearch(ob, Cinit, optimset('TolX', tol));

% Seperates out each of the coeffients from the optimized vector.

Cd = C(1);
Ct = C(2);
Cq = C(3);

%% Write the output vector

A = [LE,H,Cd,Ct,Cq];

                    
                    


