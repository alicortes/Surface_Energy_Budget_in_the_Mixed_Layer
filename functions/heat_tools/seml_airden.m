function airden = seml_airden(P,T,sph)
% CALCULATE AIR DENSITY
% This calculation is based on equation of state for air:
%     airden = P/(R*Tv)
% where R is the universal gas constant = 287.04 J/kg/øK
%     and Tv is the virtual temperature given as
%     Tv = Tx(1+0.6078xsph)
% INPUTS and UNITS:
%     sph     no units
%     T       C  NOTE: must convert to units of øK
%     P       mb  NOTE: must convert to units of Pa(N/m2) by multiplying by 100
%                       (i.e. 1 mb = 100 Pa)

a1=0.6078;	R=287.04;	ToK=273.15;
Tv = (T+ToK).*(1+0.6078*sph);
airden = (P*100)./(R*Tv);