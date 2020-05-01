function spchum = seml_spchum(T,P,RH)
% CALCULATE SPECIFIC HUMIDITY
% This calculation of specific humidity is based on the following:
% #1  vp/P=sph/(0.622+(1-0.622)xsph)
% rearranging
% #2  sph = (0.622xvp)/(P-0.378xvp)
% #2 used to calcualte saturation specific humidity (sphsat) from
%     saturation vapor pressure (svp)
% specific humidity related to relative humidity (RH) as,
% #3  RH = r/rw
% where r is the mixing ratio equal to sph/(1-sph) and rw is the
%     saturation mixing ratio which yields
% #4  RH = sph x (1 - sphsat) / ((sphsat x (1 - sph))
% rearranging
% #5  spchum = RH x sphsat / (1 - sphsat + RH x sphsat)
% #5 used to calculated specific humidity (sph) from saturation
%     specific humidity (sphsat)
% INPUTS and UNITS:
%     svp     mb
%     P       mb
%     RH      %

%copied directly from seml by FRAM 1/6/06

a1=0.62197;

svp = satvap(T,P);
sphsat = a1*svp./(P-(1-a1)*svp);
spchum = RH./100.*sphsat./(1-sphsat+RH./100.*sphsat);
