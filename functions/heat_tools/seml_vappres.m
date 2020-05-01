function vapprs = seml_vappres(P,sph)
% COMPUTE VAPOR PRESSURE (from sph)

epslon = 0.62197;
vapprs = P*sph./(epslon+(1.0-epslon)*sph);