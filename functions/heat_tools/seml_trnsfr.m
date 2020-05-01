function [CD,CE] = seml_trnsfr(z,L,CDN,CEN)
% CALCULATE TRANSFER COEFFICIENTS
% Compute transfer coefficients based on the integrated similarity 
% functions of Paulson (1970)

k=0.4;  PI=3.14159;

% compute integrated similarity functions
if (z/L) > 15.0; L = z/15.0; end
if (z/L) < -15.0; L = -z/15.0; end
if (z/L) < 0
    x = (1-16*(z/L))^0.25;
    sm  = 2.0*log((1.0+x)/2.0)+log((1+x^2.0)/2.0)-2.0*atan(x);
    sm  = sm + pi/2.0;
    se  = 2.0*log((1+x^2.0)/2.0);
elseif (z/L) < 0.5
    sm = -5*(z/L);
    se = sm;
elseif (z/L) < 10.0
    sm = 0.5/(z/L)^2.0-4.25/(z/L)-7*log(z/L)-0.852;
    se = sm;
else
    sm = log(z/L)-0.76*(z/L)-12.093;
    se = sm;
end

% compute transfer coefficients
CD = CDN/(1+CDN/k^2.0*(sm*sm-k*sm/CDN^0.5-k*sm*CDN^0.5/CDN));
CE = CEN/(1+CEN/k^2.0*(sm*se-k*se/CDN^0.5-k*sm*CDN^0.5/CEN));
