function atpres = seml_atpres(LkElev)
% CALCULATE ATMOSPHERIC PRESSURE
% This calculation is based on Berry et al (1945) Handbook of Meteorology
%     P = 1013.25x(1-0.0065xelev/288)^5.25587
% where P is in mb and elev is in m

a1=1013.25; a2=0.0065; a3=288; a4=5.25587;
atpres = a1*(1-a2*LkElev/a3)^a4;

