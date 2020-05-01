function mldepth = ml_depth(time,T,depth,dT)
% ML_DEPTH:  This program calculates the mixed layer depth 
%           when high accuracy temperature data is available (RBR) using a
%           temperature difference (dT = 0.02oC) from the surface approach
%
% USAGE:
%       mldepth = ml_depth(T,depth,dT)
% 
% INPUT:
%       T            [time x depth]
%       depth        Horizontal
%       dT            Temperature difference from the surface. Default = 0.02oC  (MacIntyre et al 2002)
%
% OUTPUT:
%       Mixed Layer Depth -- used in many SEML calculations
%
% NOTES:
%       There are still some issues with how we are calculating ml_depth in
%       that sometimes we are missing the upper layer.
%
% VERSION
%      1.0 -- Pulled together from many locations to one function - csh
%      2.0 -- ACC Jan 2018 - Clean up the function to compute Mixed Layer Depth
%      using a temperature difference from the surface approach. 
%

clear mldepth;
depth = vert(depth)';

% Find the depth of the uppermost thermistor which records a
% temperature at least d_T ºC colder than the temperature recorded
% by the thermistor closest to the surface.

diff_T      = -T(:,2:end) + repmat(T(:,1),1,length(depth)-1);
diff_T      = abs(diff_T);
[mx ix]     = max(diff_T > dT,[],2);
ix = ix;% right at the surface
%ix = ix + 1;
ix(mx == 0) = length(depth);
mldepth     = depth(ix);
nx          = any(isnan(T'));
mldepth(nx) = nan;

% % really shallow AML
% diff_TS = diff_T(:,1); % temperature differences near the surface
% dum = find(diff_TS > 1); % when the diff_TS is above 1oC (arbitrary), make the AML very shallow
% mldepth(dum) = 0.05; % 5 cm, arbitrary
% AML must be >0
dumi = find(mldepth < 0.15);
mldepth(dumi) = 0.15; %force the MLD to be 20-30cm

% % Very large temperature differences near the surface
% dT > 0.2 : the depth of the shallow AML is defined to be diway between
% them.
diff_surf = T(:,1) - T(:,2);
dum = find(diff_surf > 0.2);
if  ~isempty(dum)
    mldepth(dum) = 0.5.*(depth(2)-depth(1)); %arbitrary
end


