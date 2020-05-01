function x = interp_nan(x,varargin)

% function x = interp_nan(x,t,frac,extrap,interp)
%
% Interpolates of NaN values in a vector time series including
% extrapolating to endpoints in preparation for spectral
% processing.  Works on columns.
%
% Input:
%	x - input time series
%	t - time
%   frac - (optional) maximum fraction of NaN data/column
%   extrap - specifies method or value for extrapolation of dat
%       extrap = 'extrap'=> linear extrapolation
%       extrap = nan or 0 (places nans or 0's in data out of range)        
%       default for extrap is not to extrapolate.
% Example:
%   x = interp_nan(x,t,[],NaN,'spline');
%       This will interpolate all *internal* NaNs with spline, but will not
%       extrapolate (because extrap = NaN)
%   x = interp_nan(x,t,[],'extrap','pchip');
%       This will interpolate all NaNs with pchip, and will extrapolate.
%
% Output:
%   x - interpolated output vector

% Updated: csh Feb 05 -- changed to linear interpolation from spline
% Update: BL, 04/05 -- Added options to extrap.
% Update: csh Feb/05 -- Added options to choose type of interpolation

t = (1:length(x))';
frac = 0.0;
extrap = [];
interp = 'linear';
if ~isempty(varargin)
    if ~isempty(varargin{1})
        t = varargin{1}; t = t(:);
    end
end
if length(varargin) > 1
    if ~isempty(varargin{2})
        frac = varargin{2};
    end
end
if length(varargin) > 2
    if ~isempty(varargin{3})
        extrap = varargin{3};
    end
end

if length(varargin) > 3
    if ~isempty(varargin{4})
        interp = varargin{4};
    end
end

s = size(x);

for n = 1:s(2)

    inan = find(isnan(x(:,n)));
    igood = find(~isnan(x(:,n)));
    if ~isempty(inan)
        if length(igood)/s(1) > frac
            if isempty(extrap)
                x(inan,n) = interp1(t(igood),x(igood,n),t(inan),interp);
            else
                x(inan,n) = interp1(t(igood),x(igood,n),t(inan),interp,extrap);
            end
        else
            x(:,n) = deal(NaN);
        end
    end

end

return