function updateplotyy(ax)
%UPDATEPLOTYY updates the axes ticks for plots made with PLOTYY
% UPDATEPLOTYY(AX) updates the ticks and limits for the two axes 
% created by PLOTYY whose handles are in AX.  It is based on the code in
% PLOTYY that determines the original ticks and limits.
%
% Example:
%
%     [ax,h1,h2] = plotyy(1:10,sin(1:10),5:15,cos(5:15));
%     hold on
%     plot(-10:-1,cos(-10:-1))
%     updateplotyy(ax)


error(nargchk(1,1,nargin));

if numel(ax) ~= 2 || ~all(ishandle(ax)) || ~all(strcmp('axes',get(ax,'type')))
    error('Input must contain two axes handles')
end

% set the axes XLimMode to auto, so they both determine the best limits
% automatically
set(ax,'XLimMode','auto')

% get the new XLim values
xlim1 = get(ax(1),'XLim');
xlim2 = get(ax(2),'XLim');

% set the X limits to the widest range necessary
set(ax,'XLim',[min(xlim1(1),xlim2(1)) max(xlim1(2),xlim2(2))]);

% set the YLimMode to auto, so the Y limits are recalculated
set(ax,'YLimMode','auto');

% get the new YLim values
ylim1 = get(ax(1),'YLim');
ylim2 = get(ax(2),'YLim');

% determine whether either of the axes use log scale
islog1 = strcmp(get(ax(1),'yscale'),'log');
islog2 = strcmp(get(ax(2),'yscale'),'log');

if islog1, ylim1 = log10(ylim1); end
if islog2, ylim2 = log10(ylim2); end

% Find bestscale that produces the same number of y-ticks for both
% the left and the right.
[low,high,ticks] = bestscale(ylim1(1),ylim1(2),ylim2(1),ylim2(2),islog1 | islog2);

if ~isempty(low)
    if islog1,
      yticks1 = logsp(low(1),high(1),ticks(1));
      low(1) = 10.^low(1);
      high(1) = 10.^high(1);
      ylim1 = 10.^ylim1;
    else
      yticks1 = linspace(low(1),high(1),ticks(1));
    end

    if islog2,
      yticks2 = logsp(low(2),high(2),ticks(2));
      low(2) = 10.^low(2);
      high(2) = 10.^high(2);
      ylim2 = 10.^ylim2;
    else
      yticks2 = linspace(low(2),high(2),ticks(2));
    end

    % Set ticks on both plots the same
    set(ax(1),'ylim',[low(1) high(1)],'ytick',yticks1);
    set(ax(2),'ylim',[low(2) high(2)],'ytick',yticks2);
    set(ax,'xlim',[min(xlim1(1),xlim2(1)) max(xlim1(2),xlim2(2))])

    % Set tick labels if axis ticks aren't at decade boundaries
    % when in log mode
    if islog1
      decade1 =  abs(floor(log10(yticks1)) - log10(yticks1));
    end
    if islog2
      decade2 =  abs(floor(log10(yticks2)) - log10(yticks2));
    end
    if islog1 & any(decade1 > 0.1)
      for i=length(yticks1):-1:1
        ytickstr1{i} = sprintf('%3g',yticks1(i));
      end
      set(ax(1),'yticklabel',ytickstr1)
    end

    if islog2 & any(decade2 > 0.1)
      for i=length(yticks2):-1:1
        ytickstr2{i} = sprintf('%.3g',yticks2(i));
      end
      set(ax(2),'yticklabel',ytickstr2)
    end
    
else
    % Use the default automatic scales and turn off the box so we
    % don't get double tick marks on each side.  We'll still get
    % the grid from the left axes though (if it is on).
    set(ax,'box','off')
end


function [low,high,ticks] = bestscale(umin,umax,vmin,vmax,islog)
%BESTSCALE Returns parameters for "best" yy scale.

penalty = 0.02;

% Determine the good scales
[ulow,uhigh,uticks] = goodscales(umin,umax);
[vlow,vhigh,vticks] = goodscales(vmin,vmax);

% Find good scales where the number of ticks match
[u,v] = meshgrid(uticks,vticks);
[j,i] = find(u==v);

if islog % Filter out the cases where power of ten's don't match
  for k=length(i):-1:1
    utest = logsp(ulow(i(k)),uhigh(i(k)),uticks(i(k)));
    vtest = logsp(vlow(j(k)),vhigh(j(k)),vticks(j(k)));
    upot = abs(log10(utest)-round(log10(utest))) < 10*eps*log10(utest);
    vpot = abs(log10(vtest)-round(log10(vtest))) < 10*eps*log10(vtest);
    if ~isequal(upot,vpot),
       i(k) = [];
       j(k) = [];
    end
  end
end

if ~isempty(i)
  udelta = umax-umin;
  vdelta = vmax-vmin;
  ufit = ((uhigh(i)-ulow(i)) - udelta)./(uhigh(i)-ulow(i));
  vfit = ((vhigh(j)-vlow(j)) - vdelta)./(vhigh(j)-vlow(j));

  fit = ufit + vfit + penalty*(max(uticks(i)-6,1)).^2;

  % Choose base fit
  k = find(fit == min(fit)); k=k(1);
  low = [ulow(i(k)) vlow(j(k))];
  high = [uhigh(i(k)) vhigh(j(k))];
  ticks = [uticks(i(k)) vticks(j(k))];
else
  % Return empty to signal calling routine that we weren't able to
  % find matching scales.
  low = [];
  high = [];
  ticks = [];
end



%------------------------------------------------------------
function [low,high,ticks] = goodscales(xmin,xmax)
%GOODSCALES Returns parameters for "good" scales.
%
% [LOW,HIGH,TICKS] = GOODSCALES(XMIN,XMAX) returns lower and upper
% axis limits (LOW and HIGH) that span the interval (XMIN,XMAX) 
% with "nice" tick spacing.  The number of major axis ticks is 
% also returned in TICKS.

BestDelta = [ .1 .2 .5 1 2 5 10 20 50 ];
penalty = 0.02;

% Compute xmin, xmax if matrices passed.
if length(xmin) > 1, xmin = min(xmin(:)); end
if length(xmax) > 1, xmax = max(xmax(:)); end
if xmin==xmax, low=xmin; high=xmax+1; ticks = 1; return, end

% Compute fit error including penalty on too many ticks
Xdelta = xmax-xmin;
delta = 10.^(round(log10(Xdelta)-1))*BestDelta;
high = delta.*ceil(xmax./delta);
low = delta.*floor(xmin./delta);
ticks = round((high-low)./delta)+1;

%---------------------------------------------
function  y = logsp(low,high,n)
%LOGSP Generate nice ticks for log plots
%   LOGSP produces linear ramps between 10^k values.

y = linspace(low,high,n);

k = find(abs(y-round(y))<=10*eps*max(y));
dk = diff(k);
p = find(dk > 1);

y = 10.^y;

for i=1:length(p)
  r = linspace(0,1,dk(p(i))+1)*y(k(p(i)+1));
  y(k(p(i))+1:k(p(i)+1)-1) = r(2:end-1);
end

