function eqax(FIG);

% This program standardizes the sizes and locations of every axis in the
% ensemble of Figure handles that are passed to it.
% 
% eqax(FIG) -- FIG: array of figure handles

axhdls = flipud(findobj(FIG,'type','axes','-and','tag',''));

axpos = cell2mat(get(axhdls,'position'));

xleft = max(axpos(:,1));
xrght = min(axpos(:,3));

axpos(:,1) = xleft;
axpos(:,3) = xrght;

for i = 1:length(axhdls)
    set(axhdls(i),'position',axpos(i,:))
end