function hh=vert_line(x_value,linestyle,colour)
%function hh=vert_line(x_value,[linestyle]);
if nargin<2; linestyle='--';end;
if nargin<3; colour='k';end;
x_value=x_value(:)';
hh=line([x_value;x_value],get(gca,'ylim')'*ones(size(x_value)),...
   'color',colour,'linestyle',linestyle);
