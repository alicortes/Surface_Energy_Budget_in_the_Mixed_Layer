%pcolo by FRAM
function h=pcolo(x,y,d,swt)

if nargin>=3
    if length(size(x))==1 
        warning('bad input. This part of the code is not finished');
        d(length(d(:,1))+1,:)=NaN;
        d(:,length(d(1,:))+1)=NaN;    
    
        xl=length(x);
        x(xl+1)=2*x(xl)-x(xl-1);    
        
        yl=length(y);
        y(yl+1)=2*y(yl)-y(yl-1);    
    else
        d(length(d(:,1))+1,:)=NaN;
        d(:,length(d(1,:))+1)=NaN;    
        
        xl=length(x(:,1)); 
        if xl>1; x(xl+1,:)=2*x(xl,:)-x(xl-1,:); end   
        xl=length(x(1,:)); 
        if xl>1; x(:,xl+1)=2*x(:,xl)-x(:,xl-1); end
        if size(x,1)>1
            x=x';
        end
        yl=length(y(:,1));
        if yl>1; y(yl+1,:)=2*y(yl,:)-y(yl-1,:); end   
        yl=length(y(1,:)); 
        if yl>1; y(:,yl+1)=2*y(:,yl)-y(:,yl-1); end   
        if size(y,1)>1
            y=y';
        end
        if size(d,1)~=length(y) %size(y,1) changed 2/11/06
            d=d';
        end
    end
    if nargin==4; %Adding a fourth argument flips the axes (used for square d)
        x=x'; y=y'; d=d'; 
    end
    if size(x,1) ~= size(d,1)
            x=x';
    end
    if size(y,1) ~= size(d,1)
            y=y';
    end
    h=pcolor(x,y,d);
else
    d=x;
    d(length(d(:,1))+1,:)=NaN;
    d(:,length(d(1,:))+1)=NaN;
    h=pcolor(d);
end
set(gcf,'color','w');
shading('flat');colorbar
%set(gca,'fontsize',14);
set(gcf,'renderer','zbuffer');
