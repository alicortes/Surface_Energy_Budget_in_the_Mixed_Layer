function y = interpNaN_works(x)

y=x;

bd=isnan(x);
gd=find(~bd);

bd([1:(min(gd)-1) (max(gd)+1):end])=0;

y(bd)=interp1(gd,x(gd),find(bd),'linear');
y(bd)=interp1(gd,x(gd),find(bd),'linear');