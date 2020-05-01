function x = vert(x)
%  Turns any matrix or array into a vertical matrix or array.  No change if
%  the matrix or array is already vertical.

siz = size(x);

if siz(1) < siz(2)
    x = x';
else
    x = x;
end
