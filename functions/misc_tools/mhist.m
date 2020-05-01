function N = mhist(A,X)
% Gets histogram of the entire matrix

siz = prod(size(A));
A = reshape(A,siz,1);
if nargin == 1
N = hist(A);
bar(N)
elseif nargin == 2
N = hist(A,X);
bar(X,N)
end
end

