function matmax = max2(X)
% Gets maximum value in the entire matrix

[r c] = size(X);
Xnew = reshape(X,r*c,1);
matmax = max(Xnew);
end
