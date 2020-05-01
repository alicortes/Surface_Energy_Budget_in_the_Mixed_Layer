function matmin = min2(X)
% Gets minimum value in the entire matrix

[r c] = size(X);
Xnew = reshape(X,r*c,1);
matmin = min(Xnew);
end
