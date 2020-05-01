function matmed = median2(X)
% Gets median value for the entire matrix

[r c] = size(X);
Xnew = reshape(X,r*c,1);
matmed = median(Xnew);
end
