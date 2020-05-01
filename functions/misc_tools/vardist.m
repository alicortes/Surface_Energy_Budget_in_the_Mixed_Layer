function varargout = vardist(X)
%%% This program distributes the columns of a matrix to individual
%%% variables.  

l = length(X(1,:));

for i = 1:nargout
    varargout(:,i) = {X(:,i)};
    if i == nargout & l > nargout
        varargout(:,i) = {X(:,i:end)};
        disp('Warning: Number of colums exceeds output variables')
    end
end
