function log = iswithin(x,y,tol)
% This function returns a logical value of 1 if x and y are within a
% tolerance range of tol of one another.  
% Example:
% 
%      log = iswithin(10,11,2)
%       log =
%            1

if abs(x-y) <= tol
    log = logical(1);
else
    log = logical(0);
end
