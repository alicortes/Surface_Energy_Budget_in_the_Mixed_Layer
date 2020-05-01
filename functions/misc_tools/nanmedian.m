function y = nanmedian(x,dim)
%NANMEDIAN Median value, ignoring NaNs.
%   M = NANMEDIAN(X) returns the sample median of X, treating NaNs as
%   missing values.  For vector input, M is the median value of the non-NaN
%   elements in X.  For matrix input, M is a row vector containing the
%   median value of non-NaN elements in each column.  For N-D arrays,
%   NANMEDIAN operates along the first non-singleton dimension.
%
%   NANMEDIAN(X,DIM) takes the median along the dimension DIM of X.
%
% BE verision based on:
%   STATS_NONAN
%   [xbar,stdev,hi,lo,median]=stats_noNaN(x)
%   Computes the mean, std, and ranges of each row of a
%   matrix after removing NaN's. 
%  
%   See also mean_noNaN.m
%
% The output of stats_noNaN is transposed to mimic behavior of matlab's
% nanmedian.m

[xbar,stdev,hi,lo,medn]=stats_noNaN(x');
y=medn;


end
