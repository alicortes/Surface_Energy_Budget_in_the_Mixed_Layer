function varargout = binaveall(varlen,r)
% FUNCTION BINAVEALL
%       This function will block average ALL variables in your current 
%       workspace of length = varlen. 
%
% USAGE:
%       bindata = binaveall(varlen,r)
%        varlen: length of the variables to search for
%             r: # data points to block average together (see binave)
%       bindata: Optional output is a structure
%
% EXAMPLES:
%       Let's say you are sitting in on a SEML output workspace, and you
%       want to binaverage all of the data to run some higher-order
%       operations.  You call BINAVEALL like this:
%
%       binaveall(length(doy),12)
%
%       This will average all of the variables in you workspace and will
%       overwrite them with the new block averaged data
%
%       bindata = binaveall(length(doy),12)
%
%       This will do the same thing, but instead of overwriting your data,
%       it will dump all of the variables into a structured variable, such
%       as the following:
%     
%       binave = 
%            doy: [34651 x 1] double
%           AirT: [34651 x 1] double
%             RH: [34651 x 1] double
%               :
%               :
%               etc...

%%% Determin all variables in 'base' workspace %%%
allvars = evalin('base','whos');

%%% Loop thru variables and binave all that are of length varlen %%%
n = 0 ;
for i = 1:length(allvars)
    if any(allvars(i).size == varlen)
        n = n+1;
        dim = find(allvars(i).size == varlen);
        varname(n) = {allvars(i).name};
        bindata.(allvars(i).name) = evalin('base',['binave(' allvars(i).name ',' num2str(r) ',' num2str(dim) ');']);
    end
end

%%% Overwrite if no output arguments %%%
if ~nargout
    for i = 1:length(varname)
        assignin('base',varname{i},bindata.(varname{i}));
    end
end


            
            
        