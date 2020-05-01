function interval = commonint(time,varargin)
% COMMONINT -- This function determines the most common time step
%       (interval) in a time array.  Useful for problem data where the
%       datalogger failed or something went horribly wrong with the data
%       collection and you end up with irregular samples.
% USAGE
%         interval = commonint(time,tol)
%         tol: tolerance (in seconds)
%         time: time vector in days or seconds

tol = 1;
if nargin > 1; tol = varargin{1}; end

% If time span is < 365, we are assuming time is in units [days]
if range(time) < 365
    mult = 24*3600;
else
    mult = 1;
end

dtime = diff(time);                % dt 
udt   = unique(dtime);             % unique dts
if length(udt) == 1                % rls change- generates error if there aren't different time steps
interval = udt;
perbad = 0;
else     
n     = hist(dtime,udt);           % distribution of dts
[n I] = sort(n);                   % sort distribution
udt   = udt(I);                    % sort udts
dudt =  udt(end) - udt;            % Subtract most common dt to get deviation from it
ddtix  =  dudt > tol/mult;         % Identify which udts are outside of our tolerance
nbadtime = sum(n(ddtix));          % Count # of times tolerance is violated
perbad = nbadtime/length(time);    % Percentage of time out of tolerance
[mx ix] = max(n);                  % Index of most common udt
interval = udt(ix);                % Identify most common udt

end

if perbad > 0.40
    qda = questdlg([{'Warning:  The time interval for this data is highly irregular.'};...
        {' '}; {['The interval of ' num2str(interval*mult,'%6.2f') ' seconds occurs only '...
        num2str(mx/length(dtime)*100,'%2.0f') '% of the time']}; {' '}; {'Would you like to...'}],'Warning!',...
        'Continue','Keyboard','Continue');
    switch qda
        case 'Continue'
        case 'Keyboard'
            disp(' ')
            disp('Please determine what interval you would like to use and')
            disp('assign it to the variable ''interval''.')
            disp(['(Interval currently set at ' num2str(interval*mult,'%6.2f') ' seconds.)'])
            disp(' ')
            disp('Type ''return'' to continue this progarm')
            disp(' ')
            keyboard

    end

end
