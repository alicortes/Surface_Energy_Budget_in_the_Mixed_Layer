function savefigprops(varargin)
% SAVEFIGPROPS   Saves FIGURE PROPERTIES STRUCTURE to a *.mat file.
%       Default settings save FPS to the originating *.mat file so long as
%       that file has the variable Info.filename which contains the path
%       and name of the file.  
%
% VERSION:  
%      1.0 -- C. Helmle, 17 MAR 2005 (St. Patty's Day!!!  BEER ME!)
%                                     I blame any coding errors on the
%                                     beer.
% INPUT:
%      savefigprops('FPSname','Choose Mat File')
%      'FPSname' : The name of a new FPS variable so as not to overwrite
%                  the original.
%      'FigIx'   : Figure handle index (i.e. FPS.handles.Fig(FigIx))

%%% Get FPS from 'base' and update with current figure properties %%%
FPS = evalin('base','FPS;');

if nargin == 0
    FPSname = 'FPS';
    FigIx = 1:length(FPS);
elseif nargin > 0
    if isnumeric(varargin{1})
        FigIx = varargin{1};
        FPSname = 'FPS';
    else
    FPSname = varargin{1};
    FigIx = 1:length(FPS);
    end
end

for i = 1:length(FigIx)
FPS(FigIx(i)) = getfigprops(FPS(FigIx(i)));
end

%%% Place the updated FPS* in 'base' again %%%
assignin('base',FPSname,FPS);

%%% Try to save with Info.filename.  Otherwise, prompt for location 
try
evalin('base',['save(''-append'',Info.filename,''' FPSname ''')']);
catch
[fn pn] = uigetfile('*.mat','Save FPS to a *.mat File');
evalin('base',['save(''-append'',''' fullfile(pn,fn) ''',''' FPSname ''')']);
end





