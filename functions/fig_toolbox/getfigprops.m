function FPS = getfigprops(varargin)
% GETFIGPROPS   Gets the property values associated with the object handles
%     contained in FHS.handles, or the handles of the appropriate childern
%     of an array of figure handles and returns them to the FIGURE
%     PROPERTIES STRUCTURE.  
%
% VERSION:  
%      1.0 -- C. Helmle, 17 MAR 2005 (St. Patty's Day!!!  BEER ME!)
%                                     I blame any coding errors on the
%                                     beer.
% INPUT:
%
%     FIG:  Array of figure handles
%
%     -OR-
%
%     FHS:  Structure array containing the following fields AT THE
%           MINNIMUM (see subfunction GETFIGHDLS):
%
%     FHS.handles = 
%          Fig: <Figure Handle>
%         Axes: <Axes Handles>
%       XLabel: <XLabel Handles>
%       YLabel: <YLabel Handles>
%       ZLabel: <ZLabel Handles>
%        Title: <Title Handles>
%         Line: <Line Handles>
%       Legend: <Legend Handles>
%     Colorbar: <Colorbar Handles>
%
% OUTPUT
%      FPS = 
%      handles: <Structure of Object Handles>  (Potential input; see above)
%          Fig: <Structure of Figure Properties>
%         Axes: <Structure of Axes Properties>
%       XLabel: <Structure of XLabel Properties>
%       YLabel: <Structure of YLabel Properties>
%       ZLabel: <Structure of ZLabel Properties>
%        Title: <Structure of Title Properties>
%         Line: <Structure of Line Properties>
%       Legend: <Structure of Legend Properties>
%     Colorbar: <Structure of Colorbar Properties> 
%
% EXAMPLES:
%        FPS(N).Axes(M).XLim :  Returns the values for the XLim property 
%                               of the Mth Axes (not including Legends
%                               and Colorbars) from the top of the Nth figure.
%    FPS(N).Legend(M).String :  Returns the values for the String property
%                               of the Mth Legend from the top of the Nth
%                               figure
%      FPS(N).handles.Axes(M):  Returns the *handle* of the Mth Axes on the
%                               Nth figure.  This is best used with the SET
%                               command to change object properties.  
%
% ASSOCIATED FUNCTIONS:
%     SAVEFIGPROPS: Saves the properties to FPS variable in 'base'
%                   workspace.  
%

if ~isstruct(varargin{1})
    FIG = varargin{1};
    FPS = getfighdls(FIG);
else
    FPS = varargin{1};
end

for i = 1:length(FPS)
    %%% Assign object properties to FPS(i).<object(s)> structure %%%
    FPS(i).Fig    = get(FPS(i).handles.Fig);
    FPS(i).Axes   = get(FPS(i).handles.Axes);
    FPS(i).XLabel = get(FPS(i).handles.XLabel);
    FPS(i).YLabel = get(FPS(i).handles.YLabel);
    FPS(i).ZLabel = get(FPS(i).handles.ZLabel);
    FPS(i).Title  = get(FPS(i).handles.Title);
    FPS(i).Line   = get(FPS(i).handles.Line);
    FPS(i).Legend = get(FPS(i).handles.Legend);
    FPS(i).Colorbar = get(FPS(i).handles.Colorbar);

    %%% Remove unchangeable properties and properties that are handles %%%
    rmlabfld = {'BeingDeleted' 'Type'         'Extent'         'Parent'};
    rmlinfld = {'BeingDeleted' 'Type'         'Parent'};
    rmfigfld = {'BeingDeleted' 'CurrentPoint' 'CurrentAxes'    'CurrentObject'...
                'Type'           'Children'       'Parent'       ...
                'CurrentCharacter' 'Position'}; %AC hack - removed 'FixedColors'
    rmaxfld  = {'BeingDeleted' 'Children'     'CurrentPoint'   'Parent'       ...
                'Title'        'Type'         'XLabel'         'YLabel'       ...
                'ZLabel'       'CreateFcn'    'DeleteFcn'      'ButtonDownFcn'...
                'UIContextMenu' 'TightInset'};
    FPS(i).Fig    = rmfield(FPS(i).Fig,rmfigfld);
    FPS(i).Axes   = rmfield(FPS(i).Axes,rmaxfld);
    FPS(i).XLabel = rmfield(FPS(i).XLabel,rmlabfld);
    FPS(i).YLabel = rmfield(FPS(i).YLabel,rmlabfld);
    FPS(i).ZLabel = rmfield(FPS(i).ZLabel,rmlabfld);
    FPS(i).Title  = rmfield(FPS(i).Title,rmlabfld);
    FPS(i).Line   = rmfield(FPS(i).Line,rmlinfld);
    if ~isempty(FPS(i).Legend)   FPS(i).Legend   = rmfield(FPS(i).Legend,rmaxfld);   end
    if ~isempty(FPS(i).Colorbar) FPS(i).Colorbar = rmfield(FPS(i).Colorbar,rmaxfld); end
end


function FHS = getfighdls(fig)
% GETFIGHDLS   Gets handles to the following objects on a figure (or set of
%     figures) and returns them to the FIGURE HANDLES STRUCTURE:
%
% FHS.handles =
%          Fig: <Figure Handle>
%         Axes: <Axes Handles>
%       XLabel: <XLabel Handles>
%       YLabel: <YLabel Handles>
%       ZLabel: <ZLabel Handles>
%        Title: <Title Handles>
%         Line: <Line Handles>
%       Legend: <Legend Handles>
%     Colorbar: <Colorbar Handles>
%
% NOTES:
%     The structure is designed such that the handles are sequential from
%     top to bottom.
%
% INPUT:
%     fig:  array of figure handles

%%% Loop through each figure %%%
for i = 1:length(fig)
    %%% Assign Axes, Label, Legend, and Colorbar handles to FHS.handles %%%
    FHS(i).handles.Fig    = fig(i);
    FHS(i).handles.Axes   = findobj(FHS(i).handles.Fig,'type','axes','-and','tag','');
    %%% Sort axes handles vertically
    [xx ix] = sort(cell2mat(get(FHS(i).handles.Axes,'position')),'descend');
    FHS(i).handles.Axes   = FHS(i).handles.Axes(ix(:,2));
    
    %%% Lame Matlab quirk requires this TRY statement.  Output from GET is
    %%% a cell matrix if multiple handles are input.
    try
        FHS(i).handles.XLabel = cell2mat(get(FHS(i).handles.Axes,'XLabel'));
        FHS(i).handles.YLabel = cell2mat(get(FHS(i).handles.Axes,'YLabel'));
        FHS(i).handles.ZLabel = cell2mat(get(FHS(i).handles.Axes,'ZLabel'));
        FHS(i).handles.Title  = cell2mat(get(FHS(i).handles.Axes,'Title'));
    catch
        FHS(i).handles.XLabel = get(FHS(i).handles.Axes,'XLabel');
        FHS(i).handles.YLabel = get(FHS(i).handles.Axes,'YLabel');
        FHS(i).handles.ZLabel = get(FHS(i).handles.Axes,'ZLabel');
        FHS(i).handles.Title  = get(FHS(i).handles.Axes,'Title');
    end
    FHS(i).handles.Line   = flipud(findobj(FHS(i).handles.Axes,'type','line'));
    FHS(i).handles.Legend = flipud(findobj(FHS(i).handles.Fig,'tag','legend'));
    FHS(i).handles.Colorbar = flipud(findobj(FHS(i).handles.Fig,'tag','Colorbar'));
    % Not currently grabbing children of Colorbar and Legend Axes (i.e. XLabel,
    % YLabel,etc...).  If this is necessary, I suggest adding a variable
    % such as FHS.handles.CBXLabel and FHS.handles.LegXLabel, etc....
end