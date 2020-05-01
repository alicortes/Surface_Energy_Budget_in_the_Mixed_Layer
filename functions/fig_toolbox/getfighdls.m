function FHS = getfighdls(fig)
% GETFIGHDLS   Gets handles to the following objects on a figure (or set of
%     figures) and returns them to a structure:
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
%     top to bottom.  In other words, FHS(n).handles.Axes(1) will be the
%     uppermost axis (not including Legends and Colorbars) on figure "n".
%     Additionally, FHS(n).handles.Legend(2) will be the second uppermost
%     Legend on figure "n".  
%
%     This function is to be used in tandem with the following functions:
%     GETFIGPROPS:  Uses FHS as input to grab object properties
%     SAVEFIGPROPS: Saves the properties to FPS variable in 'base'
%                   workspace.  
%
% INPUT:
%     fig:  array of figure handles

%%% Loop through each figure %%%
for i = 1:length(fig)
    %%% Assign Axes, Label, Legend, and Colorbar handles to FHS.handles %%%
    FHS(i).handles.Fig    = fig(i);
    FHS(i).handles.Axes   = flipud(findobj(FHS(i).handles.Fig,'type','axes','-and','tag',''));
    %%% Sort axes handles vertically
    try
    [xx ix] = sort(cell2mat(get(FHS(i).handles.Axes,'position')),1,'descend');
    catch
    [xx ix] = sort((get(FHS(i).handles.Axes,'position')),1,'descend');
    end
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