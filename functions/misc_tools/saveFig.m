function saveFigPNG(figSaveDir,figNameEnd)
% SAVEFIG.M
% saveFig(figSaveDir,figNameStart)
% Loops over open figures and attempts to save them
% as publication quality in fig, png format.
%
% Set the figure name using set(gcf,'Name','blah') to include this
% information before the file extension. 

% Copyright (C) 2006-2010 Brian M. Emery
% 3feb06 

% FIG SAVING GRL (AGU) standards
% - 1200 dpi for black and white, 300-600 dpi for others
% - save as eps (vector graphics), view w/ adobe to check
% - for false-color imagery, shadding or texture use tiff
% - color figs in RGB and CMYK
% - All digital art submitted must be bitmap (Monochrome), grayscale, or
%   CMYK.
% - add labels to multi-part panels (eg a), b) etc)
%

% GET NUMBER OF OPEN FIGS
figNum=get(0,'Children');
%figNum=figNum(end:-1:1);

for i=1:length(figNum)
    
    figure(i)
    
    % try to get part of the file name from the figure name
    figNamePart = get(gcf,'Name'); 
    
    % if no name add a number
    if isempty(figNamePart)
        figNamePart = num2str(i,'%02.0f');
    end
    
    % SAVE THE FIGURE
    
    saveas(gcf,fullfile(figSaveDir,[figNameEnd figNamePart '.png']))   
    saveas(gcf,fullfile(figSaveDir,[figNameEnd figNamePart '.fig']))
    % saveas(gcf,fullfile(figSaveDir,[figNameEnd figNamePart '.jpg']))
    % print(gcf,'-deps','-r1200','-cmyk',fullfile(figSaveDir,'eps',[figNameEnd figNamePart '.eps']))
    
    disp(['SAVED: ' fullfile(figSaveDir,[figNameEnd figNamePart '.png'])])
    
    %close(figNum(i))

    %     else
    %         disp(['Fig ' num2str(figNum(i)) ' Has no user data and was not saved'])
    %         close(figNum(i))
    %     end

end

% save file creation info
%createdBy(figSaveDir)



end
