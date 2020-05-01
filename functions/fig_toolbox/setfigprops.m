function varargout = setfigprops(FHS,FPS)

%%% Assign object properties to FPS.<object(s)> structure %%%
% FHS: Figure Handle Structure
% FPS: Figure Property Structure

if length(FHS) ~= length(FPS)
    error(['Figure Handle Structure [FHS] does not match the size'...
        ' of Figure Propery Structure [FPS].']);
end

for i = 1:length(FPS)
    set(FHS(i).handles.Fig,FPS(i).Fig)
    for j = 1:length(FHS(i).handles.Axes)
        set(FHS(i).handles.Axes(j),FPS(i).Axes(j))
        set(FHS(i).handles.XLabel(j),FPS(i).XLabel(j))
        set(FHS(i).handles.YLabel(j),FPS(i).YLabel(j))
        set(FHS(i).handles.ZLabel(j),FPS(i).ZLabel(j))
        set(FHS(i).handles.Title(j),FPS(i).Title(j))
    end
    for j = 1:length(FHS(i).handles.Line)
        set(FHS(i).handles.Line(j),FPS(i).Line(j))
    end
    for j = 1:length(FHS(i).handles.Legend)
        set(FHS(i).handles.Legend(j),FPS(i).Legend(j))
    end
    for j = 1:length(FHS(i).handles.Colorbar)
    set(FHS(i).handles.Colorbar(j),FPS(i).Colorbar(j))
    end
    FPS(i).handles = FHS(i).handles;
end

varargout = {FPS};
