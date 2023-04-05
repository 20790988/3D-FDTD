classdef Num2eng_Restorer
    % Num2eng_Restorer is for num2eng internal use only
    methods (Static, Hidden)
        function obj = loadobj(obj)
            figH = gcf;
            restoreH = addlistener(figH,'FileName','PostSet',@restore_num2eng);
            restoreH.Callback = @(varargin)restore_num2eng(figH, restoreH);
        end
    end
end

function restore_num2eng(figH, listenerH)
    setappdata(figH, 'restoreNum2eng', true);
    num2eng(figH);
    delete(listenerH);
end