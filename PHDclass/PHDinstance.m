% Class representing a PHD filter (instance). Handles top level logic for
% the PHD, such as predicting & updating the intensity. 




classdef PHDinstance
    properties
        birthRFS
    end
    
    methods(Access = public)
        function predict(this)
            % TODO
        end
        
        function update(this, Z)
            % Z is set of measurements 
            % TODO
        end
    end
    
    methods(Access = private)
        function setBirthRFS(this)
           % TODO 
        end
    end
end