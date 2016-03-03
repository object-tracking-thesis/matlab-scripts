function [varargout] = plotMeasurement(X,Y,detected,Color)
% Plots a measurement (XY) differently depending on if detected or not. If
% detected, it will be illustrated with a 'x'. If not detected, by a '*' in
% the current timestep. History is illustrated with 'x' (for now, might 
% change this in the future). Color can be specified as generally accepted
% string, i.e. 'r' or as a 1x3 matrix. 
%
%     plotMeasurement(X,Y,detected,Color)
% p = plotMeasurement(X,Y,detected,Color)
%     
%     X        - x coordinates vector
%     Y        - y coordinates vector 
%     detected - designates if the current (i.e. last) element in X and Y 
%                is detected or not. 1 - detected, 0 - not detected. 
%     Color    - colorspec. 
%     p        - optional return parameter, figure handle. 

    if ~isempty(X) && ~isempty(Y)
        if detected 
            p1 = plot(X,Y,'x','Color',Color);
        else 
            plot(X(1:end-1),Y(1:end-1),'x','Color',Color);
            p1 = plot(X(end),Y(end),'*','Color',Color);
        end


        if length(nargout) > 1
            disp('Hola')
            varargout{1} = p1;
        end

    end

end