%% lidarPlot
% Plots lidar data with scatter3. Input argument is nx3 matrix.
% Optionally returns plot handle
function [varargout] = lidarPlot(lidarData)
    q = 0.3; 
    h = scatter3(lidarData(:,1),lidarData(:,2),lidarData(:,3),... 
       '.',...
       'MarkerEdgeColor', q*ones(1,3),...
       'MarkerFaceColor', q*ones(1,3));
   str = sprintf('# of Points: %d',length(lidarData));
   text(1,1,40,str);
   
   if nargout > 0
       varargout{1} = h;
   end
       
end