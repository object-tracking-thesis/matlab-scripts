%% Map subtract ver 2 - knntree algorithm
clc; close all;

limit = 0.1;

path1 = 'matlabPCapi/';
path2 = 'crossing/';

staticMapFull = pcd2matlab(strcat(path1,path2,'0000000000.pcd'));
staticMap = staticMapFull(:,1:3);
lidarPlot(staticMap)
        axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])


%% Load lidardata files from path 

lidarData = loadLidarDir(strcat(path1,path2));
lidarDataRaw = lidarData;

% Subtract static point cloud
Num = length(lidarData);
for j = 1:Num
   lidarData{j} = subPC(staticMap,lidarData{j}(:,1:3),limit); 
   
   %pcdSpheric = xyzToSpheric(lidarData{j});
   %lidarData{j} = sliceSpheric(pcdSpheric,pi/720);
end




%%
% Plot Animation 
wd = 50; % Axis size
   figure('Name','Removing static map',...
          'Position', [50 100 1200 500]);
for j = 1:Num
   
   subplot(1,2,1)
        lidarPlot(lidarData{j})
        axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
        grid off; box off;
        str = sprintf('Frame: %d / %d', [j length(lidarData)]);
        title(str);
   subplot(1,2,2)
        lidarPlot(lidarDataRaw{j})
        axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
        grid off; box off;
        pause(0.5)
   
end
