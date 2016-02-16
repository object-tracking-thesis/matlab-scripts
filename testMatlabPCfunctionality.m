path1 = '/Users/markocotra/Google Drive/';
path2 = 'Thesis Work AF - Object Tracking in Point Cloud/Data/2011_09_26_drive_0017_sync/kitti/';

lidarData = loadLidarDir(strcat(path1,path2));

%% Load scene 
ptCloud = cell(size(lidarData));

for k = 1:length(ptCloud)
    ptCloud{k} = pointCloud(lidarData{k}(:,1:3));
end

%% Play scene
player = pcplayer([-40 40],[-40 40], [-3 40]);
while isOpen(player)
    for k = 1:length(ptCloud)
       view(player,ptCloud{k}) 
       pause(0.2)
    end

end

%% Test denoise

ptCloud1de = pcdenoise(ptCloud{1},'NumNeighbors',10,'Threshold',0.1);

figure
pcshow(ptCloud1de)
figure
pcshow(ptCloud{1})

%% Test downsample
close all
maxNumPoints = 20;

ptCloud1ds = pcdownsample(ptCloud{1},'nonuniformGridSample',maxNumPoints);

figure
pcshow(ptCloud1ds)
figure
pcshow(ptCloud{1})

%% Subtract static map

filteredLidar = cell(size(lidarData));

for k = 2:length(lidarData)
    filteredLidar{k} = subPC(lidarData{1},lidarData{k},0.1);
end

%% Port to ptCloud class 
filteredLidar{1} = lidarData{1};

ptCloudFiltered = cell(size(lidarData));
for k = 1:length(ptCloudFiltered)
    ptCloudFiltered{k} = pointCloud(filteredLidar{k}(:,1:3));
end

%% Display filtered frames 

player = pcplayer([-40 40],[-40 40], [-3 40]);
findobj('name', 'Point Cloud Player')

while isOpen(player)
    for k = 1:length(ptCloudFiltered)
       view(player,ptCloudFiltered{k})  
       pause(0.2)
    end

end

%% Try Denoise

ptCloudFDE = cell(size(ptCloudFiltered));

for k = 1:length(ptCloudFDE)
    ptCloudFDE{k} = pcdenoise(ptCloudFiltered{k},'NumNeighbors',50,'Threshold',0.005);
end

%% 

player = pcplayer([-40 40],[-40 40], [-3 40]);

while isOpen(player)
    for k = 1:length(ptCloudFDE)
       view(player,ptCloudFDE{k})  
       pause(0.2)
    end

end


















