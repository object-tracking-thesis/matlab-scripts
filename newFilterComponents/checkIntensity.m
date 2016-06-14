




% Load Raw Lidar Data 

pth = '/Users/markocotra//Google Drive/Thesis Work AF - Object Tracking in Point Cloud/Data/kitti_campus_frames/';


data = loadLidarDir(pth);


%%
N = 1;
plot(data{N}(:,1), data{N}(:,2),'k.'); axis equal

datData = data{N};

lbx = 12;
ubx = 17;
lby = 1;
uby = 3;

datData = datData(datData(:,1)<ubx,:);
datData = datData(lbx < datData(:,1),:);

datData = datData(datData(:,2)<uby,:);
datData = datData(lby < datData(:,2),:);

datData = datData(:,1:4);

for j = 1:length(datData)
    
    op = 1-double(datData(j,4));
    plot3(datData(j,1), datData(j,2), datData(j,3),'x','Color',op*[1 1 1]); axis equal; hold on;
end