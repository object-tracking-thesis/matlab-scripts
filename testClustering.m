%% test points
points = [1 1; 0.9 0.9; 1.4 1; 3 3; 2.9 2.9; 1.1 1.1];

%% real points
%use ground-removed points here
lidarData = ptCloudFiltered{2};
size(lidarData.Location)

lidarData = pcdownsample(lidarData,'gridAverage', 0.5);

size(lidarData.Location)

points = lidarData.Location;%{8};

%% find clusters
%clustering(pcdXYZ,max-distance,min-points-per-cluster)
tic 
clusters = clustering(points,1,50);
toc

%% Port clusters to pointcloud class with separate color
clusteredPC = [];
for k = 1:length(clusters)
   [r, ~] = size(clusters{k});
   clusteredPC = [clusteredPC; [clusters{k} repmat((rand(1,3).*255),r,1)]]; 
end
%%
myJam = pointCloud(clusteredPC(:,1:3),'Color',uint8(clusteredPC(:,4:6)));

player = pcplayer([-40 40],[-40 40], [-3 40]);
view(player,myJam)  

%% plot all clusters
figure;
wd = 50;

for i = 1:length(clusters)
    h = lidarPlot(clusters{i});
    axis(wd*[-1 1 -1 1 -1/(0.5*wd) 1])
    grid off; box off;  
    color = rand(1,3);
    set(h,'MarkerEdgeColor', color,'MarkerFaceColor', color);
    hold on
end