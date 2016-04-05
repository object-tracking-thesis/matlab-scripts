% 
clustersDup = load('test.mat');
clustersDup = clustersDup.clusters;
size(clustersDup)

tic
num = 1:length(clustersDup);
A = cellfun('length',clustersDup);
smaller = A < 5;

remove = smaller.*num';
remove = remove(remove > 0);

clustersDup(remove) = [];
size(clustersDup)
toc

%% 

% Alpha = pdist2(points,points);
% size(Alpha)
% Alpha = Alpha < 1;
% % remove mini clusters
% 
% 
% 
% tic
% Alpha = unique(Alpha, 'rows');
% toc
% size(Alpha)


%% 
lidarData = ptCloudFiltered{10};
size(lidarData.Location)

%lidarData = pcdownsample(lidarData,'gridAverage', 0.5);

size(lidarData.Location)

points = lidarData.Location;%{8};
scatter3(points(:,1),points(:,2),points(:,3),10,'.')

 wd = 50;
   axis(wd*[-1 1 -1 1 -1/(0.25*wd) 1])
%
%tic
Z = linkage(points,'single','euclidean');
figure(1)
dendrogram(Z)
%
c = cluster(Z,'Cutoff',0.5,'Criterion','distance');
figure(2)
histogram(c,20000)
%toc
%scatter3(points(:,1),points(:,2),points(:,3),10,c,'.')

%%

y = zeros(size(c));

for i = 1:length(c)
    y(i) = sum(c==c(i));
end

y = y>200;

col = y.*c;

%

lom = col == 0;
points(lom,:) = [];
clas = col(col > 0);
%%
scatter3(points(:,1),points(:,2),points(:,3),1,myColors,'.')

 wd = 50;
   axis(wd*[-1 1 -1 1 -1/(0.25*wd) 1])

%%

% Load lidarDir 
path1 = '/Users/markocotra/Google Drive/';
path2 = 'Thesis Work AF - Object Tracking in Point Cloud/Data/2011_09_26_drive_0017_sync/kitti/';

lidarData = loadLidarDir(strcat(path1,path2));

%% 
nomLidar = cell(lidarData);

for k = 1:length(lidarData)
    nomLidar{k} = pointCloud(lidarData{k}(:,1:3));
end


%% Subtract static map

filteredLidar = cell(size(lidarData));

for k = 2:length(lidarData)
    filteredLidar{k} = subPC(lidarData{1}(:,1:3),lidarData{k}(:,1:3),0.1);
end
filteredLidar{1} = lidarData{1};



%% Cluster Scene 
clusteredLidar = cell(size(filteredLidar));

for k = 2:length(filteredLidar)
    k
    [points, col] = clusterLidar(filteredLidar{k},0.5,200);
    myColors = class2color(col);
    myPC = pointCloud(points, 'Color', uint8(myColors));
    
    clusteredLidar{k} = myPC;
end

%%

player = pcplayer([-40 40],[-40 40], [-3 40]);
findobj('name', 'Point Cloud Player')

while isOpen(player)
    for k = 2:length(clusteredLidar)
       view(player,clusteredLidar{k})  
       pause(0.2)
    end

end
clusteredLidar{2}.Location
%% 

y = 2*mat2gray(filteredLidar{k}(:,1:3))-1;

tic 
[points, col] = clusterLidar(y,1,200);
toc

%%
player1 = pcplayer([-40 40],[-40 40], [-3 40]);
player2 = pcplayer([-40 40],[-40 40], [-3 40]);
player3 = pcplayer([-40 40],[-40 40], [-3 40]);



while isOpen(player)
    for k = 2:length(clusteredLidar)
       view(player1,clusteredLidar{k})  
       view(player2,filteredLidar{k})  
       view(player3,nomLidar{k})  
       pause(0.2)
    end

end














