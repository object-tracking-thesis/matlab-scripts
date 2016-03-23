%% clear all
clear all;
clc;

%% load geometrical data for static
pathWalls = 'static/walls.txt';
pathPoles = 'static/poles.txt';
pathRoadEdges = 'static/roadEdges.txt';
walls = walls2matlab(pathWalls);
poles = poles2matlab(pathPoles);
roadEdges = roadEdges2matlab(pathRoadEdges);

staticZeroGeo = [106473.25371600001; 6406253.6631990001; 133.126];

%% test plot all walls and poles
figure;
hold on
for i = 1:length(walls)
    plotCubes(walls{i}(1:3)',walls{i}(4),walls{i}(5),walls{i}(6),walls{i}(7:9),0,0)
end
for i = 1:length(poles)
    plotCylinder(poles{i}(2:4)',poles{i}(1),poles{i}(5),poles{i}(6))
end

%% load lidar and oxts data for live
path1 = '~/Downloads/thesis/share/pcap_scenarios/';
path2 = 'car/';
path3 = 'oxts/';
%path4 = 'pcd/2000to2200/';
path4 = 'pcd/1550to1650/';

oxts = loadOxtsDir(strcat(path1,path2,path3));
lidarData = loadLidarDir(strcat(path1,path2,path4));

Num = length(lidarData);

for i=1:Num
    lidarData{i}(:,4:7) = [];
end

%% plot lidarData
figure;
for j=1:Num
    test = lidarData{j};
    pcshow(test)
    pause(0.3)
end

%% test to cut away the forest part of the lidarData
for j=1:1
    lidarData{j}(lidarData{j}(:,2) < 0,:) = [];
end

%% rotate and translate the live frames according to their gps data
k = 1550; %frame offset to the gps data
liveFrames = cell(1,Num);
offset = cell(1,Num);
for i=1:Num
    transmat = [rotationMatrixZYX(-oxts{i+k}(4),-oxts{i+k}(5),-oxts{i+k}(6)) zeros(3,1); zeros(1,3) 1];
    liveFrames{i} = transformFrameTransMat(lidarData{i}(:,1:3), transmat);
    [easting, northing] = latlonToSweref991330(oxts{i+k}(1),oxts{i+k}(2));
    offset{i} = [easting; northing; oxts{i+k}(3)+0.5] - staticZeroGeo;
    angle = deg2rad(2);
    rotmat = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
    liveFrames{i} = transformFrameTransMat(liveFrames{i}(:,1:3), [rotmat offset{i};zeros(1,3) 1]);
end

%% plot live frames vs. static geo map
figure;
for j=1:Num
    test = liveFrames{j};
    pcshow(test)
    hold on
    for i = 1:length(walls)
        plotCubes(walls{i}(1:3)',walls{i}(4),walls{i}(5),walls{i}(6),walls{i}(7:9),-0.5,walls{i}(3)-2)
    end
    for i = 1:length(poles)
        plotCylinder(poles{i}(2:4)',poles{i}(1),poles{i}(5)+poles{i}(4),poles{i}(6))
    end
    hold off
    pause(0.3)
end

%% ground removal
cleanedFrames = cell(1,Num);
for j = 1:Num
    tic
    cleanedFrames{j} = gridGroundRemoval(liveFrames{j}, 200, 0.3);
    toc
    size(cleanedFrames{j},1)
    fprintf('compression: %6.2f\n', size(cleanedFrames{j},1)/size(liveFrames{j},1));
end

%% remove static points
mWalls = cell2mat(walls);
for i=1:Num
    i
    tic
    % remove points that are within the walls
    % sort walls by closest to current position first
    currentPos = offset{1};
    currentPos = repmat(currentPos',16,1);
    [trash idx] = sort([sum(abs(mWalls(:,1:3)-currentPos),2)], 'ascend');
    mWalls = mWalls(idx,:);
    % remove all points that are within walls
    for j = 1:size(mWalls,1)
        cleanedFrames{i} = removePointsInsideCube(mWalls(j,:), cleanedFrames{i});
    end
    toc
end

%% plot cleaned frames vs. static geo map
figure;
for j=1:Num
    pcshow(cleanedFrames{j})
%     hold on
%     for i = 1:length(walls)
%         plotCubes(walls{i}(1:3)',walls{i}(4),walls{i}(5),walls{i}(6),walls{i}(7:9),-0.5,walls{i}(3)-2)
%     end
%     for i = 1:length(poles)
%         plotCylinder(poles{i}(2:4)',poles{i}(1),poles{i}(5)+poles{i}(4),poles{i}(6))
%     end
%     hold off
    pause(0.3)
end