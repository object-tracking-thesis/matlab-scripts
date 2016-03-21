%% clear
clear all;
clc;

%% load lidar and oxts data for live
path1 = '~/Downloads/thesis/share/pcap_scenarios/';
path2 = 'car/';
path3 = 'oxts/';
%path4 = 'pcd/2000to2200/';
path4 = 'pcd/1550to1650/';

oxts = loadOxtsDir(strcat(path1,path2,path3));
lidarData = loadLidarDir(strcat(path1,path2,path4));
Num = length(lidarData);

%% load lidar data for static
path1 = 'static/wsp_square_1.txt';
path2 = 'static/wsp_square_2.txt';
path3 = 'static/wsp_square_3.txt';
path4 = 'static/wsp_square_4.txt';
staticFrame1 = xyzirgb2matlab(path1);
staticFrame2 = xyzirgb2matlab(path2);
staticFrame3 = xyzirgb2matlab(path3);
staticFrame4 = xyzirgb2matlab(path4);
staticFrame = [staticFrame1; staticFrame2; staticFrame3; staticFrame4];
clear staticFrame1;
clear staticFrame2;
clear staticFrame3;
clear staticFrame4;
staticFrame(:,4:7) = [];

staticZero = [106349.981; 6406149.9800000004; 133.126];

%% downsample the static map
staticCloud = pointCloud(staticFrame(:,1:3));
gridStep = 0.5;
staticCloudDownsampled = pcdownsample(staticCloud, 'gridAverage', gridStep);

%% rotate and translate the live frames according to their gps data
k = 1550; %frame offset to the gps data
liveFrames = cell(1,Num);
offset = cell(1,Num);
for i=1:Num
    transmat = [rotationMatrixZYX(-oxts{i+k}(4),-oxts{i+k}(5),-oxts{i+k}(6)) zeros(3,1); zeros(1,3) 1];
    liveFrames{i} = transformFrameTransMat(lidarData{i}(:,1:3), transmat);
    [easting, northing] = latlonToSweref991330(oxts{i+k}(1),oxts{i+k}(2));
    offset{i} = [easting; northing; oxts{i+k}(3)+0.5] - staticZero;
    angle = deg2rad(2);
    rotmat = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
    liveFrames{i} = transformFrameTransMat(liveFrames{i}(:,1:3), [rotmat offset{i};zeros(1,3) 1]);
end

%% show live frames
figure
for frame=1:100
    %plot3(0,0,0,'.','MarkerSize',20);
    %hold on
    live = pointCloud(liveFrames{frame}(:,1:3));
    pcshow(live)
    %hold off
    str = sprintf('Frame: %d', [frame]);
    title(str);
    zoom(2)
    pause(0.2)
end

%% plot difference between transformed and untransformed liveframe
figure
for frame=1:Num
    plot3(0,0,0,'.','MarkerSize',20);
    hold on
    live = pointCloud(lidarData{frame}(:,1:3));
    liveTrans = pointCloud(liveFrames{frame}(:,1:3));
    pcshowpair(live, liveTrans)
    hold off
    zoom(2)
    pause(0.5)
end

%% plot difference between static and live map
figure
for frame=1:Num
    plot3(0,0,0,'.','MarkerSize',20);
    hold on
    liveCloud = pointCloud(liveFrames{frame}(:,1:3));
    pcshowpair(staticCloudDownsampled, liveCloud)
    hold off
    zoom(4)
    %set(gca, 'CameraPosition', [18.5823 -566.5555   63.7898]);
    %set(gca, 'CameraViewAngle', 2.7243);
    pause(0.1)
end

%% knnsearch to remove the static map
cleanedFrames = cell(1,Num);
scalingz = 1.5;
scalingxy = 1;
off = 56;
staticFrame(:,3) = staticFrame(:,3).*(scalingz);
staticFrame(:,1:2) = staticFrame(:,1:2).*(scalingxy);
for i = 1:Num
   i
   tic
   limit = 1;
   test = liveFrames{i};
   test(:,3) = test(:,3).*(scalingz);
   test(:,1:2) = test(:,1:2).*(scalingxy);
   cleanedFrames{i} = subPC(staticFrame(:,1:3),test(:,1:3),limit);
   cleanedFrames{i}(:,3) = cleanedFrames{i}(:,3).*(1/(scalingz));
   cleanedFrames{i}(:,1:2) = cleanedFrames{i}(:,1:2).*(1/(scalingxy));
   toc
end
staticFrame(:,3) = staticFrame(:,3).*(1/scalingz);
staticFrame(:,1:2) = staticFrame(:,1:2).*(1/scalingxy);

%% plot difference between cleaned and original liveframes
figure
for i=1:Num
    orig = pointCloud(liveFrames{i}(:,1:3));
    clean = pointCloud(cleanedFrames{i}(:,1:3));
    subplot(2,2,1)
    pcshow(orig)
    zoom(3)
    subplot(2,2,2)
    pcshow(clean)
    zoom(3)
    %pcshowpair(orig, clean)
    pause(0.5)
end